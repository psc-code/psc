----------------------------------------------------------
-- SCRIPT FOR 2D 5M TIME-STEPPING THROUGH runTimeStep   --
--                                                      --
-- ALL RELEVANT PARAMETERS ARE OBTAINED FROM THE SCRIPT -- 
--   FOR INITIALIZATION.                                --
-- THIS SCRIPT SHOULD BE UNIVERSAL FOR DIFFERENT RUNS.  --
----------------------------------------------------------

------------------------------------------
-- LOAD LOCAL DOMAIN PARAMS FROM LIBMRC --
------------------------------------------
ggcm_mhd, rank, mx, my, mz, lx, ly, lz, hx, hy, hz, lxg, lyg, lzg, hxg, hyg, hzg, common_script = ...

-----------------------------------------------
-- LOAD PARAMS SHARED WITH INIT              --
-----------------------------------------------
-- Time-stepping needs gasGamma, lightSpeed, --
--   epsilon0, elc/mgnErrorSpeedFactor,      --
--   elcCharge/Mass, ionCharge/Mass, cfl     --
-----------------------------------------------
showlog = true
showlocallog = false
dofile(common_script)

elcIdx = 0
ionIdx = elcIdx + nr_moments
emIdx = nr_fluids*nr_moments

-----------------------------------------
-- HANDLE 1,2,3-D 5M/10M AUTOMATICALLY --
-----------------------------------------
if (nr_dims == 1) then
   myFluidUpdater = Updater.WavePropagation1D
   if (nr_moments == 5) then
      mySourceUpdater = Updater.ImplicitFiveMomentSrc1D
   elseif (nr_moments == 10) then
      mySourceUpdater = Updater.ImplicitTenMomentSrc1D
   end
elseif (nr_dims == 2) then
   myFluidUpdater = Updater.WavePropagation2D
   if (nr_moments == 5) then
      mySourceUpdater = Updater.ImplicitFiveMomentSrc2D
   elseif (nr_moments == 10) then
      mySourceUpdater = Updater.ImplicitTenMomentSrc2D
   end
elseif (nr_dims == 3) then
   myFluidUpdater = Updater.WavePropagation3D
   if (nr_moments == 5) then
      mySourceUpdater = Updater.ImplicitFiveMomentSrc3D
   elseif (nr_moments == 10) then
      mySourceUpdater = Updater.ImplicitTenMomentSrc3D
   end
end

-------------------------------------
-- TIME-STEPPING FOR LIMBRC-GKEYLL --
-------------------------------------
grid = createGrid()
-- create sol. needed by dimensional-splitting algorithms
-- e.g. qX is sol. after sweeping in x
-- q/qNew are the initial/final sol.
q = createData(grid)
qX = createData(grid)
if (nr_dims == 2 or nr_dims == 3) then qY = createData(grid) end
if (nr_dims == 3) then qZ = createData(grid) end
-- qNew as aliases
if (nr_dims == 1) then qNew = qX end
if (nr_dims == 2) then qNew = qY end
if (nr_dims == 3) then qNew = qZ end

-- duplicates in case that the 
-- time step should be re-taken
qDup = createData(grid)
qNewDup = createData(grid)

-- aliases for i/o of solvers
elc = q:alias(elcIdx, elcIdx + nr_moments)
ion = q:alias(ionIdx, ionIdx + nr_moments)
emf = q:alias(emIdx, emIdx + nr_em_comps)
elcX = qX:alias(elcIdx, elcIdx + nr_moments)
ionX = qX:alias(ionIdx, ionIdx + nr_moments)
emfX = qX:alias(emIdx, emIdx + nr_em_comps)
if (nr_dims == 2 or nr_dims == 3) then
   elcY = qY:alias(elcIdx, elcIdx + nr_moments)
   ionY = qY:alias(ionIdx, ionIdx + nr_moments)
   emfY = qY:alias(emIdx, emIdx + nr_em_comps)
end
if (nr_dims == 3) then
   elcZ = qZ:alias(elcIdx, elcIdx + nr_moments)
   ionZ = qZ:alias(ionIdx, ionIdx + nr_moments)
   emfZ = qZ:alias(emIdx, emIdx + nr_em_comps)
end
elcNew = qNew:alias(elcIdx, elcIdx + nr_moments)
ionNew = qNew:alias(ionIdx, ionIdx + nr_moments)
emfNew = qNew:alias(emIdx, emIdx + nr_em_comps)

------------------------
-- Boundary Condition --
------------------------
-- TODO: reuse temp_mrc_fld
temp_mrc_fld = ggcm_get_3d_fld(ggcm_mhd, 18)
temp_cptr = mrc_fld_get_arr(temp_mrc_fld)

function applyBc(fld, tCurr, myDt, ggcm_mhd)
 fld:copy_to_cptr(temp_cptr)
 ggcm_fill_ghosts(ggcm_mhd, temp_mrc_fld, tCurr)
 fld:copy_from_cptr(temp_cptr)
end

function finalize()
 ggcm_put_3d_fld(ggcm_mhd, temp_mrc_fld)
end

--------------------------
-- HYPERBOLIC EQUATIONS --
--------------------------
-- Euler equations with regular, high-
-- resolution fluxes
elcEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
ionEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
-- Euler equations with low-resolution, 
-- positivity-preserving fluxes
elcLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",   
}
ionLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
}
-- Maxwell equations
maxEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

-------------------------------------------
-- EQUATION SOLVERS ALONG EACH DIMENSION --
-------------------------------------------
function createSlvr(equation, direction, limiter)
   return myFluidUpdater {
      onGrid = grid,
      equation = equation,
      -- one of no-limiter, min-mod, superbee, 
      -- van-leer, monotonized-centered, beam-warming
      limiter = limiter,
      cfl = cfl,
      cflm = 1.1*cfl,
      updateDirections = {direction} -- directions to update
   }
end

-- dimensional-splitting solvers for regular equations
elcSlvrDir0 = createSlvr(elcEqn, 0, "van-leer")
ionSlvrDir0 = createSlvr(ionEqn, 0, "van-leer")
maxSlvrDir0 = createSlvr(maxEqn, 0, "van-leer")
elcSlvrDir0:setIn( {elc} )
elcSlvrDir0:setOut( {elcX} )
ionSlvrDir0:setIn( {ion} )
ionSlvrDir0:setOut( {ionX} )
maxSlvrDir0:setIn( {emf} )
maxSlvrDir0:setOut( {emfX} )

if (nr_dims == 2 or nr_dims == 3) then
   elcSlvrDir1 = createSlvr(elcEqn, 1, "van-leer")
   ionSlvrDir1 = createSlvr(ionEqn, 1, "van-leer")
   maxSlvrDir1 = createSlvr(maxEqn, 1, "van-leer")
   elcSlvrDir1:setIn( {elcX} )
   elcSlvrDir1:setOut( {elcY} )
   ionSlvrDir1:setIn( {ionX} )
   ionSlvrDir1:setOut( {ionY} )
   maxSlvrDir1:setIn( {emfX} )
   maxSlvrDir1:setOut( {emfY} )
end

if (nr_dims == 3) then
   elcSlvrDir2 = createSlvr(elcEqn, 2, "van-leer")
   ionSlvrDir2 = createSlvr(ionEqn, 2, "van-leer")
   maxSlvrDir2 = createSlvr(maxEqn, 2, "van-leer")
   elcSlvrDir2:setIn( {elcY} )
   elcSlvrDir2:setOut( {elcZ} )
   ionSlvrDir2:setIn( {ionY} )
   ionSlvrDir2:setOut( {ionZ} )
   maxSlvrDir2:setIn( {emfY} )
   maxSlvrDir2:setOut( {emfZ} )
end

-- dimensional-splitting solvers for positivity-preserving equations
elcLaxSlvrDir0 = createSlvr(elcEqn, 0, "zero")
ionLaxSlvrDir0 = createSlvr(ionEqn, 0, "zero")
maxLaxSlvrDir0 = createSlvr(maxEqn, 0, "zero")
elcLaxSlvrDir0:setIn( {elc} )
elcLaxSlvrDir0:setOut( {elcX} )
ionLaxSlvrDir0:setIn( {ion} )
ionLaxSlvrDir0:setOut( {ionX} )
maxLaxSlvrDir0:setIn( {emf} )
maxLaxSlvrDir0:setOut( {emfX} )

if (nr_dims == 2 or nr_dims == 3) then
   elcLaxSlvrDir1 = createSlvr(elcEqn, 1, "zero")
   ionLaxSlvrDir1 = createSlvr(ionEqn, 1, "zero")
   maxLaxSlvrDir1 = createSlvr(maxEqn, 1, "zero")
   elcLaxSlvrDir1:setIn( {elcX} )
   elcLaxSlvrDir1:setOut( {elcY} )
   ionLaxSlvrDir1:setIn( {ionX} )
   ionLaxSlvrDir1:setOut( {ionY} )
   maxLaxSlvrDir1:setIn( {emfX} )
   maxLaxSlvrDir1:setOut( {emfY} )
end

if (nr_dims == 3) then
   elcLaxSlvrDir2 = createSlvr(elcEqn, 2, "zero")
   ionLaxSlvrDir2 = createSlvr(ionEqn, 2, "zero")
   maxLaxSlvrDir2 = createSlvr(maxEqn, 2, "zero")
   elcLaxSlvrDir2:setIn( {elcY} )
   elcLaxSlvrDir2:setOut( {elcZ} )
   ionLaxSlvrDir2:setIn( {ionY} )
   ionLaxSlvrDir2:setOut( {ionZ} )
   maxLaxSlvrDir2:setIn( {emfY} )
   maxLaxSlvrDir2:setOut( {emfZ} )
end

---------------------
-- SOURCE UPDATERS --
---------------------
sourceSlvr = mySourceUpdater {
   onGrid = grid,
   numFluids = nr_fluids,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
   epsilon0 = epsilon0,
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

function updateSource(elcIn, ionIn, emIn, tCurr, t)
   sourceSlvr:setOut( {elcIn, ionIn, emIn} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
end

-----------------------------------------------------
-- COMPLETE EQUATION SOLUTION WITHOUT SOURCE TERMS --
-----------------------------------------------------
elcSlvr = {elcSlvrDir0, elcSlvrDir1, elcSlvrDir2}
ionSlvr = {ionSlvrDir0, ionSlvrDir1, ionSlvrDir2}
maxSlvr = {maxSlvrDir0, maxSlvrDir1, maxSlvrDir2}
elcOut = {elcX, elcY, elcZ}
ionOut = {ionX, ionY, ionZ}
qOut = {qX, qY, qZ}
function updateFluidsAndField(tCurr, t, ggcm_mhd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = False

   for dir = 1,nr_dims do
      for i,slvr in ipairs({elcSlvr[dir], ionSlvr[dir], maxSlvr[dir]}) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(t)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end
      myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
      myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)
 
      if ((elcEqn:checkInvariantDomain(elcOut[dir]) == false)
         or (ionEqn:checkInvariantDomain(ionOut[dir]) == false)) then
         useLaxSolver = true
      end
      useLaxSolver = ggcm_mhd_reduce_boolean(ggcm_mhd, useLaxSolver, false)
 
      if ((myStatus == false) or (useLaxSolver == true)) then
         return myStatus, myDtSuggested, useLaxSolver
      end
 
      applyBc(qOut[dir], tCurr, t-tCurr, ggcm_mhd)
   end

   return myStatus, myDtSuggested, useLaxSolver
end

elcLaxSlvr = {elcLaxSlvrDir0, elcLaxSlvrDir1, elcLaxSlvrDir2}
ionLaxSlvr = {ionLaxSlvrDir0, ionLaxSlvrDir1, ionLaxSlvrDir2}
maxLaxSlvr = {maxLaxSlvrDir0, maxLaxSlvrDir1, maxLaxSlvrDir2}
function updateFluidsAndFieldLax(tCurr, t, ggcm_mhd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)

   for dir = 1,nr_dims do
      for i,slvr in ipairs({elcLaxSlvr[dir], ionLaxSlvr[dir], maxLaxSlvr[dir]}) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(t)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end
      myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
      myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)

      if (myStatus == false) then
         return myStatus, myDtSuggested
      end

      applyBc(qOut[dir], tCurr, t-tCurr, ggcm_mhd)
   end

   return myStatus, myDtSuggested
end

-----------------------
-- COMPLETE TIMESTEP --
-----------------------
-- timestepping with regular fluxes
function solveTwoFluidSystem(tCurr, t, ggcm_mhd)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr, ggcm_mhd)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t, ggcm_mhd)

   -- update source terms
   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr, ggcm_mhd)

   return status, dtSuggested,useLaxSolver
end

-- timestepping with positivity-preserving fluxes
function solveTwoFluidLaxSystem(tCurr, t, ggcm_mhd)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elc, ion, emf, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr, ggcm_mhd)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t, ggcm_mhd)

   -- update source terms
   updateSource(elcNew, ionNew, emfNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr, ggcm_mhd)

   return status, dtSuggested
end

----------------------------------
-- GRAND TIME-STEPPING FUNCTION --
----------------------------------
function runTimeStep(myDt, tCurr, step, cptr, ggcm_mhd)
   q:copy_from_cptr(cptr)
   useLaxSolver = false

   local loop = 1
   while true do
      -- copy q and qNew in case we need to take this step again
      qDup:copy(q)
      qNewDup:copy(qNew)

      -- advance fluids and fields
      if (useLaxSolver) then
         status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt, ggcm_mhd)
         useLaxSolver = false
      else
         status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt, ggcm_mhd)
      end

      if (useLaxSolver == true) then
         mprint (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr+myDt))
         q:copy(qDup)
         qNew:copy(qNewDup)
      elseif (status == false) then
         mprint (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
         myDt = dtSuggested
         q:copy(qDup)
         qNew:copy(qNewDup)
      else
         local hasNan = qNew:hasNan()
         hasNan = ggcm_mhd_reduce_boolean(ggcm_mhd, hasNan, false)
         if (hasNan and rank == 0) then
            error (string.format(" ** NaN occured at %g! Stopping simulation", tCurr))
         end

         -- copy updated solution back
         q:copy(qNew)
         tCurr = tCurr + myDt
         myDt = dtSuggested
         break
      end 
      mprint("loop", loop)
   end -- end of retry loop

   q:copy_to_cptr(cptr)

   return myDt
end

-------------------------------------------------
-- SCRIPT FOR 1D/2D/3D MULTIPLE-SPECIES 5M/10M --
-- TIME-STEPPING THROUGH runTimeStep           --
-------------------------------------------------

------------------------------------------
-- LOAD LOCAL DOMAIN PARAMS FROM LIBMRC --
------------------------------------------
ggcm_mhd, rank, mx, my, mz, lx, ly, lz, hx, hy, hz, lxg, lyg, lzg, hxg, hyg, hzg, gasGamma, common_script = ...

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

-----------------------------------------
-- HANDLE 1,2,3-D 5M/10M AUTOMATICALLY --
-----------------------------------------
if (nr_moments == 5) then
   fluidEqnBase = HyperEquation.Euler
elseif (nr_moments == 10) then
   fluidEqnBase = HyperEquation.TenMoment
end

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

-------------------------
-- SETUP GRID AND DATA --
-------------------------
-- create grid for local domain
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
qNewDup = createData(grid) --FIXME needed?

-- begining index of the sth species
-- s ranges from 0 to nr_fluids-1!!!
function fluidIdx(s)
   return s * nr_moments
end
-- begining index of the EM field (and correction potentials)
emfIdx = nr_fluids * nr_moments

-- aliases for i/o of solvers
function getAliases(myQ)
   local myFluids = {}
   for s=0,nr_fluids-1 do
      myFluids[s] = myQ:alias(fluidIdx(s), fluidIdx(s) + nr_moments)
   end
   local myEmf = myQ:alias(emfIdx, emfIdx + 8)
   return myFluids,myEmf
end

fluids,emf = getAliases(q)
fluidsX,emfX = getAliases(qX)
if (nr_dims == 2 or nr_dims == 3) then
   fluidsY,emfY = getAliases(qY)
end
if (nr_dims == 3) then
   fluidsZ,emfZ = getAliases(qZ)
end
fluidsNew,emfNew = getAliases(qNew)

------------------------
-- Boundary Condition --
------------------------
-- create mrc_fld to be cached and reused
temp_mrc_fld = ggcm_get_3d_fld(ggcm_mhd, nr_moments * nr_fluids + 8)
temp_cptr = mrc_fld_get_arr(temp_mrc_fld)

function applyBc(myQ, tCurr, myDt)
   -- copy current solution to the cache field
   myQ:copy_to_cptr(temp_cptr)
   -- call ggcm's boundary condition to work on the cache field
   ggcm_fill_ghosts(ggcm_mhd, temp_mrc_fld, tCurr)
   -- copy the updated cache field back
   myQ:copy_from_cptr(temp_cptr)
end

-- free the space of the temporary field
function finalize()
 ggcm_put_3d_fld(ggcm_mhd, temp_mrc_fld)
end

--------------------------
-- HYPERBOLIC EQUATIONS --
--------------------------
-- Euler equations with regular, high-resolution fluxes
fluidEqns = {}
-- Euler equations with low-resolution, positivity-preserving fluxes
fluidLaxEqns = {}
for s=0,nr_fluids-1 do
   fluidEqns[s] = fluidEqnBase {
      gasGamma = gasGamma,
   }
   fluidLaxEqns[s] = fluidEqnBase {
      gasGamma = gasGamma,
      numericalFlux = "lax",
   }
end
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

function createFluidEMSlvrs(eqnType, dir, limiter)
   if (eqnType == "regular") then
      myFluidEqns = fluidEqns
   elseif (eqnType == "lax") then
      myFluidEqns = fluidLaxEqns
   end
   if (dir == 0) then
      myFluidsIn = fluids
      myFluidsOut = fluidsX
      myEmfIn = emf
      myEmfOut = emfX
   elseif (dir == 1) then
      myFluidsIn = fluidsX
      myFluidsOut = fluidsY
      myEmfIn = emfX
      myEmfOut = emfY
   elseif (dir == 2) then
      myFluidsIn = fluidsY
      myFluidsOut = fluidsZ
      myEmfIn = emfY
      myEmfOut = emfZ
   end
   local myFluidSlvrsDir = {}
   for s=0,nr_fluids-1 do
      myFluidSlvrsDir[s] = createSlvr(myFluidEqns[s], dir, limiter)
      myFluidSlvrsDir[s]:setIn({myFluidsIn[s]})
      myFluidSlvrsDir[s]:setOut({myFluidsOut[s]})
   end
   local myMaxSlvrDir = createSlvr(maxEqn, dir, limiter)
   myMaxSlvrDir:setIn( {myEmfIn} )
   myMaxSlvrDir:setOut( {myEmfOut} )
   return myFluidSlvrsDir,myMaxSlvrDir
end

-- dimensional-splitting solvers for regular equations
fluidSlvrsDir0,maxSlvrDir0 = createFluidEMSlvrs("regular", 0, "van-leer")
if (nr_dims == 2 or nr_dims == 3) then
   fluidSlvrsDir1,maxSlvrDir1 = createFluidEMSlvrs("regular", 1, "van-leer")
end
if (nr_dims == 3) then
   fluidSlvrsDir2,maxSlvrDir2 = createFluidEMSlvrs("regular", 2, "van-leer")
end

-- dimensional-splitting solvers for positivity-preserving equations
fluidLaxSlvrsDir0,maxLaxSlvrDir0 = createFluidEMSlvrs("lax", 0, "zero")
if (nr_dims == 2 or nr_dims == 3) then
   fluidLaxSlvrsDir1,maxLaxSlvrDir1 = createFluidEMSlvrs("lax", 1, "zero")
end
if (nr_dims == 3) then
   fluidLaxSlvrsDir2,maxLaxSlvrDir2 = createFluidEMSlvrs("lax", 2, "zero")
end

---------------------
-- SOURCE UPDATERS --
---------------------
sourceSlvr = mySourceUpdater {
   onGrid = grid,
   numFluids = nr_fluids,
   charge = charge,
   mass = mass,
   epsilon0 = epsilon0,
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

-- output of source updater using the initial sol.
sourceOut = {}
-- output of source updater using the updated sol.
sourceOutNew = {}
for s=0,nr_fluids-1 do
   -- FIXME index matters?
   sourceOut[s+1] = fluids[s]
   sourceOutNew[s+1] = fluidsNew[s]
end
sourceOut[nr_fluids+1] =  emf
sourceOutNew[nr_fluids+1] = emfNew

function updateSource(mySourceOut, tCurr, t)
   sourceSlvr:setOut(mySourceOut)
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
end

-----------------------------------------------------
-- COMPLETE EQUATION SOLUTION WITHOUT SOURCE TERMS --
-----------------------------------------------------
-- output of fluid solvers in each direction
-- for each species
fluidsOut = {}
fluidsOut[0] = fluidsX
if (nr_dims ==2 or nr_dims == 3) then fluidsOut[1] = fluidsY end
if (nr_dims == 3) then fluidsOut[2] = fluidsZ end
-- output of all solvers in each direction
qOut = {}; qOut[0] = qX; qOut[1] = qY; qOut[2] = qZ

-- helper function to get slvrs in along a specific direction
function getSlvrsDir(mySlvrs, slvrType, dir)
   if (slvrType == "regular") then
      if (dir == 0) then
         fluidSlvrsDir = fluidSlvrsDir0
         maxSlvrDir = maxSlvrDir0
      elseif (dir == 1) then
         fluidSlvrsDir = fluidSlvrsDir1
         maxSlvrDir = maxSlvrDir1
      elseif (dir == 2) then
         fluidSlvrsDir = fluidSlvrsDir2
         maxSlvrDir = maxSlvrDir2
      end
   elseif (slvrType == "lax") then
      if (dir == 0) then
         fluidSlvrsDir = fluidLaxSlvrsDir0
         maxSlvrDir = maxLaxSlvrDir0
      elseif (dir == 1) then
         fluidSlvrsDir = fluidLaxSlvrsDir1
         maxSlvrDir = maxLaxSlvrDir1
      elseif (dir == 2) then
         fluidSlvrsDir = fluidLaxSlvrsDir2
         maxSlvrDir = maxLaxSlvrDir2
      end
   end
   slvrsDir = {}
   for s=0,nr_fluids-1 do
      -- note that lua iterates from [1] even when
      -- the table has value at [0]
      slvrsDir[s+1] = fluidSlvrsDir[s]
   end
   slvrsDir[nr_fluids+1] = maxSlvrDir
   mySlvrs[dir] = slvrsDir
end

-- regular solvers to be executed along each direction,
-- including both fluid and Maxwell equations
slvrs = {}
getSlvrsDir(slvrs, "regular", 0)
if (nr_dims == 2 or nr_dims == 3) then getSlvrsDir(slvrs, "regular", 1) end
if (nr_dims == 3) then getSlvrsDir(slvrs, "regular", 2) end

function updateFluidsAndField(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = false

   for dir = 0,nr_dims-1 do
      for i,slvr in ipairs(slvrs[dir]) do
         slvr:setCurrTime(tCurr)
         local status, dtSuggested = slvr:advance(t)
         myStatus = status and myStatus
         myDtSuggested = math.min(myDtSuggested, dtSuggested)
      end
      myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
      myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)

      for s=0,nr_fluids-1 do
         if (fluidEqns[s]:checkInvariantDomain((fluidsOut[dir])[s]) == false) then
            useLaxSolver = true -- TODO breaka for loop
         end
      end
      useLaxSolver = ggcm_mhd_reduce_boolean(ggcm_mhd, useLaxSolver, false)
 
      if ((myStatus == false) or (useLaxSolver == true)) then
         return myStatus, myDtSuggested, useLaxSolver
      end
 
      applyBc(qOut[dir], tCurr, t-tCurr)
   end

   return myStatus, myDtSuggested, useLaxSolver
end

-- positive-preserving solvers to be executed along each direction
laxSlvrs = {}
getSlvrsDir(laxSlvrs, "lax", 0)
if (nr_dims == 2 or nr_dims == 3) then getSlvrsDir(laxSlvrs, "lax", 1) end
if (nr_dims == 3) then getSlvrsDir(laxSlvrs, "lax", 2) end

function updateFluidsAndFieldLax(tCurr, t)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)

   for dir = 0,nr_dims-1 do
      for i,slvr in ipairs(laxSlvrs[dir]) do
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

      applyBc(qOut[dir], tCurr, t-tCurr)
   end

   return myStatus, myDtSuggested
end

-----------------------
-- COMPLETE TIMESTEP --
-----------------------
-- timestepping with regular fluxes
function solveTwoFluidSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(sourceOut, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(sourceOutNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested,useLaxSolver
end

-- timestepping with positivity-preserving fluxes
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(sourceOut, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms
   updateSource(sourceOutNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr)

   return status, dtSuggested
end

----------------------------------
-- GRAND TIME-STEPPING FUNCTION --
----------------------------------
function runTimeStep(myDt, tCurr, step, cptr)
   q:copy_from_cptr(cptr)
   useLaxSolver = false

   --FIXME qDup and qNewDup?

   while true do
      -- copy q and qNew in case we need to take this step again
      qDup:copy(q)
      qNewDup:copy(qNew)

      -- advance fluids and fields
      if (useLaxSolver) then
         status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt)
         if (status == true) then
            useLaxSolver = false
         end
      else
         status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt)
      end

      if (status == false) then
         mprint (string.format(" ** Time step %g too large! Will retake with dt %g; useLaxSolver = %s", myDt, dtSuggested, tostring(useLaxSolver)))
         myDt = dtSuggested
         q:copy(qDup)
         qNew:copy(qNewDup)
      elseif (useLaxSolver == true) then
         mprint (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr))
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
   end -- end of retry loop

   q:copy_to_cptr(cptr)

   return myDt
end

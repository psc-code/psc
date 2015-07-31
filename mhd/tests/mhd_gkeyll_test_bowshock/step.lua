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
rank, mx, my, mz, lx, ly, lz, hx, hy, hz, lxg, lyg, lzg, hxg, hyg, hzg = ...

-----------------------------------------------
-- LOAD PARAMS SHARED WITH INIT              --
-----------------------------------------------
-- Time-stepping needs gasGamma, lightSpeed, --
--   epsilon0, elc/mgnErrorSpeedFactor,      --
--   elcCharge/Mass, ionCharge/Mass, cfl     --
-----------------------------------------------
showlog = true
showlocallog = true
dofile("common.lua")

function mprint(s)
   if (rank == 0) then
      print(s)
   end
end

-------------------------------------
-- TIME-STEPPING FOR LIMBRC-GKEYLL --
-------------------------------------
grid = createGrid()
q = createFields(grid)
qX = createFields(grid)
qNew = createFields(grid)
qDup = createFields(grid)
qNewDup = createFields(grid)

-- TODO: use nr_comps to calculate indices
elcFluid = q:alias(0, 5)
ionFluid = q:alias(5, 10)
emField = q:alias(10, 18)

elcFluidX = qX:alias(0, 5)
ionFluidX = qX:alias(5, 10)
emFieldX = qX:alias(10, 18)

elcFluidNew = qNew:alias(0, 5)
ionFluidNew = qNew:alias(5, 10)
emFieldNew = qNew:alias(10, 18)


------------------------
-- Boundary Condition --
------------------------
-- function to apply boundary conditions to specified field
-- FIXME universal?
function applyBc(fld, tCurr, myDt, ggcm_mhd)
 temp_mrc_fld = ggcm_get_3d_fld(ggcm_mhd, 18)
 temp_cptr = mrc_fld_get_arr(temp_mrc_fld)
 fld:copy_to_cptr(temp_cptr)
 ggcm_fill_ghosts(ggcm_mhd, temp_mrc_fld, tCurr)
 fld:copy_from_cptr(temp_cptr)
 ggcm_put_3d_fld(ggcm, temp_mrc_fld)
end

----------------------
-- EQUATION SOLVERS --
----------------------
-- regular Euler equations
elcEulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
ionEulerEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
}
-- (Lax equations are used to fix negative pressure/density)
elcEulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",   
}
ionEulerLaxEqn = HyperEquation.Euler {
   gasGamma = gasGamma,
   numericalFlux = "lax",
}
maxwellEqn = HyperEquation.PhMaxwell {
   lightSpeed = lightSpeed,
   elcErrorSpeedFactor = elcErrorSpeedFactor,
   mgnErrorSpeedFactor = mgnErrorSpeedFactor
}

-- ds solvers for regular Euler equations along X
elcFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   -- one of no-limiter, min-mod, superbee, 
   -- van-leer, monotonized-centered, beam-warming
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0} -- directions to update
}
ionFluidSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for regular Euler equations along Y
elcFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionFluidSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
maxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- ds solvers for Lax Euler equations along X
elcLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
ionLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}
maxLaxSlvrDir0 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {0}
}

-- ds solvers for Lax Euler equations along Y
elcLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = elcEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
ionLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = ionEulerLaxEqn,
   limiter = "zero",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}
maxLaxSlvrDir1 = Updater.WavePropagation2D {
   onGrid = grid,
   equation = maxwellEqn,
   limiter = "van-leer",
   cfl = cfl,
   cflm = 1.1*cfl,
   updateDirections = {1}
}

-- updater for source terms
sourceSlvr = Updater.ImplicitFiveMomentSrc2D {
   onGrid = grid,
   numFluids = 2,
   charge = {elcCharge, ionCharge},
   mass = {elcMass, ionMass},
   epsilon0 = epsilon0,
   -- linear solver to use: one of partialPivLu or colPivHouseholderQr
   linearSolver = "partialPivLu",
   hasStaticField = false,
}

-- function to update source terms
function updateSource(elcIn, ionIn, emIn, tCurr, t)
   sourceSlvr:setOut( {elcIn, ionIn, emIn} )
   sourceSlvr:setCurrTime(tCurr)
   sourceSlvr:advance(t)
   -- momentum drag relaxation
end

-- function to update the fluid and field using dimensional splitting
function updateFluidsAndField(tCurr, t, ggcm_mhd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   local useLaxSolver = False
   -- X-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir0, ionFluidSlvrDir0, maxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end
   myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
   myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)

   if ((elcEulerEqn:checkInvariantDomain(elcFluidX) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidX) == false)) then
      useLaxSolver = true
   end
   useLaxSolver = ggcm_mhd_reduce_boolean(ggcm_mhd, useLaxSolver, false)

   if ((myStatus == false) or (useLaxSolver == true)) then
      return myStatus, myDtSuggested, useLaxSolver
   end

   -- apply BCs to intermediate update after X sweep
   applyBc(qX, tCurr, t-tCurr, ggcm_mhd)

   -- Y-direction updates
   for i,slvr in ipairs({elcFluidSlvrDir1, ionFluidSlvrDir1, maxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end
   myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
   myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)

   if ((elcEulerEqn:checkInvariantDomain(elcFluidNew) == false)
    or (ionEulerEqn:checkInvariantDomain(ionFluidNew) == false)) then
       useLaxSolver = true
   end
   useLaxSolver = ggcm_mhd_reduce_boolean(ggcm_mhd, useLaxSolver, false)

   return myStatus, myDtSuggested, useLaxSolver
end

-- function to take one time-step with Euler solver
function solveTwoFluidSystem(tCurr, t, ggcm_mhd)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr, ggcm_mhd)

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t, ggcm_mhd)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr, ggcm_mhd)

   return status, dtSuggested,useLaxSolver
end

-- function to update the fluid and field using dimensional splitting Lax scheme
function updateFluidsAndFieldLax(tCurr, t, ggcm_mhd)
   local myStatus = true
   local myDtSuggested = 1e3*math.abs(t-tCurr)
   for i,slvr in ipairs({elcLaxSlvrDir0, ionLaxSlvrDir0, maxLaxSlvrDir0}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end
   myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
   myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)

   applyBc(qX, tCurr, t-tCurr, ggcm_mhd)

   -- Y-direction updates
   for i,slvr in ipairs({elcLaxSlvrDir1, ionLaxSlvrDir1, maxLaxSlvrDir1}) do
      slvr:setCurrTime(tCurr)
      local status, dtSuggested = slvr:advance(t)
      myStatus = status and myStatus
      myDtSuggested = math.min(myDtSuggested, dtSuggested)
   end
   myStatus = ggcm_mhd_reduce_boolean(ggcm_mhd, myStatus, true)
   myDtSuggested = ggcm_mhd_reduce_double_min(ggcm_mhd, myDtSuggested)

   return myStatus, myDtSuggested
end

-- function to take one time-step with Lax Euler solver
function solveTwoFluidLaxSystem(tCurr, t, ggcm_mhd)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(elcFluid, ionFluid, emField, tCurr, tCurr+dthalf)
   applyBc(q, tCurr, t-tCurr, ggcm_mhd)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t, ggcm_mhd)

   -- update source terms
   updateSource(elcFluidNew, ionFluidNew, emFieldNew, tCurr, tCurr+dthalf)
   applyBc(qNew, tCurr, t-tCurr, ggcm_mhd)

   return status, dtSuggested
end

----------------------------
-- DIAGNOSIS AND DATA I/O --
----------------------------
-- write data to H5 files
function writeFields(frame, t)
   qNew:write( string.format("q_%d.h5", frame), t )
end

----------------------------
-- TIME-STEPPING FUNCTION --
----------------------------
function runTimeStep(myDt, tCurr, step, cptr, ggcm_mhd)
   -- print(rank, "start", "myDt", myDt, "tCurr", tCurr)

   q:copy_from_cptr(cptr)

   while true do
      -- copy q and qNew in case we need to take this step again
      qDup:copy(q)
      qNewDup:copy(qNew)

      -- advance fluids and fields
      if (useLaxSolver) then
        -- call Lax solver if positivity violated
        -- mprint (string.format(" Taking step %5d at time %6g with dt %g (using Lax solvers)", step, tCurr, myDt))
        status, dtSuggested = solveTwoFluidLaxSystem(tCurr, tCurr+myDt, ggcm_mhd)
        useLaxSolver = false
      else
        -- mprint (string.format(" Taking step %5d at time %6g with dt %g", step, tCurr, myDt))
        status, dtSuggested, useLaxSolver = solveTwoFluidSystem(tCurr, tCurr+myDt, ggcm_mhd)
      end

      if (useLaxSolver == true) then
        -- negative density/pressure occured
        mprint (string.format(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr+myDt))
        q:copy(qDup)
        qNew:copy(qNewDup)
      elseif (status == false) then
        -- time-step too large
        mprint (string.format(" ** Time step %g too large! Will retake with dt %g", myDt, dtSuggested))
        myDt = dtSuggested
        qNew:copy(qNewDup)
        q:copy(qDup)
      else
        -- check if a nan occured
        if (qNew:hasNan()) then
           mprint (string.format(" ** NaN occured at %g! Stopping simulation", tCurr))
  	   assert(0)
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

-- set input/output arrays for various solvers
elcFluidSlvrDir0:setIn( {elcFluid} )
elcFluidSlvrDir0:setOut( {elcFluidX} )
ionFluidSlvrDir0:setIn( {ionFluid} )
ionFluidSlvrDir0:setOut( {ionFluidX} )
maxSlvrDir0:setIn( {emField} )
maxSlvrDir0:setOut( {emFieldX} )

elcFluidSlvrDir1:setIn( {elcFluidX} )
elcFluidSlvrDir1:setOut( {elcFluidNew} )
ionFluidSlvrDir1:setIn( {ionFluidX} )
ionFluidSlvrDir1:setOut( {ionFluidNew} )
maxSlvrDir1:setIn( {emFieldX} )
maxSlvrDir1:setOut( {emFieldNew} )

elcLaxSlvrDir0:setIn( {elcFluid} )
elcLaxSlvrDir0:setOut( {elcFluidX} )
ionLaxSlvrDir0:setIn( {ionFluid} )
ionLaxSlvrDir0:setOut( {ionFluidX} )
maxLaxSlvrDir0:setIn( {emField} )
maxLaxSlvrDir0:setOut( {emFieldX} )

elcLaxSlvrDir1:setIn( {elcFluidX} )
elcLaxSlvrDir1:setOut( {elcFluidNew} )
ionLaxSlvrDir1:setIn( {ionFluidX} )
ionLaxSlvrDir1:setOut( {ionFluidNew} )
maxLaxSlvrDir1:setIn( {emFieldX} )
maxLaxSlvrDir1:setOut( {emFieldNew} )


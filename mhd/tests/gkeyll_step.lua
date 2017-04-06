
-------------------------------------------------
-- SCRIPT FOR 1D/2D/3D MULTIPLE-SPECIES 5M/10M --
-- TIME-STEPPING THROUGH runTimeStep           --
-------------------------------------------------

-----------------------------------------------
-- Time-stepping needs gasGamma, lightSpeed, --
--   epsilon0, elc/mgnErrorSpeedFactor,      --
--   elcCharge/Mass, ionCharge/Mass, cfl     --
-----------------------------------------------

----------------------
-- HELPER FUNCTIONS --
----------------------
function log(...)
   if (rank == 0) then print(string.format(...)) end
end
function HERE()
   local info = debug.getinfo(2)
   local str = string.format("HERE: %d %s: %s",
       info.currentline, info.source, tostring(info.name))
   if (rank == 0) then print(str) end
end

-- print("qFld", qFld, "qFldX", qFldX, "qFldY", qFldY, "qFldZ", qFldZ)
-- print("ymask", ymask, "b0", b0)

hasStairSteppedBoundary = tostring(ymask) ~= "userdata: (nil)"
hasStaticField = tostring(b0) ~= "userdata: (nil)"

mu0 = 1.0
epsilon0 =1.0/mu0/(lightSpeed^2)
elcErrorSpeedFactor = 0.0
mgnErrorSpeedFactor = 1.0

nr_comps = nr_fluids * nr_moments + 8
nr_ghosts = 2
if (mz == 1) then
   if (my == 1) then
      nr_dims = 1
   else
      nr_dims = 2
   end
else
   nr_dims = 3
end

charge = {getCArray(charge_array, nr_fluids)}
mass = {getCArray(mass_array, nr_fluids)}

log("==================================================== GKEYLL")
log("%-19s | %s", "has_ymask", tostring(hasStairSteppedBoundary))
log("%-19s | %s", "has_b0", tostring(hasStaticField))

-------------------
-- GRID AND DATA --
-------------------
if nonuniform then

   if (nr_dims == 1) then
      myGrid = Grid.NonUniformRectCart1D
      lower = {0}
      upper = {1}
      cells = {mx}
      vertices = {crdx}
      myDataStruct = DataStruct.Field1D
   elseif (nr_dims == 2) then
      myGrid = Grid.NonUniformRectCart2D
      lower = {0, 0}
      upper = {1, 1}
      cells = {mx, my}
      vertices = {crdx, crdy}
      myDataStruct = DataStruct.Field2D
   elseif (nr_dims == 3) then
      myGrid = Grid.NonUniformRectCart3D
      lower = {0, 0, 0}
      upper = {1, 1, 1}
      cells = {mx, my, mz}
      vertices = {crdx, crdy, crdz}
      myDataStruct = DataStruct.Field3D
   end

   function createGrid()
      return myGrid {
         lower = lower,
         upper = upper,
         cells = cells,
         vertices = vertices,
      }
   end

else -- uniform

   if (nr_dims == 1) then
      myGrid = Grid.RectCart1D
      lower = {lx}
      upper = {hx}
      cells = {mx}
      myDataStruct = DataStruct.Field1D
   elseif (nr_dims == 2) then
      myGrid = Grid.RectCart2D
      lower = {lx, ly}
      upper = {hx, hy}
      cells = {mx, my}
      myDataStruct = DataStruct.Field2D
   elseif (nr_dims == 3) then
      myGrid = Grid.RectCart3D
      lower = {lx, ly, lz}
      upper = {hx, hy, hz}
      cells = {mx, my, mz}
      myDataStruct = DataStruct.Field3D
   end

   function createGrid()
      return myGrid {
         lower = lower,
         upper = upper,
         cells = cells,
      }
   end

end

function createData(numComponents, cptr)
   if cptr then
      return myDataStruct {
         onGrid = grid,
         numComponents = numComponents,
         ghost = {nr_ghosts, nr_ghosts},
         rawPointer = cptr,
      }
   else
      return myDataStruct {
         onGrid = grid,
         numComponents = numComponents,
         ghost = {nr_ghosts, nr_ghosts},
      }
   end
end

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

if (nr_dims == 1) then
   FieldFunction = Updater.FieldFunction1D
elseif (nr_dims == 2) then
   FieldFunction = Updater.FieldFunction2D
elseif (nr_dims == 3) then
   FieldFunction = Updater.FieldFunction3D
end

-------------------------
-- SETUP GRID AND DATA --
-------------------------
-- create grid for local domain
grid = createGrid()
-- create sol. needed by dimensional-splitting algorithms
-- e.g. qX is sol. after sweeping in x
-- q/qNew are the initial/final sol.
q = createData(nr_comps, mrc_fld_get_arr(qFld))
qX = createData(nr_comps, mrc_fld_get_arr(qFldX))
if (nr_dims == 2 or nr_dims == 3) then qY = createData(nr_comps, mrc_fld_get_arr(qFldY)) end
if (nr_dims == 3) then qZ = createData(nr_comps, mrc_fld_get_arr(qFldZ)) end
-- qNew as aliases
if (nr_dims == 1) then qNew = qX end
if (nr_dims == 2) then qNew = qY end
if (nr_dims == 3) then qNew = qZ end

-- duplicates in case that the 
-- time step should be re-taken
qDup = createData(nr_comps)
qNewDup = createData(nr_comps) --FIXME needed?

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

qFlds = {}
if qFld then qFlds[q] = qFld end
if qFldX then qFlds[qX] = qFldX end
if qFldY then qFlds[qY] = qFldY end
if qFldZ then qFlds[qZ] = qFldZ end

------------------------
-- Boundary Condition --
------------------------
function applyBc(myQ, tCurr, myDt)
   -- call ggcm's boundary condition to work on the cache field
   ggcm_fill_ghosts(ggcm_mhd, qFlds[myQ], tCurr)
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
if (hasStairSteppedBoundary) then
   inOutField = createData(1)
end

-- FIXME: inOutField not correctly set yet
function createSlvr(equation, direction, limiter)
   return myFluidUpdater {
      onGrid = grid,
      equation = equation,
      -- one of no-limiter, min-mod, superbee, 
      -- van-leer, monotonized-centered, beam-warming
      limiter = limiter,
      cfl = cfl,
      cflm = 1.1*cfl,
      updateDirections = {direction}, -- directions to update
      hasStairSteppedBoundary = hasStairSteppedBoundary,
      inOutField = inOutField,
      zeroLimiterSsBnd = true,
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
if (hasStaticField) then
   staticField = createData(8)
   staticBField = staticField:alias(3,6)
end

sourceSlvr = mySourceUpdater {
   onGrid = grid,
   numFluids = nr_fluids,
   charge = charge,
   mass = mass,
   epsilon0 = epsilon0,
   linearSolver = "partialPivLu",
   hasStaticField = hasStaticField,
}
if (hasStaticField) then
   sourceSlvr:setIn({staticField})
end

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
qIn = {}; qIn[0] = q; qIn[1] = qX; qIn[2] = qY
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
      applyBc(qIn[dir], tCurr, t-tCurr)
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
      applyBc(qIn[dir], tCurr, t-tCurr)
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

   -- update fluids and fields
   local status, dtSuggested, useLaxSolver = updateFluidsAndField(tCurr, t)

   -- update source terms
   updateSource(sourceOutNew, tCurr, tCurr+dthalf)

   return status, dtSuggested,useLaxSolver
end

-- timestepping with positivity-preserving fluxes
function solveTwoFluidLaxSystem(tCurr, t)
   local dthalf = 0.5*(t-tCurr)

   -- update source terms
   updateSource(sourceOut, tCurr, tCurr+dthalf)

   -- update fluids and fields
   local status, dtSuggested = updateFluidsAndFieldLax(tCurr, t)

   -- update source terms
   updateSource(sourceOutNew, tCurr, tCurr+dthalf)

   return status, dtSuggested
end

copyStaticField = true
copyInOutField = true
ymask2InOutField = FieldFunction {
   onGrid = grid,
   inpComponents = {0},
   outComponents = {0},
   func = function(x,y,z,t,val)
      if val == 0 then return -1
      else return val end
   end,
}

----------------------------------
-- GRAND TIME-STEPPING FUNCTION --
----------------------------------
function runTimeStep(myDt, tCurr, step, cptr)
   -- FIXME: make sure ghost points on b0 and ymask are synced?
   if (copyStaticField) then
      if (hasStaticField) then
         staticBField:copy_from_cptr(mrc_fld_get_arr(b0))
         -- staticBField:write("staticBField.h5", tCurr)
      end
      copyStaticField = false
   end
   if (copyInOutField) then
      if (hasStairSteppedBoundary) then
         inOutField:copy_from_cptr(mrc_fld_get_arr(ymask))
         --[[ FIXME FieldFunction not working
         -- lua_rawgeti(*L, LUA_REGISTRYINDEX, fnRef);
         ymask2InOutField:setIn( {inOutField} )
         ymask2InOutField:setOut( {inOutField} )
         -- ggcm ymask value 0 to gkeyll inOutField value -1
         ymask2InOutField:advance(0)
         --]]
         local shift = createData(1)
         shift:set(function () return -.5 end)
         inOutField:accumulate(1., shift)
         -- inOutField:write("inOutField.h5", tCurr)
      end
      copyInOutField = false
   end

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
         log(" ** Time step %g too large! Will retake with dt %g; useLaxSolver = %s", myDt, dtSuggested, tostring(useLaxSolver))
         myDt = dtSuggested
         q:copy(qDup)
         qNew:copy(qNewDup)
      elseif (useLaxSolver == true) then
         log(" ** Negative pressure or density at %8g! Will retake step with Lax fluxes", tCurr)
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

   return myDt
end


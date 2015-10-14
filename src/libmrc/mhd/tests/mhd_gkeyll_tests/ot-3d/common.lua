-- PROBLEM: 2D FIVE-MOMENT, ORSZAG-TANG --

-----------------------------------------------
-- PARAMETERS AND FUNCTIONS COMMONLY USED BY --
-- INITIALIZATION AND TIME-STEPPING          --
-----------------------------------------------
-- Assumes knowledge of global domain sizes  --
-- i.e., lyg, lzg, hxg, hyg, hzg             --
-----------------------------------------------

skip_execute = ...

--------------------------
-- PARAMETERS TO LIBMRC --
--------------------------
nr_moments = 5
nr_fluids = 2
nr_em_comps = 8
nr_comps = nr_moments*nr_fluids + nr_em_comps
nr_ghosts = 2
nr_dims = 3

-- if we only need the parameters above and want
-- to skip executing the remaining codes, do not
-- specify skip_execute (nil) or set it false
if (skip_execute) then
   return
end

----------------------
-- HELPER FUNCTIONS --
----------------------
pi = math.pi
sqrt = math.sqrt
cos = math.cos
sin = math.sin
tan = math.tan

function mprint(s)
   if (rank == 0) then
      print(s)
   end
end

-----------------------------------------------
-- PHYSICAL PARAMETERS                       --
-----------------------------------------------
-- Time-stepping needs gasGamma, lightSpeed, --
--   epsilon0, elc/mgnErrorSpeedFactor,      --
--   elcCharge/Mass, ionCharge/Mass, cfl     --
-----------------------------------------------
lightSpeed = 20.0
mu0 = 1.0
epsilon0 =1.0/mu0/(lightSpeed^2)
elcErrorSpeedFactor = 0.0
mgnErrorSpeedFactor = 1.0

rho0 = 25. / 9.
pr0 = 5. / 3.
v0 = 1.
B0 = 1.

elcCharge = -1.0
ionCharge = 1.0
charge = {elcCharge, ionCharge}

ionElcMassRatio = 25.
ionElcPressRatio = 1.
ionInertiaLength0 = 0.25

rhoe0 = rho0 / (1.+ionElcMassRatio)
rhoi0 = rho0 - rhoe0
pre0 = pr0 / (1.+ionElcPressRatio)
pri0 = pr0 - pre0
ionMass = ionInertiaLength0 * sqrt(mu0*rhoi0*ionCharge*ionCharge)
elcMass = ionMass / ionElcMassRatio
mass = {elcMass, ionMass}

cfl = 0.9

-----------------------
-- INITIAL CONDITION --
-----------------------
function init(x,y,z)
   local Lx = hxg - lxg
   local Ly = hyg - lyg
   
   local vx = -v0*sin(2.*pi*y/Ly)
   local vy = v0*sin(2.*pi*x/Lx)

   local momxe = rhoe0*vx
   local momye = rhoe0*vy
   local ere = pre0/(gasGamma-1.) + 0.5*(momxe^2+momye^2)/rhoe0
   
   local momxi = rhoi0*vx
   local momyi = rhoi0*vy
   local eri = pri0/(gasGamma-1.) + 0.5*(momxi^2+momyi^2)/rhoi0

   local Bx = -B0*sin(2.*pi*y/Ly) 
   local By = B0*sin(4.*pi*x/Lx)
   local Ez = -vx*By + vy*Bx

   return rhoe0, momxe, momye, 0.0, ere, rhoi0, momxi, momyi, 0.0, eri, 0.0, 0.0, Ez, Bx, By, 0.0, 0.0, 0.0
end

------------------------------------
-- LOGS TO BE DISPLAYED ON SCREEN --
------------------------------------
if (showlog) then
   mprint(string.format("===================================================="))
   mprint(string.format("nr_fluids = %d  nr_moments = %d", nr_fluids, nr_moments))
   mprint(string.format("nr_comps = %d  nr_ghosts = %d nr_dims = %d", nr_comps, nr_ghosts, nr_dims))
   mprint(string.format("lightSpeed = %g  mu0 = %g  epsilon0 = %g", lightSpeed, mu0, epsilon0))
   mprint(string.format("elcErrorSpeedFactor = %g  mgnErrorSpeedFactor = %g", elcErrorSpeedFactor, mgnErrorSpeedFactor))
   mprint(string.format("ionMass = %g  elcMass = ionMass/%g",ionMass, ionMass/elcMass))
   mprint(string.format("ionCharge = %g  elcCharge = %g",ionCharge, elcCharge))
   mprint(string.format("cfl = %g", cfl))
   mprint(string.format("===================================================="))
end

if (showlocallog) then
   print(string.format("[%03d] dims = [%d,%d,%d] l = [%g,%g,%g] h = [%g,%g,%g]", rank,mx,my,mz,lx,ly,lz,hx,hy,hz))
end

------------------------
-- COMMON CODES       --
------------------------
-- 1,2,3-D COMPATIBLE --
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

function createData(grid)
   return myDataStruct {
      onGrid = grid,
      numComponents = nr_comps,
      ghost = {nr_ghosts, nr_ghosts},
   }
end

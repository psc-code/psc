-- PROBLEM: 2D FIVE-MOMENT, ORSZAG-TANG --

-----------------------------------------------
-- PARAMETERS AND FUNCTIONS COMMONLY USED BY --
-- INITIALIZATION AND TIME-STEPPING          --
-----------------------------------------------
-- Assumes knowledge of global domain sizes  --
-- i.e., lyg, lzg, hxg, hyg, hzg             --
-----------------------------------------------

--------------------------
-- PARAMETERS TO LIBMRC --
--------------------------
nr_moments = 5
nr_fluids = 2
nr_comps = nr_moments*nr_fluids + 8
nr_ghosts = 2
-- TODO: skip the remaining codes if not needed

pi = math.pi
sqrt = math.sqrt
cos = math.cos
sin = math.sin
tan = math.tan
-----------------------------------------------
-- PHYSICAL PARAMETERS                       --
-----------------------------------------------
-- Time-stepping needs gasGamma, lightSpeed, --
--   epsilon0, elc/mgnErrorSpeedFactor,      --
--   elcCharge/Mass, ionCharge/Mass, cfl     --
-----------------------------------------------
gasGamma = 5./3.

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
ionElcMassRatio = 25.
ionElcPressRatio = 1.
ionInertiaLength0 = 0.25

rhoe0 = rho0 / (1.+ionElcMassRatio)
rhoi0 = rho0 - rhoe0
pre0 = pr0 / (1.+ionElcPressRatio)
pri0 = pr0 - pre0
ionMass = ionInertiaLength0 * sqrt(mu0*rhoi0*ionCharge*ionCharge)
elcMass = ionMass / ionElcMassRatio

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
   local Bz = 0.0

   return rhoe0, momxe, momye, 0.0, ere, rhoi0, momxi, momyi, 0.0, eri, 0.0, 0.0, 0.0, Bx, By, Bz, 0.0, 0.0
end

------------------------------------
-- LOGS TO BE DISPLAYED ON SCREEN --
------------------------------------
if (showlog and rank==0) then
   print(string.format("===================================================="))
   print(string.format("nr_comps = %g  nr_ghosts = %g", nr_comps, nr_ghosts))
   print(string.format("lightSpeed = %g  mu0 = %g  epsilon0 = %g", lightSpeed, mu0, epsilon0))
   print(string.format("elcErrorSpeedFactor = %g  mgnErrorSpeedFactor = %g", elcErrorSpeedFactor, mgnErrorSpeedFactor))
   print(string.format("ionMass = %g  elcMass = ionMass/%g",ionMass, ionMass/elcMass))
   print(string.format("ionCharge = %g  elcCharge = %g",ionCharge, elcCharge))
   print(string.format("cfl = %g", cfl))
   print(string.format("===================================================="))
end

if (showlocallog) then
   print(string.format("[%03d] dims = [%d,%d,%d] l = [%g,%g,%g] h = [%g,%g,%g]", rank,mx,my,mz,lx,ly,lz,hx,hy,hz))
end

------------------
-- COMMON CODES --
------------------
function createGrid()
   return Grid.RectCart2D {
      lower = {lx, ly},
      upper = {hx, hy},
      cells = {mx, my},
   }
end

function createFields(grid)
   return DataStruct.Field2D {
      onGrid = grid,
      numComponents = nr_comps,
      ghost = {nr_ghosts, nr_ghosts},
   }
end


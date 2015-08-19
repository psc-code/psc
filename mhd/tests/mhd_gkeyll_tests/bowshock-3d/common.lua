-- PROBLEM: 2D FIVE-MOMENT, BOWSHOCK FORMATION --

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

-- NOTE: lua array indices start from 1, not 0
mass_ratios = {1./26., 25./26.} 
momentum_ratios = {1./26, 25./26.}
pressure_ratios = {0.5, 0.5}

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
lightSpeed = 10.0
mu0 = 1.0
epsilon0 =1.0/mu0/(lightSpeed^2)
elcErrorSpeedFactor = 0.0
mgnErrorSpeedFactor = 1.0

scale = 1.0
-- background
rho0 = 0.01*scale
pr0 = 0.0015
vx0 = 1.0
-- dense core
rhoCore = 100.0/scale
prCore = 0.0015
vxCore = 0.0
xCore = lxg + 0.25 * (hxg - lxg)
yCore = 0.0
zCore = 0.0
radCore = 0.0625

elcCharge = -1.0
ionCharge = 1.0
charge = {elcCharge, ionCharge}

ionElcMassRatio = mass_ratios[2]/mass_ratios[1]
ionElcPressRatio = pressure_ratios[2]/pressure_ratios[1] 
ionElcMomentumRatio = momentum_ratios[2]/momentum_ratios[1]
ionInertiaLength0 = 0.06

BxIn = 0.0
ByIn = 0.001
BzIn = 0.0

--- derived params
rhoe0 = rho0 / (ionElcMassRatio + 1.0)
rhoi0 = rho0 - rhoe0
momxe0 = rhoe0*vx0
momxi0 = rhoi0*vx0
pre0 = pr0 / (ionElcPressRatio + 1.0)
pri0 = pr0 - pre0
ere0 = pre0/(gasGamma-1.0) + 0.5*momxe0^2/rhoe0
eri0 = pri0/(gasGamma-1.0) + 0.5*momxi0^2/rhoi0

rhoeCore = rhoCore / (ionElcMassRatio + 1.0)
rhoiCore = rhoCore - rhoeCore
momxeCore = rhoeCore*vxCore
momxiCore = rhoiCore*vxCore
preCore = prCore / (ionElcPressRatio + 1.0)
priCore = prCore - preCore
ereCore = preCore/(gasGamma-1.0) + 0.5*momxeCore^2/rhoeCore
eriCore = priCore/(gasGamma-1.0) + 0.5*momxiCore^2/rhoiCore

ionMass = ionInertiaLength0 * sqrt(mu0*rhoi0*ionCharge*ionCharge)
elcMass = ionMass / ionElcMassRatio
mass = {elcMass, ionMass}

rhoeIn = rhoe0
momxeIn = momxe0
ereIn = ere0
rhoiIn = rhoi0
momxiIn = momxi0
eriIn = eri0

cfl = 0.9

-- diagnostic params
BIn = sqrt(BxIn^2+ByIn^2+BzIn^2)
pMagIn = BIn*BIn/2.0/mu0

vA0 = BIn/sqrt(mu0*rho0)
cs0 = sqrt(gasGamma*pr0/rho0)
plasmaBeta0 = pr0 / pMagIn

ionInertiaLengthCore = sqrt(ionMass*ionMass/mu0/rhoiCore/(ionCharge^2))
csCore = sqrt(gasGamma*prCore/rhoCore)

-----------------------
-- INITIAL CONDITION --
-----------------------
-- initial conditions
function linearChange(x, x1, val1, x2, val2)
   -- x1<x2
   local val = val1
   if (x > x2) then
      val = val2
   elseif (x > x1) then
      val = val1+(val2-val1)*(x-x1)/(x2-x1)
   end
   return val
end
function init(x,y,z)
   local Lx = hxg - lxg
   local Ly = hyg - lyg

   local rhoe = rhoe0
   local rhoi = rhoi0
   local vx = vx0
   local pre = pre0
   local pri = pri0

   rad = sqrt((x-xCore)^2+(y-yCore)^2+(z-zCore)^2)
   if (rad < radCore) then
      rhoe = rhoeCore
      rhoi = rhoiCore
      vx = vxCore
      pre = preCore
      pri = priCore
   end
   if ((x > xCore) and sqrt((y-yCore)^2+(z-zCore)^2) < radCore) then
      vx = vxCore
   end
   local momx = (rhoe + rhoi) * vx
   local momxe = momx / (1 + ionElcMomentumRatio)
   local momxi = momx - momxe
   local ere = pre/(gasGamma-1) + 0.5*momxe^2/rhoe
   local eri = pri/(gasGamma-1) + 0.5*momxi^2/rhoi
   local Bx = 0.0
   local By = 0.0
   local Bz = 0.0
   if (x < xCore-radCore) then
      Bx = linearChange(x, -Lx/2, BxIn, xCore-radCore, 0.0)
      By = linearChange(x, -Lx/2, ByIn, xCore-radCore, 0.0)
      Bz = linearChange(x, -Lx/2, BzIn, xCore-radCore, 0.0)
   end

   return rhoe, momxe, 0.0, 0.0, ere, rhoi, momxi, 0.0, 0.0, eri, 0.0, 0.0, 0.0, Bx, By, Bz, 0.0, 0.0
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

   Lx = hxg - lxg
   Ly = hyg - lyg
   mprint(string.format("Lx=%gdi0=%gdiCore", Lx/ionInertiaLength0, Lx/ionInertiaLengthCore))
   mprint(string.format("Ly=%gdi0=%gdiCore", Ly/ionInertiaLength0, Ly/ionInertiaLengthCore))
   mprint(string.format("                                "))
   mprint(string.format("Background, |(x,y,z) - (%g,%g,%g)| > %g:", xCore, yCore, zCore, radCore))
   mprint(string.format("  rho=%g  di=%g", rho0, ionInertiaLength0))
   mprint(string.format("  cs=%g=%gc  vA=%g=%gc", cs0, cs0/lightSpeed, vA0, vA0/lightSpeed))
   mprint(string.format("  vx=%g=%gc=%gcs=%gvA", vx0, vx0/lightSpeed, vx0/cs0, vx0/vA0))
   mprint(string.format("  plasmaBeta=%g", plasmaBeta0))
   mprint(string.format("  pr=%g", pr0))
   mprint(string.format("                                "))
   mprint(string.format("Core, |(x,y,z) - (%g,%g,%g)| < %g:", xCore, yCore, zCore, radCore))
   mprint(string.format("  rho=%g  di=%g", rhoCore, ionInertiaLengthCore))
   mprint(string.format("  cs=%g=%gc", csCore, csCore/lightSpeed))
   mprint(string.format("  vx=%g=%gc=%gcs", vxCore, vxCore/lightSpeed, vxCore/csCore))
   mprint(string.format("  pr=%g", prCore))
   mprint(string.format("                                "))
   mprint(string.format("Inflow"))
   mprint(string.format("  rhoeIn=%g  momxeIn=%g ereIn=%g", rhoeIn, momxeIn, ereIn))
   mprint(string.format("  rhoiIn=%g  momxiIn=%g eriIn=%g", rhoiIn, momxiIn, eriIn))
   mprint(string.format("  BxIn=%g ByIn=%g, BzIn=%g", BxIn, ByIn, BzIn))
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

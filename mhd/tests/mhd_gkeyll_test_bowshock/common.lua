-----------------------------------------------
-- PARAMETERS AND FUNCTIONS COMMONLY USED BY --
-- INITIALIZATION AND TIME-STEPPING          --
-----------------------------------------------

--------------------------
-- PARAMETERS TO LIBMRC --
--------------------------
nr_comps = 18
nr_ghosts = 2

Pi = Lucee.Pi
sqrt = math.sqrt
cos = math.cos
sin = math.sin
tan = math.tan

-----------------------------------------------
-- PHYSICAL PARAMETERS                       --
-- Time-stepping needs gasGamma, lightSpeed, --
--   epsilon0, elc/mgnErrorSpeedFactor,      --
--   elcCharge/Mass, ionCharge/Mass, cfl     --
-----------------------------------------------
gasGamma = 5./3.

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
xCore = -0.75
yCore = 0.0
radCore = 0.0625

elcCharge = -1.0
ionCharge = 1.0
ionElcMassRatio = 25.
ionElcPressRatio = 1.
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

   rad = sqrt((x-xCore)^2+(y-yCore)^2)
   if (rad < radCore) then
      rhoe = rhoeCore
      rhoi = rhoiCore
      vx = vxCore
      pre = preCore
      pri = priCore
   end
   if ((x > xCore) and sqrt((y-yCore)^2) < radCore) then
      vx = vxCore
   end
   local momxe = rhoe*vx
   local momxi = rhoi*vx
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
   Lx = hxg - lxg
   Ly = hyg - lyg
   Lucee.logInfo(string.format("====================================="))
   Lucee.logInfo(string.format("lightSpeed = %g  mu0 = %g  epsilon0 = %g", lightSpeed, mu0, epsilon0))
   Lucee.logInfo(string.format("elcErrorSpeedFactor = %g  mgnErrorSpeedFactor = %g", elcErrorSpeedFactor, mgnErrorSpeedFactor))
   Lucee.logInfo(string.format("ionMass = %g  elcMass = ionMass/%g",ionMass, ionMass/elcMass))
   Lucee.logInfo(string.format("ionCharge = %g  elcCharge = %g",ionCharge, elcCharge))
   Lucee.logInfo(string.format("cfl = %g", cfl))

Lucee.logInfo(string.format("Lx=%gdi0=%gdiCore", Lx/ionInertiaLength0, Lx/ionInertiaLengthCore))
Lucee.logInfo(string.format("Ly=%gdi0=%gdiCore", Ly/ionInertiaLength0, Ly/ionInertiaLengthCore))
   Lucee.logInfo(string.format("                                "))
   
   Lucee.logInfo(string.format("Background, |(x,y) - (%g,%g)| < %g:", xCore, yCore, radCore))
   Lucee.logInfo(string.format("  rho=%g  di=%g", rho0, ionInertiaLength0))
   Lucee.logInfo(string.format("  cs=%g=%gc  vA=%g=%gc", cs0, cs0/lightSpeed, vA0, vA0/lightSpeed))
   Lucee.logInfo(string.format("  vx=%g=%gc=%gcs=%gvA", vx0, vx0/lightSpeed, vx0/cs0, vx0/vA0))
   Lucee.logInfo(string.format("  plasmaBeta=%g", plasmaBeta0))
   Lucee.logInfo(string.format("  pr=%g", pr0))
   Lucee.logInfo(string.format("                                "))
   Lucee.logInfo(string.format("Core, |(x,y) - (%g,%g)| > %g:", xCore, yCore, radCore))
   Lucee.logInfo(string.format("  rho=%g  di=%g", rhoCore, ionInertiaLengthCore))
   Lucee.logInfo(string.format("  cs=%g=%gc", csCore, csCore/lightSpeed))
   Lucee.logInfo(string.format("  vx=%g=%gc=%gcs", vxCore, vxCore/lightSpeed, vxCore/csCore))
   Lucee.logInfo(string.format("  pr=%g", prCore))
   Lucee.logInfo(string.format("                                "))
   Lucee.logInfo(string.format("Inflow"))
   Lucee.logInfo(string.format("  rhoeIn=%g  momxeIn=%g ereIn=%g", rhoeIn, momxeIn, ereIn))
   Lucee.logInfo(string.format("  rhoiIn=%g  momxiIn=%g eriIn=%g", rhoiIn, momxiIn, eriIn))
   Lucee.logInfo(string.format("  BxIn=%g ByIn=%g, BzIn=%g", BxIn, ByIn, BzIn))
   Lucee.logInfo(string.format("====================================="))
end

if (showlocallog) then
   print("     rank = ", rank)
   print("nr_ghosts = ", nr_ghosts)
   print("nr_comps  = ", nr_comps)
   print("    dims  = ", mx, my, mz)
   print("       l  = ", lx, ly, lz)
   print("       h  = ", hx, hy, hz)
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


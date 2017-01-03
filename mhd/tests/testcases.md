
# Test Cases

* t001_ot: Ideal MHD Orszag-Tang vortex with ggcm_mhd_step c3_double

* t002_kh: Ideal MHD Kelvin Helmholtz test with ggcm_mhd_step c3_double

* t003_ic_mirdip: Test the "mirdip" initial condition on staggered B
  grid. This tests bnd_sphere and bnd_inoutflow to some extent, too.

* t004_ic_mirdip: Test the "mirdip" initial condition all cell
  centered grid. This tests bnd_sphere and bnd_inoutflow to some
  extent, too.

* t005_ic_mirdip: Test the "mirdip" initial condition on
  ggcm-staggered B grid. This tests bnd_sphere and bnd_inoutflow to
  some extent, too.

* t006_alfven_wing: Ganymede-based Alfven wing test without dipole
  field.

* t007_cpaw_mhd: MHD-based circular polarized Alfven wave test, 1d,
  Jimmy-MHD staggered scheme

* t008_cpaw_mhd: MHD-based circular polarized Alfven wave test, 1d,
  cell-centered scheme

* t009_ot_mhdcc: Ideal MHD Orzag-Tang vortex with mhdcc

* t010_ot_vlct: Ideal MHD Orzag-Tang vortex with VLCT

* t011_sod_c3: Sod shocktube (hydro) with c3_double

* t012_sod_mhdcc: Sod shocktube (hydro) with mhdcc

* t013_sod_vl: Sod shocktube (hydro) with VL

* t014_briowu_c3: Brio-Wu shocktube with c3_double

* t015_briowu_mhdcc: Brio-Wu shocktube with mhdcc

* t016_gem_c3: GEM challenge ideal MHD, c3_double

* t017_gem_mhdcc: GEM challenge ideal MHD, mhdcc

* t018_gem_mhdcc_res: GEM challenge resistive MHD, mhdcc

* t019_gem_mhdcc_hall: GEM challenge resistive Hall-MHD, mhdcc

* t020_gem_mhdcc_res_divbglm: GEM challenge resistive MHD, mhdcc, GLM divB cleaning

* t021_openggcm_step: Full OpenGGCM test (run the runme), hi-order,
  const resistivity. ref.sh / run.sh run step "fortran" / "c_float",
  respectively, and should produce identical results.

* t024_gem_mhdcc_bg: Same as t017_gem_mhdcc, but with a (fake) background field to test
  background option in mhdcc
















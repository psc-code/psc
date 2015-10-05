-------------------------------
-- SCRIPT FOR INITIALIZATION --
-------------------------------

------------------------------------------
-- LOAD LOCAL DOMAIN PARAMS FROM LIBMRC --
------------------------------------------
rank, mx, my, mz, lx, ly, lz, hx, hy, hz, lxg, lyg, lzg, hxg, hyg, hzg, gasGamma, common_script, cptr = ...

----------------------------------------
-- LOAD COMMON CODES SHARED WITH STEP --
----------------------------------------
showlog = false
showlocallog = false
dofile(common_script)

----------------------------
-- I.C. FOR LIBMRC-GKEYLL --
----------------------------
grid = createGrid()
q = createData(grid)
q:set(init)
q:copy_to_cptr(cptr)



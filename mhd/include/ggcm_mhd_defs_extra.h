
#ifndef GGCM_MHD_DEFS_EXTRA_H
#define GGCM_MHD_DEFS_EXTRA_H

// ----------------------------------------------------------------------
// openggcm fields
// These are dependent on what openggcm does in its Fortran code
// indices based on mhd-cliches.for

enum {
  _RR1,
  _RV1X,
  _RV1Y,
  _RV1Z,
  _UU1,
  _B1X,
  _B1Y,
  _B1Z,

  _RR2 = 8,
  _RV2X,
  _RV2Y,
  _RV2Z,
  _UU2,
  _B2X,
  _B2Y,
  _B2Z,

  _YMASK = 16,
  _ZMASK = 17,
  _CMSV = 18,

  _RR = 19,
  _PP,
  _VX,
  _VY,
  _VZ,
  _BX,
  _BY,
  _BZ,

  _TMP1 = 27,
  _TMP2,
  _TMP3,
  _TMP4,

  _FLX = 31,
  _FLY,
  _FLZ,

  _CX = 34,
  _CY,
  _CZ,
  _XTRA1 = 37,
  _XTRA2,
  _RESIS = 39,
  _CURRX = 40,
  _CURRY,
  _CURRZ,
  _RMASK = 43,

  _NR_FLDS,
};

#endif


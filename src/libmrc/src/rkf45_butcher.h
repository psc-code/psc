#ifndef RKF_BUTCHER_H
#define RKF_BUTCHER_H

// Fehlburg Butcher tableau (from wikipedia);

// C1 | A10
// C2 | A20 A21
// C3 | A30 A31 A32
// C4 | A40 A41 A42 A43
// C5 | A50 A51 A52 A53 A54
// -------------------------
//    | B5_0 B5_1 B5_2 B5_3 B5_4 B5_5
//    | B4_0 B4_1 B4_2 B4_3 B4_4 B4_5

#define C1 (1./4.)
#define C2 (3./8.)
#define C3 (12./13.)
#define C4 (1.)
#define C5 (1./2.)

#define A10 (1./4.)

#define A20 (3./32.)
#define A21 (9./32.)

#define A30 (1932./2197.)
#define A31 (-7200./2197.)
#define A32 (7296./2197.)

#define A40 (439./216.)
#define A41 (-8.)
#define A42 (3680./512.)
#define A43 (-845./4104.)


#define A50 (-8./27.)
#define A51 (2.)
#define A52 (-3544./2565.)
#define A53 (1859./4104.)
#define A54 (-11./40.)

//--------------------------------------------------------------------------
// Fifth order coefficients (b5_2 = 0)
#define B5_0 (16./135.)
#define B5_1 (0.0)
#define B5_2 (6656./12825.)
#define B5_3 (28561./56430.)
#define B5_4 (-9./50.)
#define B5_5 (2./55.)

// Fourth order coefficients (b4_2 = b4_6 = 0)
#define B4_0 (25./216.)
#define B4_1 (0.0)
#define B4_2 (1408./2565.)
#define B4_3 (2197./4104.)
#define B4_4 (-1./5.)
#define B4_5 (0.0)

#endif

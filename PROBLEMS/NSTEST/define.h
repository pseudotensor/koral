#define MYCOORDS MSPH1COORDS

#define MKS1R0 0.
#define RMIN 5.
#define RMAX 30.
#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(30.-MKS1R0))
//#define MINX RMIN
//#define MAXX RMAX
#define MINY 1.e-6
#define MAXY (M_PI-1.e-6)
#define MINZ -M_PI
#define MAXZ M_PI

#define TNX 80
#define TNY 80
#define TNZ 1

#define NTX 1
#define NTY 1
#define NTZ 1

#define INT_ORDER 1
#define TSTEPLIM .6

#define SPECIFIC_BC

#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define MAGNFIELD
#define GDETIN 1
//#define VECPOTGIVEN
#define BETANORMFACTOR 10.
#define B2RHORATIOMAX 20.


#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2

#define DTOUT1 1.

#define RHO_AMB 1.e-3
#define U_AMB 1.e-3
#define RHO_BLOB 1.e3

#define SILOOUTPUT 1
#define SILO2D_XZPLANE


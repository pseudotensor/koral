#define RADIATION
//#define RADSOURCEOFF
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.

#define myMKS1COORDS

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS KSCOORDS
#endif

#define OUTCOORDS KERRCOORDS
#define OUTVEL VELR
#define DTOUT1 10.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
//#define CGSOUTPUT
//#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF
#define PRINTGC_LEFT
//#define PRINTGC_RIGHT

#define NX 100
#define NY 1
#define NZ 1

#ifdef myMKS1COORDS
#define MKS1R0 .7
#define MINX 0.
#define MAXX 3.
#else
#define MINX (.8*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#endif

#define MINY 0.*Pi/4.
#define MAXY Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define SPECIFIC_BC

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5
#define UTPOT 1.
//#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define RK3_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define NODONUT 1
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

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
#define OUTVEL VEL4
#define DTOUT1 5.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
//#define CGSOUTPUT
#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT

#ifdef myMKS1COORDS
#define MKS1R0 0.
#define MINX (log(1.7-MKS1R0))
#define MAXX (log(50.-MKS1R0))
#define NX 100
#else
#define MINX (1.8*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#define NX 60
#endif

#define NY 60
#define NZ 1


#define MINY 0.*Pi/4.
#define MAXY Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define SPECIFIC_BC

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5
#define UTPOT .98
//#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e8))

#define INT_ORDER 1
#define RK3_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40
//#define RADIATION

//#define RADSOURCEOFF
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE
//#define ALLOW_EXPLICIT_RAD_SOURCE 1

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define myMKS1COORDS

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS KSCOORDS
#endif

#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define DTOUT1 1.e1
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
#define RADOUTPUTINZAMO

#ifdef myMKS1COORDS
#define MKS1R0 -2.
#define MINX (log(1.5-MKS1R0))
#define MAXX (log(50.-MKS1R0))
#define NX 64
#else
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 30.
#define NX 70
#endif

#define NY 48
#define NZ 1


#define MINY (0.005*Pi/4.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.

#define SPECIFIC_BC

#define GAMMA (4./3.)
#define KKK 9.e-4
#define ELL 4.5
#define UTPOT .975

#define RHOATMMIN  3.e-5
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e8))

#define INT_ORDER 2
#define TIMESTEPPING RK2

#define ENFORCEENTROPY
//#define FLUXMETHOD HLL_FLUX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

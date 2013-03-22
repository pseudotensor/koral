#define BHSPIN 0.

#define myMKS1COORDS

#define MYCOORDS MKS1COORDS
#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS KSCOORDS//SPHCOORDS//KERRCOORDS
#endif

#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define DTOUT1 1.e1
#define ALLSTEPSOUTPUT 0

#ifdef myMKS1COORDS
#define MKS1R0 -2.
#define MINX (log(1.9-MKS1R0))
#define MAXX (log(50.-MKS1R0))
#define NX 30
#else
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#define NX 80
#endif

#define NY 50
#define NZ 1

#define MINY 0.*Pi/4.
#define MAXY Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define SPECIFIC_BC
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT

#define GAMMA (4./3.)
#define KKK 0.03
#define ELL 4.5
#define UTPOT 1.
#define RHOATMMIN  1.e-4
#define UINTATMMIN  1.e-6

#define INT_ORDER 1
#define RK3_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define NODONUT 0
#define INFLOWING 0
#define NSTEPSTOP 50e10
#define NOUTSTOP 1000.

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40


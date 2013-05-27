#define RADIATION
#define RADSOURCEOFF
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define myMKS1COORDS

#define OMSCALE 0.

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS SPHCOORDS//KSCOORDS//SPHCOORDS//KERRCOORDS
#endif

//#define WIDENPRESSURE
#define WIDENPRESSUREPOWER 5.

#define IMAGETYPE "jpg"
#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
//#define CGSOUTPUT
#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF
//#define RADOUTPUTVELS
#define SKIPRADWAVESPEEDLIMITER
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT
//#define PRINTYGC_LEFT

#ifdef myMKS1COORDS
#define MKS1R0 0.
#define MINX (log(1.-MKS1R0))
#define MAXX (log(50.-MKS1R0))
#define NX 20
#else
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#define NX 80
#endif

#define NY 20
#define NZ 1


#define MINY (0.02*Pi/4.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5
#define UTPOT .98
//#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define RK2_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

#define RADIATION

//#define RADSOURCEOFF
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE
#define ALLOW_EXPLICIT_RAD_SOURCE 0

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define myMKS1COORDS

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS SCHWCOORDS
#endif

#define VELPRIM VELR
//#define BLOB

#define VISCOSITY
#define SIMPLEVISCOSITY
#define ALPHAVISC .1
#define ALPHATOTALPRESSURE
#define RMINVISC 2.

#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 5000.
#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT

#define PAR_D 1.e0
#define PAR_E 1.e-8

#ifdef myMKS1COORDS
#define MKS1R0 -2.
#define MINX (log(1.5-MKS1R0))
#define MAXX (log(16.-MKS1R0))//(log(16.-MKS1R0))
#define NX 72
#else
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 16.
#define NX 48
#endif

#define NY 32
#define NZ 1


#define MINY (0.005*Pi/2.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC
//#define PUREAXISOUTFLOW

#define GAMMA (4./3.)
#define ELL 4.5

 
#ifdef RADIATION

#define URIN (5.23e8/CCC)

#define KKK 7845//1.e3 //the higher KKK the hotter disk i.e. the lower density - the larger prad/pgas
#define UTPOT .983//.9715//.9715
#define RHOATMMIN  1.e-23
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e6))
#define DTOUT1 1.e1
#define CGSOUTPUT

#else
#define URIN 0.5

#define KKK 9.e-4//1.e-4
#define UTPOT .99
#define RHOATMMIN  3.e-1
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define DTOUT1 10.e0

#endif


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

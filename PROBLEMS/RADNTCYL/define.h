#define RADIATION
#define SKIPRADSOURCE
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

//#define RADVISCOSITY SHEARVISCOSITY
#define ALPHARADVISC 1.
#define ZEROTIMEINSHEAR

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 100.

#define NUMRADWAVESPEEDS

//#define myMCYL1COORDS

#define OMSCALE 1.

#ifdef myMKS1COORDS
#define MYCOORDS MCYL1COORDS
#else
#define MYCOORDS CYLCOORDS//KSCOORDS//SPHCOORDS//KERRCOORDS
#endif

#define IMAGETYPE "gif"
#define OUTCOORDS CYLCOORDS
#define OUTVEL VEL4
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
//#define CGSOUTPUT
#define RADOUTPUTINZAMO
//#define WAVESPEEDSATFACES
//#define FULLRADWAVESPEEDS
//#define RADOUTPUTINFF
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT
//#define PRINTYGC_LEFT

#ifdef myMCYL1COORDS
#define MKS1R0 0.
#define MINX (log(1.9-MKS1R0))
#define MAXX (log(50.-MKS1R0))
#define NX 30
#else
#define MINX (0.01)
#define MAXX 15.
#define NX 80
#endif

#define NY 40 
#define NZ 1

#define MINY 0
#define MAXY MAXX
#define MINZ 0.
#define MAXZ 1.
#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5
#define UTPOT .98
//#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define TIMESTEPPING RK2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

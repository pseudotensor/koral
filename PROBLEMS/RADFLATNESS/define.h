#define TMAX 1.e10
#define RADIATION
#define COMPTONIZATION
#define NCOMPTONIZATION

#define MYCOORDS MINKCOORDS
#define NX 5
#define NY 1
#define NZ 1
#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1

//#define RADOUTPUTINZAMO
#define RESOUTPUT_ASCII
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define ALLSTEPSOUTPUT 0
#define DTOUT1 1.0
#define SCAOUTPUT 1
#define TESTCOMPTINSCALARS4FLAT

#define GAMMA (4./3.)
#define NSTEPSTOP 1e200 //max n of steps
#define NOUTSTOP 20 //max n of steps

//#define EXPLICIT_RAD_SOURCE

#define MINX 0
#define MAXX 1
#define MINY 0
#define MAXY 1
#define MINZ 0
#define MAXZ 1

#define MASSCM 1.e11

#define LTEFACTOR .01
#define URFX 0.

#define KAPPA .1
#define KAPPAES 0.
#define PERIODIC_XBC
#define PERIODIC_YBC
#define PERIODIC_ZBC

#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-4
#define RADFORCEPREC 1.e-5

#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40

#define TMAX 1.e10
#define RADIATION

#define MYCOORDS MINKCOORDS
#define NX 120
#define NY 21
#define NZ 1
#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1
#define RK3STEPPING

#define RADOUTPUTINZAMO

//#define VELPRIMRAD VEL4

#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0
#define DTOUT1 .5
#define NOUTSTOP 100
#define GAMMA (4./3.)
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MINX -6
#define MAXX 3
#define MINY 0
#define MAXY 1.5
#define MINZ 0
#define MAXZ 1

#define MASSCM 1.

#define LTEFACTOR 1.
#define URFX 13.
#define URFY -13.

#define KAPPA 0.
#define KAPPAES 0.

#define SPECIFIC_BC

//#define PERIODIC_XBC
//#define PERIODIC_YBC
//#define PERIODIC_ZBC

#define PRINTGC_LEFT

#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-4
#define RADFORCEPREC 1.e-5

#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40

#define TMAX 1.e10
#define RADIATION

//#define LABRADFLUXES
//#define EDDINGTON_APR

#define MYCOORDS MINKCOORDS
#define NX 41
#define NY 41
#define NZ 1 
//#define ZSLICE 20

#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1
#define RK3STEPPING

#define RADOUTPUTINZAMO

#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0
#define DTOUT1 .05
#define GAMMA (4./3.)
#define NOUTSTOP 60
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MINX 0
#define MAXX 1
#define MINY 0
#define MAXY 1
#define MINZ 0
#define MAXZ 1

#define IXDOTMIN NX/2
#define IYDOTMIN NY/3
#define IZDOTMIN 0
#define IXDOTMAX NX/2
#define IYDOTMAX NY/3
#define IZDOTMAX 0
#define FXDOT 0.
#define FYDOT 0.9
#define FZDOT 0.

#define MASSCM 1.

#define LTEFACTOR 1.
#define URFX 0.

#define KAPPA 0.
#define KAPPAES 0.

#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC

#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-4
#define RADFORCEPREC 1.e-5

#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40

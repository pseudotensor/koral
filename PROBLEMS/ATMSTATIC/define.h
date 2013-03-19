#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-7
#define U2PPREC 1.e-6

#define TMAX 1.e10
#define GAMMA 1.4
//#define RADIATION
#define SCHWARZSCHILD
//#define KERR
#define BHSPIN .1
#define RK3STEPPING
//#define SPHERICAL
//#define MINKOWSKI

#define NX 400
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define INT_ORDER 1
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 1.e2

#define ALLSTEPSOUTPUT 1
#define NSTEPSTOP 20
#define VERBOSE0 0
//#define EXPLICIT_RAD_SOURCE

#define MASS 1.
#define MDOTEDD 2.23/16.*1e18*MASS //cm/s
#define LUMEDD 1.25e38*MASS //erg/s


#define MINX 1.e6
#define MAXX 2.e6

//#define LOGXGRID
//#define LOGPAR1 2.2
//#define LOGPAR2 2.

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#define MINZ -1.
#define MAXZ 1.

#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define EFLOOR 1.e-40
//#define CGSOUTPUT

#define PRINTGC_LEFT
#define PRINTGC_RIGHT
#define OUTPUTINZAMO

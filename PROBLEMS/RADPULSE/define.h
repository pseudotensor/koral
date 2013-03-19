#define U2PPREC 1.e-5
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define MYCOORDS MINKCOORDS

//#define LABRADFLUXES

//#define RK3STEPPING
#define RK2STEPPING
#define INT_ORDER 1
#define NX 101
#define NY 1
#define NZ 1
#define TSTEPLIM .8
#define INITTSTEPLIM (TSTEPLIM/10.)
#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0 // ?
#define MINMOD_THETA 1.
#define DTOUT1 100.
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 100
#define GAMMA (ldouble)(5./3.)
#define MINX -50.
#define MAXX 50.
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.

#define RHO_AMB 1.e0
#define T_AMB 1.e6

#define EFLOOR 1.e-50
#define BLOBP 100.
#define BLOBW 5.

#define RADIATION
//#define GASRADOFF

#define KAPPA 0.
#define KAPPAES 1.e3
//#define KAPPA 5.90799

//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT

//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE
//#define EXPLICIT_SUBSTEP_RAD_SOURCE
#define RADOUTPUTINZAMO
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

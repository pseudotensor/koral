#define U2PPREC 1.e-5
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define MYCOORDS MINKCOORDS

#define RK2STEPPING
#define INT_ORDER 1
#define NX 201
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 2.
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 100
#define GAMMA (ldouble)(5./3.)
#define MINX -200.
#define MAXX 200.
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.

#define RHO_AMB 1.e0
#define T_AMB 1.e6

#define EFLOOR 1.e-50
#define BLOBP 100.
#define BLOBW 1.

#define BLOBX1 (-20.)
#define BLOBX2 (20.)

#define RADIATION
#define MULTIRADFLUID

#define KAPPA 0.
#define KAPPAES 0.

#define EXPLICIT_RAD_SOURCE

//#define IMPLICIT_FF_RAD_SOURCE
//#define EXPLICIT_SUBSTEP_RAD_SOURCE
#define RADOUTPUTINZAMO

#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

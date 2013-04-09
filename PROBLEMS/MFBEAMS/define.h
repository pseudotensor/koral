#define U2PPREC 1.e-5
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define MYCOORDS MINKCOORDS

#define RADIATION
#define MULTIRADFLUID
#define REFLECT


#define RK2STEPPING
#define INT_ORDER 1
#define NX 201
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)

#define SPECIFIC_BC
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT

#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 .175
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 60
#define GAMMA (ldouble)(5./3.)
#define MINX -1.
#define MAXX 1.
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.

#define RHO_AMB 1.e0
#define T_AMB 1.e6

#define EEAMB 1.e-4
#define EEBEAM1 1.
#define EEBEAM2 1.
#define FRATIO1 .99
#define FRATIO2 .99


#define KAPPA 0.
#define KAPPAES 0.

#define EXPLICIT_RAD_SOURCE

//#define IMPLICIT_FF_RAD_SOURCE
//#define EXPLICIT_SUBSTEP_RAD_SOURCE
#define RADOUTPUTINZAMO

#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

#define TMAX 1.e10
//#define FLUXDISSIPATIONFULL
#define RADIATION

#define MYCOORDS MINKCOORDS
#define MINKOWSKI
#define NX 30
#define NY 1
#define NZ 30
#define YZXDUMP
#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1
#define RK3STEPPING

#define RADPUTINZAMO

#define INITTSTEPLIM (TSTEPLIM/10.)//for the 1st time step
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0
#define GAMMA (4./3.)
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE
#define IFBEAM 1
//#define GASRADOFF
//#define RADSOURCEOFF

#define MINX 0
#define MAXX 1
#define MINY 0
#define MAXY 1
#define MINZ 0
#define MAXZ 1

#define BEAML .4
#define BEAMR .6
#define DTOUT1 .1 //dt for basic output

#define RADBEAMFLAT_FRATIO 0.995
#define RADBEAMFLAT_ERAD 1.
#define RADBEAMFLAT_RHO 1.
#define RADBEAMFLAT_UU 0.1


#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40

#define KAPPA 0.
#define SPECIFIC_BC

#define PRINTGC_LEFT
#define PRINTGC_RIGHT

#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-4
#define RADFORCEPREC 1.e-5

#define TMAX 1.e10
//#define FLUXDISSIPATIONFULL
#define RADIATION
//#define myMKS1COORDS

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS KSCOORDS
#endif

#define MYCOORDS2 KERRCOORDS //suplementary, only for tetrads used when initializing problems - not yet used
#define OUTCOORDS KERRCOORDS //coordinates for output

#define NX 30
#define NY 1
#define NZ 30
#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1
#define RK3STEPPING

#define INITTSTEPLIM (TSTEPLIM/10.)//for the 1st time step
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0
#define GAMMA (1.4)
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define BEAMNO 1

#define IFBEAM 1
//#define GASRADOFF
//#define RADSOURCEOFF
#define FLATBACKGROUND


#if (BEAMNO==1)
#define MKS1R0 .5
#define MINX 2.6
#define MAXX 3.5
#define BEAML 2.9
#define BEAMR 3.1
#define DTOUT1 .1 //dt for basic output
#elif (BEAMNO==2)
#define MKS1R0 .5

#define MINX 5.5
#define MAXX 7.5
//#define MAXX 12.5
#define BEAML 5.8
#define BEAMR 6.2
#define DTOUT1 .4 //dt for basic output
#elif (BEAMNO==3)
#define MKS1R0 .5
#define BEAML 15.5
#define BEAMR 16.5
#ifdef myMKS1COORDS
#define MINX 2.5
#define MAXX 3.
#else
#define MINX 14.5
#define MAXX 20.5
#endif
#define DTOUT1 1. //dt for basic output
#define NLEFT 0.99
#elif (BEAMNO==4)
#define MINX 30
#define MAXX 50
#define BEAML 37
#define BEAMR 43
#define DTOUT1 .25 //dt for basic output
#endif

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.

#define MINZ 0.
#define MAXZ 2.*Pi/4.//8.*Pi/4.

//#define PLOTFULLPHI
#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40
#define TAMB 1e7
#define TLEFT 1e9

#ifndef NLEFT
#define NLEFT 0.995
#endif

#define RHOAMB 1.e0
#define KAPPA 0.
#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define YZXDUMP
#define RADOUTPUTINZAMO
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT
//#define BLOB
#define BLOBW .1
#define BLOBP 100000.
#define BLOBX 10.
#define BLOBZ Pi/20.
#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-4
#define RADFORCEPREC 1.e-5
#define PAR_D 1.
#define PAR_E 1e-4

//#define LOGXGRID

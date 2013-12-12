#define TMAX 1.e10
#define RADIATION

#define HDVISCOSITY NOVISCOSITY
#define RADVISCOSITY SHEARVISCOSITY
#define MAXRADVISCVEL 1./3.
#define ALPHARADVISC 1.
#define MAXRADVISCMFP (1./NX)
#define ZEROTIMEINSHEAR

//#define EDDINGTON_APR
//#define MULTIRADFLUID
#define MFSKEW 20.
#define MFMINVEL 1.e-3
#define MFREDISTRIBUTEMETHOD 3
#define MFWEDGESTYPE 1
#define NRF 4

//#define FORGETDOTS

//#define LABRADFLUXES

#define MYCOORDS MINKCOORDS
#define NX 40
#define NY 40
#define NZ 1
//#define ZSLICE 20

#define TSTEPLIM .5//kind of courant limiter
#define INT_ORDER 1
#define RK2STEPPING

#define RADOUTPUTINZAMO

#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0

#ifdef FORGETDOTS
#define DTOUT1 .005
#else
#define DTOUT1 .025
#endif

#define GAMMA (4./3.)
#define NOUTSTOP 100

#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MINX 0
#define MAXX 1
#define MINY 0
#define MAXY 1
#define MINZ 0
#define MAXZ 1


#define IXDOT1 (1./4.*NX)
#define IYDOT1 (NY/2)//10//10//20
#define IZDOT1 0
#define FXDOT1 0.
#define FYDOT1 0.
#define FZDOT1 0.

#define IXDOT2 (3./4.*NX)
#define IYDOT2 (NY/2)//30//20
#define IZDOT2 0
#define FXDOT2 0.
#define FYDOT2 0.
#define FZDOT2 0.

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

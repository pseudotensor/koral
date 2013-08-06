/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 11

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE
//#define SKIPRADWAVESPEEDLIMITER
//#define ALLOW_EXPLICIT_RAD_SOURCE 0

/************************************/
//hydro choices
/************************************/
#define ALLOWENTROPYU2P 1
#define DOFIXUPS 0

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2K1K2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 2.
//#define FLUXMETHOD HLL_FLUX
//#define WAVESPEEDSATFACES 
#define GDETIN 0

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
//#define ALPHATOTALPRESSURE
//#define RMINVISC 2.
#define RADVISCOSITY NOVISCOSITY
//#define TAUSUPPRESSPARAM 100. //the larger the less prad
//#define ALPHARADVISC 1.
//#define ENFORCERADWAVESPEEDS

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e6
#define EEUURATIOMIN 1.e-15
#define EEUURATIOMAX 1.e6
#define ERADLIMIT 1.e-50
#define RHOFLOOR 1.e-50
#define GAMMAMAXRAD 50.

/************************************/
//blackhole
/************************************/
#define BHSPIN 0.0
#define MASS 4.3e6

/************************************/
//coordinates / resolution
/************************************/

#define MYCOORDS MINKCOORDS
#define MINX -1.
#define MAXX 1. 
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1. 
#define NX 50 
#define NY 50
#define NZ 1
#define SPECIFIC_BC

/************************************/
//output
/************************************/
//#define SILOOUTPUT
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT
#define PRINTYGC_LEFT
#define PRINTYGC_RIGHT
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define DTOUT1 1.e1

/*
#ifdef VERTPLANEOUTPUT
#define ZSLICE (NZ*0.1)
#endif

#ifdef EQPLANEOUTPUT
#define YSLICE NY/2
#define PRINTZONEMORE
#endif
*/

/************************************/
//problem params
/************************************/
#define GAMMA (5./3.)
#define RHOZERO 1.
#define RHOBEAM 1.e2
#define TEMPZERO 1.e7
#define TEMPBEAM (1.e5*2.)
#define VELBULK .01
#define VELBEAM .02
#define BEAMX1 -.2
#define BEAMX2 .2

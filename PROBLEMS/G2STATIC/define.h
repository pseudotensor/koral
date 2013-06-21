/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 1

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE
//#define SKIPRADWAVESPEEDLIMITER
#define ALLOW_EXPLICIT_RAD_SOURCE 0

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
#define MINMOD_THETA 1.5
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
#define MASS 10.

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MKS1COORDS
#define MKS1R0 -8.
#define BHSPIN 0.0
#define MINX (log(1000.-MKS1R0))
#define MAXX (log(6000.-MKS1R0))
#define MINY (0.35*Pi)
#define MAXY (.65*Pi)
#define MINZ 0.
#define MAXZ (1.*Pi)
#define NX 30
#define NY 20
#define NZ 60
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
//#define RADOUTPUTINZAMO
//#define PRINTINSIDEBH
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define CGSOUTPUT
#define DTOUT1 5.e3
#define YSLICE NY/2
//#define PRINTZONEMORE

/************************************/
//common physics 
/************************************/
#define GAMMA (5./3.)
#define RHOATMMIN 1.
#define UINTATMMIN 1.e-4
#define TRACER
#define MINTRACE 0.01
#define CLMAG 100.
#define CLWIDTH 400.

//test
//#define MAGNFIELD

/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
#define RADIATION
#define RADCLOSURE VETCLOSURE
#define RADSTARTWITHM1INTENSITIES
#define RADCLOSURECOORDS MINKCOORDS
//#define NUMRADWAVESPEEDS
#define SKIPRADWAVESPEEDLIMITER
#define SKIPRADSOURCE
//#define NCOMPTONIZATION
//#define COMPTONIZATION

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define ERADLIMIT 1.e-50
#define RHOFLOOR 1.e-50
#define GAMMAMAXRAD 200.

/************************************/
//hydro choices
/************************************/
#define ALLOWENTROPYU2P 1
#define DOFIXUPS 0

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
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
//blackhole
/************************************/
#define BHSPIN 0.0
#define MASS 1.e10

/************************************/
//coordinates / resolution
/************************************/

#define MYCOORDS MINKCOORDS
#define MINX .1
#define MAXX (M_PI-0.1)
#define MINY .1
#define MAXY (M_PI-0.1)
#define MINZ -1.
#define MAXZ 1. 

#define TNX 200
#define TNY 200
#define TNZ 1

//# of tiles
#define NTX 1
#define NTY 1
#define NTZ 1

//#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define SHUFFLELOOPS 0

/************************************/
//output
/************************************/
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
//#define PRINTYGC_LEFT
//#define PRINTYGC_RIGHT

#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1e50
#define NOUTSTOP 5000
#define DTOUT1 0.1
#define RADOUTPUTINZAMO
#define SILOOUTPUT 1

/************************************/
//problem params
/************************************/
#define RRR
#define GAMMA (5./3.)
#define RHOZERO 0.01
#define BLOBMAG1 1.e4
#define BLOBMAG2 1.e4
#define TEMPZERO 1.e7
#define TEMPBLOB1 1.e8
#define TEMPBLOB2 1.e8
#define VELXBLOB1 0.
#define VELXBLOB2 -0.
#define VELYBLOB1 (VELXBLOB1*0.5)
#define VELYBLOB2 (VELXBLOB2*0.25)
#define XBLOB1 (MINX+0.35*(MAXX-MINX))
#define XBLOB2 (MINX+0.65*(MAXX-MINX))
#define YBLOB1 (MINY+0.35*(MAXY-MINY))
#define YBLOB2 (MINY+0.65*(MAXY-MINY))
#define SIZEBLOB1 (.025*(MAXX-MINX))
#define SIZEBLOB2 (.025*(MAXX-MINX))




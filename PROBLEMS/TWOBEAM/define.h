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
#define USEDUALNEIGHBOR
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES
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
#define GAMMAMAXRAD 20.

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
#define TSTEPLIM .3
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
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
#define RADCLOSURECOORDS MYCOORDS
//#define SILO2D_XZPLANE
#define MINX 0.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define TNX 50
#define TNY 50
#define TNZ 1

//# of tiles
#define NTX 4
#define NTY 1
#define NTZ 1

#define SPECIFIC_BC
//#define COPY_XBC
//#define COPY_YBC
//#define COPY_ZBC
#define SHUFFLELOOPS 0

/************************************/
//output
/************************************/
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
//#define PRINTYGC_LEFT
//#define PRINTYGC_RIGHT

#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e20
#define NOUTSTOP 5000
#define DTOUT1 0.02
#define RADOUTPUTINZAMO
#define SILOOUTPUT 1

/************************************/
//problem params
/************************************/





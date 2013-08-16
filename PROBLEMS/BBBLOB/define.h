/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 7

/************************************/
//radiation choices
/************************************/
#define RADIATION
//#define SKIPRADSOURCE

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
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define ERADLIMIT 1.e-50
#define RHOFLOOR 1.e-50
#define GAMMAMAXRAD 100.

/************************************/
//blackhole
/************************************/
#define BHSPIN 0.0
#define MASS 1.e10

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
#define NX 60 
#define NY 60
#define NZ 1

//#define SPECIFIC_BC
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC

/************************************/
//output
/************************************/
//#define SILOOUTPUT
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
//#define PRINTYGC_LEFT
//#define PRINTYGC_RIGHT

#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define DTOUT1 .5
#define RADOUTPUTINZAMO
#define SILOOUTPUT

/************************************/
//problem params
/************************************/
#define RRR
#define GAMMA (4./3.)
#define RHOZERO 0.01
#define BLOBMAG1 1.e4
#define BLOBMAG2 1.e4
#define TEMPZERO 1.e7
#define TEMPBLOB1 1.e7
#define TEMPBLOB2 1.e7
#define VELXBLOB1 0.1
#define VELXBLOB2 -0.1
#define VELYBLOB1 (VELXBLOB1*0.25)
#define VELYBLOB2 (VELXBLOB2*0.)
#define XBLOB1 -.5
#define XBLOB2 .5
#define SIZEBLOB1 .1
#define SIZEBLOB2 .1




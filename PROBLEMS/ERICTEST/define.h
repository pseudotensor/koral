/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
//#define RADIATION

/************************************/
//magnetic fields
/************************************/
//#define MAGNFIELD

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define GDETIN 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY NOVISCOSITY

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
//coordinates / resolution
/************************************/

#define MYCOORDS MINKCOORDS //KERRCOORDS //SPHCOORD //MKS1COORDS
#define MINX 0.
#define MAXX 1. 
#define MINY 0.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1. 
#define TNX 30 
#define TNY 30
#define TNZ 30
#define NTX 2
#define NTY 2
#define NTZ 2
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define SILOOUTPUT 1
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define DTOUT1 1.e0

/************************************/
//problem params
/************************************/
#define GAMMA (5./3.)
#define RHOAMB 1.
#define UUAMB 1.e-2
#define RHOCLOUD 1.e2
#define UUCLOUD 1.e-2
#define VELAMB  .1
#define SIGMACLOUD .1



/************************************/
//restart
/************************************/
//#define RESTART
//#define RESTARTNUM -1

/************************************/
//radiation
/************************************/
//#define RADIATION

/************************************/
//magnetic fields
/************************************/
//#define MAGNFIELD

/************************************/
//coordinates / resolution
/************************************/

#define MYCOORDS MINKCOORDS
#define MINX 0.
#define MAXX 1. 
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1. 
#define TNX 30 
#define TNY 30
#define TNZ 1
#define NTX 1
#define NTY 1
#define NTZ 1
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define SILOOUTPUT 1 //to silo file
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1.e10 //stop after this number of steps
#define NOUTSTOP 5000 //stop after this number of outputs
#define DTOUT1 1.e0 //res
#define DTOUT2 1.e50 //avg

/************************************/
//reconstruction / stepping
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

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
//physics
/************************************/
#define GAMMA (5./3.)

/************************************/
//problem params
/************************************/
#define RHOAMB 1.
#define UUAMB 1.e-3



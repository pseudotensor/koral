/************************************/
//restart
/************************************/
//#define RESTART 
//#define RESTARTNUM 39

#define BHDISK_PROBLEMTYPE

/************************************/
//radiation
/************************************/
//#define RADIATION

/************************************/
//magnetic fields
/************************************/
#define MAGNFIELD
#define VECPOTGIVEN  //we provide vector potential
#define MAXBETA .01 //target max pmag/pgas 
#define BETANORMFULL  //normalize everywhere


/************************************/
//coordinates / resolution
/************************************/

#define MYCOORDS KERRCOORDS
#define RMIN 5.
#define RMAX 20.
#define MINX RMIN
#define MAXX RMAX 
#define DTH .5
#define MINY (M_PI/2.-DTH)
#define MAXY (M_PI/2.+DTH)
#define MINZ -1.
#define MAXZ 1. 
#define TNX 60 // Total number of cells in X 
#define TNY 60
#define TNZ 1
#define NTX 1 //number of tiles in X 
#define NTY 1
#define NTZ 1
#define SPECIFIC_BC
//#define PERIODIC_XBC
//#define PERIODIC_YBC
//#define PERIODIC_ZBC
/************************************/
//output
/************************************/
#define SILOOUTPUT 1 //to silo file
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define SILO2D_XZPLANE
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1.e10 //stop after this number of steps
#define NOUTSTOP 5000 //stop after this number of outputs
#define DTOUT1 10. //res
#define DTOUT2 100. //avg

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
#define UURHORATIOMIN 1.e-8
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
#define UUAMB 1.e-2 //temp ~ uu/rho
#define RBLOB 15.
#define BW 3.
#define BRHO 1000. //density 

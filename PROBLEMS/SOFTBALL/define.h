/************************************/
//restart
/************************************/
#define RESTART 
#define RESTARTNUM -1
#define RESTARTGENERALINDICES
#define BHDISK_PROBLEMTYPE
#define DOFIXUPS 0
//#define UNPERTURBED

/************************************/
//radiation
/************************************/
//#define RADIATION
//#define BALANCEENTROPYWITHRADIATION

/************************************/
//magnetic fields
/************************************/
//#define MAGNFIELD
//#define VECPOTGIVEN //we provide vector potential
//#define MAXBETA .01 //target max pgas/pgas
//#define BETANORMFULL //normalize everywhere

/************************************/
//coordinates / resolution
/************************************/

#define MYCOORDS KERRCOORDS
#define RMIN 3.5
#define RMAX 17
#define MINX RMIN
#define MAXX RMAX 
#define DTH .6
#define MINY (M_PI/2.-DTH)
#define MAXY (M_PI/2.+DTH)
#define MINZ 0.
#define MAXZ (2.*M_PI) 
#define TNX 40 // Total number of cells in X 
#define TNY 40
#define TNZ 40
#define NTX 1 //number of tiles in X 
#define NTY 1
#define NTZ 4
#define SPECIFIC_BC
#define PERIODIC_ZBC

/************************************/
//Output
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
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//problem params
/************************************/
#define MASS 10.
#define TORUSENTR 5.e2
#define RHOAMB 1.e-30
#define UUAMB 1.e-2*RHOAMB //temp ~ uu/rho
#define ERADATMMIN 1.e-25

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e6
#define EEUURATIOMIN 1.e-15
#define EEUURATIOMAX 1.e6
#define ERADLIMIT 1.e-50
#define RHOFLOOR (RHOAMB/10.)
#define GAMMAMAXRAD 20.

/************************************/
//physics
/************************************/
#define GAMMA (4./3.)
#define FRACMICHEL 0.06

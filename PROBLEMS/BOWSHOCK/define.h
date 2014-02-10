/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 11

/************************************/
//radiation choices
/************************************/
//#define RADIATION

/************************************/
//hydro choices
/************************************/
#define ALLOWENTROPYU2P 1
#define DOFIXUPS 0

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define GDETIN 1

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
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

#define MYCOORDS MINKCOORDS
#define MINX 0.
#define MAXX 1. 
#define MINY 0.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 2. 
#define NX 30 
#define NY 30
#define NZ 30

#define SPECIFIC_BC
/*
#define PERIODIC_XBC
#define PERIODIC_YBC
#define PERIODIC_ZBC
*/


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
#define SIGMACLOUD .25



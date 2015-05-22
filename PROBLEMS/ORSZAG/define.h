#define NONRELMHD

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1

/************************************/
//radiation
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE

/************************************/
//magn. field
/************************************/
#define MAGNFIELD

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define TIMESTEPPING RK2

/************************************/
//viscosity choices
/************************************/
//#define HDVISCOSITY NOVISCOSITY
//#define RADVISCOSITY NOVISCOSITY

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e2
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.6

/************************************/
//coordinates / resolutio
/************************************/
#define MYCOORDS MINKCOORDS
#define MINX 0.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.
#define TNX 64
#define TNY 64
#define TNZ 1
#define PERIODIC_XBC
#define PERIODIC_YBC
#define PERIODIC_ZBC

/************************************/
//output
/************************************/
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e50
#define NOUTSTOP 50
#define SILOOUTPUT 1
#define AVGOUTPUT 1
#define OUTOUTPUT 0
#define DTOUT1 100.
#define DTOUT2 3000.

/************************************/
//common physics / atmosphere
/************************************/

#define GAMMA (4./3.)
#define CSCALE 100.
#define BZERO (1./sqrt(4.*M_PI))
#define VECPOTGIVEN

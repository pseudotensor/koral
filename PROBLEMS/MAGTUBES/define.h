/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 30

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
#define MINMOD_THETA 1.5

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
//coordinates / resolution
/************************************/
#define MYCOORDS MINKCOORDS//KERRCOORDS
#define MINX -0.5
#define MAXX 0.5
#define NX 800
#define NY 1
#define NZ 1
#define MINY (0.95*Pi/2.)
#define MAXY (1.05*Pi/2.)
#define MINZ -1.
#define MAXZ 1.
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC

/************************************/
//output
/************************************/
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 50
#define DTOUT1 0.02

/************************************/
//common physics / atmosphere
/************************************/


#define TUBE 4

#if(TUBE==1) //Sod
#define GAMMA (1.4)
#define BX 0.
#define RHOL 1.
#define VXL 0.
#define VYL 0.
#define VZL 0.
#define PL 1.0
#define BYL 0.
#define BZL 0.
#define RHOR 0.125
#define VXR 0.
#define VYR 0.
#define VZR 0.
#define PR 0.1
#define BYR 0.
#define BZR 0.
#endif

#if(TUBE==4) //Brio & Wu
#define GAMMA (2.)
#define BX 0.75
#define RHOL 1.
#define VXL 0.
#define VYL 0.
#define VZL 0.
#define PL 1.0
#define BYL 1.
#define BZL 0.
#define RHOR 0.125
#define VXR 0.
#define VYR 0.
#define VZR 0.
#define PR 0.1
#define BYR -1.
#define BZR 0.
#endif

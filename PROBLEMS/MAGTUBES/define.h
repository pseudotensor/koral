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
#define MINMOD_THETA 2.
#define FLUXMETHOD LAXF_FLUX
//#define WAVESPEEDSATFACES
#define TIMESTEPPING RK2 //time stepping

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
#define MINX 0.
#define MAXX 1.
#define NY 1
#define NZ 1
#define MINY (0.95*Pi/2.)
#define NX 512
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

/************************************/
//common physics / atmosphere
/************************************/


#define TUBE 3

#if(TUBE==1) //Sod
#define DTOUT1 0.02
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

#if(TUBE==2) //Brio & Wu, RJ 5a
#define CSCALE 100.
#define DTOUT1 0.5
#define GAMMA (5./3.)
#define BX (0.75/CSCALE)
#define RHOL 1.
#define VXL 0.
#define VYL 0.
#define VZL 0.
#define PL (1.0/CSCALE/CSCALE)
#define BYL (1./CSCALE)
#define BZL 0.
#define RHOR 0.125
#define VXR 0.
#define VYR 0.
#define VZR 0.
#define PR (0.1/CSCALE/CSCALE)
#define BYR (-1./CSCALE)
#define BZR 0.
#endif

#if(TUBE==3) //Brio & Wu, RJ 2a
#define CSCALE 100.
#define FAC sqrt(4.*M_PI)
#define DTOUT1 0.5
#define GAMMA (5./3.)
#define BX (2./FAC/CSCALE)
#define RHOL 1.08
#define VXL (1.2/CSCALE)
#define VYL (0.01/CSCALE)
#define VZL (0.5/CSCALE)
#define PL (0.95/CSCALE/CSCALE)
#define BYL (3.6/FAC/CSCALE)
#define BZL (2./FAC/CSCALE)
#define RHOR 1.
#define VXR 0.
#define VYR 0.
#define VZR 0.
#define PR (1./CSCALE/CSCALE)
#define BYR (4./FAC/CSCALE)
#define BZR (2./FAC/CSCALE)
#endif

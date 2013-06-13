/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 0

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE
//#define SKIPRADWAVESPEEDLIMITER
#define ALLOW_EXPLICIT_RAD_SOURCE 0

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
//#define FLUXMETHOD HLL_FLUX
//#define WAVESPEEDSATFACES 
#define GDETIN 0

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY SIMPLEVISCOSITY
#define ALPHATOTALPRESSURE
#define ALPHAHDVISC .1
#define RMINVISC 3.
#define RADVISCOSITY SHEARVISCOSITY
#define TAUSUPPRESSPARAM 100. //the larger the less prad
#define ALPHARADVISC 1.
//#define ENFORCERADWAVESPEEDS

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
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MKS1COORDS
#define MINY (0.01*Pi/2.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define RADOUTPUTINZAMO
//#define PRINTINSIDEBH
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define CGSOUTPUT

/************************************/
//common physics 
/************************************/
#define GAMMA (4./3.)

/************************************/
//model choice
/************************************/
#define NDISK 100

/************************************/
#if (NDISK==100) //mdot = 10, alpha=0.1
/************************************/
#define MKS1R0 -2.
#define RKEP 15.
#define ROUT 30.
#define INJTHETA (0.45*M_PI)
#define MINX (log(1.1-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define NX 40
#define NY 20
#define NZ 1
#define DTOUT1 5.e0
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

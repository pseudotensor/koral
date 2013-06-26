/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 30

/************************************/
//radiation choices
/************************************/
//#define RADIATION
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
//#define SKIPRADSOURCE
//#define SKIPRADWAVESPEEDLIMITER
//#define PUREAXISOUTFLOW

/************************************/
//viscosity choices
/************************************/
//#define HDVISCOSITY SIMPLEVISCOSITY
//#define ALPHATOTALPRESSURE
//#define HDVISCOSITY SHEARVISCOSITY
//#define SHEARVISCOSITYONLYRPHI
//#define ALPHAHDVISC .1
//#define RMINVISC 3.
//#define ZEROTIMEINSHEAR
//#define RADVISCOSITY SHEARVISCOSITY
//#define TAUSUPPRESSPARAM 100. //the larger the less prad
//#define ALPHARADVISC 1.
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
#define MKS1R0 -2.
#define MYCOORDS MKS1COORDS
#define MINX (log(5.1-MKS1R0))
#define MAXX (log(55.3-MKS1R0))
#define NX 40
#define NY 20
#define NZ 1
#define MINY (0.05*Pi/2.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
//#define PRINTINSIDEBH
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define CGSOUTPUT
#define CBAUTOSCALE

/************************************/
//common physics / atmosphere
/************************************/
#define ELL 4.5
#define GAMMA (5./3.)
#define URIN 0.5
#define KKK 9.e-4//1.e-4
#define UTPOT .99
#define RHOATMMIN  3.e-3
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e9,RHOATMMIN))
#define DTOUT1 10.e-1






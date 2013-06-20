/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 5

/************************************/
//radiation choices
/************************************/
#define RADIATION
#define SKIPRADSOURCE

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
#define MINMOD_THETA 1.
#define GDETIN 0

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
//#define ALPHATOTALPRESSURE
//#define RMINVISC 2.
#define RADVISCOSITY NOVISCOSITY
//#define TAUSUPPRESSPARAM 100. //the larger the less prad
//#define ALPHARADVISC 1.
//#define ENFORCERADWAVESPEEDS

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-55
#define EERHORATIOMAX 1.e55
#define EEUURATIOMIN 1.e-55
#define EEUURATIOMAX 1.e55
#define ERADLIMIT 1.e-50
#define RHOFLOOR 1.e-50
#define GAMMAMAXRAD 50.

/************************************/
//blackhole
/************************************/
#define MASS 1.

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MKS1COORDS
#define MKS1R0 -10.
#define BHSPIN 0.0
#define MINX (log(10.-MKS1R0))
#define MAXX (log(100.-MKS1R0))
#define MINY (.5*Pi/2.)
#define MAXY (1.5*Pi/2.)
#define MINZ 0.
#define MAXZ (2.*Pi)
#define NX 200
#define NY 1
#define NZ 1
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define RADOUTPUTINZAMO
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define CGSOUTPUT
#define DTOUT1 20.

/************************************/
//common physics 
/************************************/
#define GAMMA (5./3.)
#define RHOATMMIN 1.
#define UINTATMMIN 1.e-3
#define ERADATMMIN 1.e-14//(1.e-10*calc_LTE_EfromT(calc_PEQ_Tfromurho(UINTATMMIN, RHOATMMIN)))
#define BLOBR 60.
#define BLOBW 10.
#define BLOBP 0.

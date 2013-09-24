/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 1

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
#define BHSPIN 0.0
#define MASS 4.3e6

/************************************/
//coordinates / resolution
/************************************/
#define USEMKS1COORDS
#ifdef USEKSCOORDS
#define MYCOORDS KSCOORDS
#define MINX 500.
#define MAXX 6000.
#define MKS1R0 -20.
#endif
#ifdef USEMKS1COORDS
#define MYCOORDS MKS1COORDS
#define MINX (log(1000.-MKS1R0))
#define MAXX (log(50000.-MKS1R0))
#define MKS1R0 -3000.
#endif


#define MINY (0.05*Pi)
#define MAXY (.95*Pi)
#define MINZ 0.
#define MAXZ (2.*Pi)
#define FULLPHI //for output - to print extra cell in phi
#define NZ 80
#define NX (NZ/2)
#define NY (NZ/8+1)//(NZ/2+1)
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define SILOOUTPUT 1
#define SIMOUTPUT 1
#define OUTCOORDS KERRCOORDS                                                                    
//#define RADOUTPUTINZAMO
//#define PRINTINSIDEBH
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define CGSOUTPUT
#define DTOUT1 1.e3

#define EQPLANEOUTPUT
//#define VERTPLANEOUTPUT

#ifdef VERTPLANEOUTPUT
#define ZSLICE (NZ*0.1)
#endif

#ifdef EQPLANEOUTPUT
#define YSLICE NY/2
#define PRINTZONEMORE
#endif

/************************************/
//problem params
/************************************/
#define GAMMA (5./3.)
#define RHOATMMIN 1.
#define UINTATMMIN 1.e-4
#define TRACER
#define MINTRACE 1.e-4 //edge for enforcing stationary atmosphere

//#define SPHCLOUD
#define G2CLOUD

#define IANGLE (M_PI/3.)
#define OMANGLE 0.//M_PI/2.
#define EARTHMASS (5.97219e27)
#define MASSCLOUD (3.*EARTHMASS)
#define RMINFORPART 5000.
#define KERNELWIDTH 2000.
#define DISKDAMPPAR1 0.5
#define DISKDAMPPAR2 0.1

#define CLRHO 1.e-9
#define CLWIDTH 800.

//#define DONUT
#define IMPOSEDRHO
#define TORUSRIN 10.
#define TORUSKAPPA 0.01
#define TORUSXI 0.708
#define TORUSRBREAK1 45.
#define TORUSRBREAK2 100000.

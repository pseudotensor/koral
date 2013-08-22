/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 36

/************************************/
//radiation choices
/************************************/
#define RADIATION
//#define SKIPRADSOURCE

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define VECPOTGIVEN
#define MAXBETA 0.3 //close to the target pmag/pgas

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2K1K2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
#define RADVISCOSITY NOVISCOSITY
//#define RADVISCOSITY SHEARVISCOSITY
#define ALPHARADVISC 1.

/************************************/
//rmhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 1000.
#define GAMMAMAXRAD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define MKS1R0 -1.

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MINX (log(1.4-MKS1R0))
#define MAXX (log(40.-MKS1R0))
#define NX 60
#define NY 60
#define NZ 1
#endif

#define MINY (0.05*Pi/2.)
#define MAXY (Pi-0.05*Pi/2.)
//#define MAXY (Pi/2.)
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define SILOOUTPUT 
#define SILO2D_XZPLANE
#define CBAUTOSCALE


/************************************/
//common physics / atmosphere
/************************************/
#define GAMMA (4./3.)
#define NODONUT 0
#define INFLOWING 0
#define ELL 3.85
#define URIN (0.)
#define KKK 3000.
#define UTPOT .965
#define DTOUT1 1.e1
#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

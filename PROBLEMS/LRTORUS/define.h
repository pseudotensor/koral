/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
#define RADIATION
//#define SKIPRADSOURCE

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
#define MAXBETA .01 //target pmag/pgas int the midplane

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2K1K2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY NOVISCOSITY
#define RADVISCOSITY NOVISCOSITY
//#define RADVISCOSITY SHEARVISCOSITY
//#define ALPHARADVISC 1.

/************************************/
//rmhd floors
/************************************/
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.9

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define MKS1R0 0.

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MINX (log(1.2-MKS1R0))
#define MAXX (log(100.-MKS1R0))
#define NX 200
#define NY 120
#define NZ 1
#endif

#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)
//#define MAXY (Pi/2.) //change in postinit.c
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
#define SILOOUTPUT 1
#define OUTOUTPUT 1
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 10.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (4./3.)

#define LT_KAPPA 2.e3
#define LT_XI 0.835
#define LT_R1 25.
#define LT_R2 1000.
#define LT_GAMMA 4./3.
#define LT_RIN 10.

#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

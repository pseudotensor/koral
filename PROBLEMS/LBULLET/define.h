/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
//#define RADIATION

/************************************/
//magnetic choices
/************************************/
//#define MAGNFIELD
#define GDETIN 1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX //test IMEX with radiation etc!!!
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//viscosity choices
/************************************/
//#define RADVISCOSITY SHEARVISCOSITY

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
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
#define MINX (log(10.-MKS1R0))
#define MAXX (log(1000.-MKS1R0))
//total resolution
#define TNX 164
#define TNY 92
#define TNZ 1
//number of tiles
#define NTX 4
#define NTY 2
#define NTZ 1
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
//#define OUTPUTPERCORE
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define AVGOUTPUT 0
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 50.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)
#define LT_GAMMA (4./3.)

#define LT_KAPPA 1.e-2

#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#define LT_RIN 22.

#define RHOATMMIN 1.e-10
#define UINTATMMIN 1.e-12

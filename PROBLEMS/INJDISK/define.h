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
#define RADIATION
//#define SKIPRADSOURCE

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define DOFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
#define ALPHARADVISC 1.
#define MAXRADVISCVEL 1.

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 1
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
#define myMKS2COORDS
#define ROUT 200.

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MKS1R0 0.
#define MYCOORDS MKS1COORDS
#define MINX (log(1.2-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild with more cells towards the eq.plane
#define MKS2R0 0.
#define MKS2H0 0.7
#define MYCOORDS MKS2COORDS
#define MINX (log(1.2-MKS2R0))
#define MAXX (log(ROUT-MKS2R0))
#define MINY (0.0025)
#define MAXY (1.-0.0025)
#endif

//total resolution
#define TNX 198
#define TNY 192
#define TNZ 1
//number of tiles
#define NTX 9
#define NTY 4
#define NTZ 1


//#define HALFY
#ifdef HALFY
#undef NY
#undef MAXY
#define NY 20
#define MAXY (Pi/2.)
#endif

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

#define OUTPUTPERCORE
//#define RESOUTPUT_ASCII //by default binary
#define SCAOUTPUT 1
#define SILOOUTPUT 1
#define RADOUTPUT 1
#define OUTOUTPUT 0
#define AVGOUTPUT 1
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 50.
#define DTOUT2 250.
//#define PRINTXGC_RIGHT

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#define DISKHR (1./3.)
#define DISKH (DISKHR*ROUT)
#define DISKSIGMA surfdensCGS2GU(5.e3)
#define DISKRHO DISKSIGMA/DISKH/(96./35.)
#define DISKRCIR 50.
#define DISKVR (-0.02*pow(ROUT/50.,-1.5)) //normalized to 0.01 at R=50
#define MAGNOMEGA (1.5e-2*pow(200./50.,-1.5)) //omega should follow the free fall time?
#define MAGBETA 0.02

//KTORUS only:
#define VSZERO sqrt(1./ROUT/ROUT/ROUT*DISKH*DISKH*3.3) //determines H/R
#define NPOLI 3.
#define ELLA 0.2

//atmosphere
#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

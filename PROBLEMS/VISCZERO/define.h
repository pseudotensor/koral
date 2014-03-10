/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM 10
#define RESTARTGENERALINDICES

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
#define RADVISCMFPCONST 2.
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 1.
#define MAXRADVISCVEL .5

/************************************/
//radiation choices
/************************************/
#define RADIATION
//#define COMPTONIZATION
//#define NCOMPTONIZATION
//#define RADIMPCONV 1.e-5
//#define MAXDIFFTRADS 1000.

/************************************/
//magnetic choices
/************************************/
//#define MAGNFIELD
//#define VECPOTGIVEN
//#define MAXBETA .01 //target pmag/pgas int the midplane

/************************************/
//dynamo
/************************************/
//#define MIMICDYNAMO
//#define ALPHAFLIPSSIGN                                                        
//#define ALPHADYNAMO 0.03

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX //test IMEX with radiation etc!!!
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define GDETIN 1

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
#define BHSPIN 0.

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define MKS1R0 0.

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MINX (log(4.6-MKS1R0))
#define MAXX (log(1000.-MKS1R0))
//total resolution
#define TNX 120
#define TNY 50
#define TNZ 1
//number of tiles
#define NTX 8
#define NTY 8
#define NTZ 1
#endif

#define MINY (0.0025*Pi/2.)
//#define MAXY (Pi-0.0025*Pi/2.)
#define HALFTHETA
#define MAXY (Pi/2.) //change in postinit.c


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
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define AVGOUTPUT 1
#define SIMOUTPUT 0
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 10.
#define DTOUT2 50.
#define NOUTSTOP 50 //max n of outputs                                                                                                                                                                                        

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#define NTORUS 2

#if(NTORUS==2) //for Yucong
#define LT_KAPPA 5.e2
#define EXPECTEDHR 0.4
#define LT_XI 0.95
#define LT_R1 12.5
#define LT_R2 500.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 8.
#endif

#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

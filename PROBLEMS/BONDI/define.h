/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1

//#define CURVETEST
//#define MODYFIKUJKRZYSIE 0
#define GDETIN 1

/************************************/
//radiation
/************************************/
#define RADIATION
//#define RADOUTPUTINFF
#define RADOUTPUTVELS

/************************************/
//coordinates / resolution
/************************************/
#define MKS1R0 0.
//#define MYCOORDS MKS1COORDS
#define MYCOORDS MKS1COORDS
#define OUTCOORDS BLCOORDS
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT
#define PRINTINSIDEBH
#define RMIN 1.6
#define RMAX 20000.

#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#define MINZ -1.
#define MAXZ 1.

#define TNX 128
#define TNY 1
#define TNZ 1
#define NTX 1
#define NTY 1
#define NTZ 1

#define SPECIFIC_BC

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-12
#define RADIMPMAXITER 25

/************************************/
//output
/************************************/
#define SILOOUTPUT 0 //to silo file
#define OUTOUTPUT 1 //to out file
#define AVGOUTPUT 1
#define RADOUTPUT 1
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1.e10 //stop after this number of steps
#define NOUTSTOP 1000 //stop after this number of outputs
#define DTOUT1 1.e2 //res
#define DTOUT2 1.e2 //avg
#define TMAX 1.e10 //time to stop

/************************************/
//test specific
/************************************/
#define TESTNO 0
#define ERADRES 1.e-10
#define PRADGASINIT 1.e-10 
#define RBONDI 0.

#if (TESTNO==0)
#define MDOT 1.e3
#define TGAS0 1.e8
#endif

#if (TESTNO==1)
#define PRADGAS 1.2e-7
#define TGAS0 1e5
#define MDOT 10.
#endif

#if (TESTNO==2)
#define PRADGAS 1.2e-4
#define TGAS0 1.e6
#define MDOT 10.
#endif

#if (TESTNO==3)
//#define LIKEINFRAGILE
#define PRADGAS 1.2e-1
#define TGAS0 1e7
#define MDOT 10.
#endif

#if (TESTNO==4)
#define PRADGAS 1.2e-5
#define TGAS0 1e6
#define MDOT 100.
#endif 

#define MASS 3.
#define BHSPIN 0.
#define MDOTEDD 2.23/16.*1.e18*MASS //cm/s
#define RHOAMB 1.e-25
#define TAMB 1.e5
#define BONDIGAMMA (0.9*5./3.)
#define GAMMA (5./3.)
//#define GAMMA (long double)(1.+1./3.*((1.+PRADGAS)/(.5+PRADGAS)))
#define MUGAS .5


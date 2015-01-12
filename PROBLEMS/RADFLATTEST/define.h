/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM 6
#define MODYFIKUJKRZYSIE 0

/************************************/
//radiation
/************************************/
#define RADIATION
#define SKIPRADSOURCE
#define TWOBEAMS

#define myVET

#ifdef myVET

#define RADCLOSURE VETCLOSURE
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES
#define RADCLOSURECOORDS MINKCOORDS
#define SOCCERBALL 0
#if (SOCCERBALL==0) 
#define NUMANGLES 80
#endif
#if(SOCCERBALL==1)
#define USE3ANGLELOOKUP
#define NUMANGLES 80
#endif
#if(SOCCERBALL==2)
#define USE3ANGLELOOKUP
#define NUMANGLES 160
#endif
#if(SOCCERBALL==3)
#define USE3ANGLELOOKUP
#define NUMANGLES 48
#endif
#endif

#define ALLOWRADCEILINGINIMPLICIT
#define BASICRADIMPLICIT
//#define NCOMPTONIZATION 
//#define RADOUTPUTVELS
#define RADOUTPUTINFF

//#define RADIMPLICITFIXVEL

/************************************/
//coordinates / resolution
/************************************/
#define MKS1R0 0.
#define MYCOORDS MINKCOORDS
#define OUTCOORDS MINKCOORDS

#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT
//#define PRINTYGC_LEFT
//#define PRINTYGC_RIGHT
#define MINX 0.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define TNX 40
#define TNY 5
#define TNZ 1
#define NTX 1 //for MPI and OMP
#define NTY 1
#define NTZ 1

#define SPECIFIC_BC
#define PERIODIC_YBC
#define PERIODIC_ZBC
#define VELRAD 2.


/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .3
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0      

#define DOFIXUPS 0
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 15
#define GAMMAMAXRAD 5.

/************************************/
//output
/************************************/
#define SILOOUTPUT 0 //to silo file
#define OUTOUTPUT 1 //to out file
#define AVGOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 0
#define ALLSTEPSOUTPUT 0 //whether to output every step
//#define NSTEPSTOP 5 //stop after this number of steps
#define NOUTSTOP 10000 //stop after this number of outputs
#define DTOUT1 .1 //res
#define DTOUT2 1.e40 //avg
#define TMAX 1.e100 //time to stop

/************************************/
//test specific
/***********************************/
#define GAMMA (5./3.)


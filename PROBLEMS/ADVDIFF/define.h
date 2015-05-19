/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM 90
#define MODYFIKUJKRZYSIE 0

/************************************/
//radiation
/************************************/
#define RADIATION
#define MASS 1.e1
//#define SKIPRADSORCE
//#define TWOBEAMS
#define RHOAMB 1.
#define PULSEMAG 10.
#define UUAMB (RHOAMB*1.e-4)
#define TAMB (1.e6)
#define VELX 0.05

//#define myVET

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
//#define BASICRADIMPLICIT
#define RADIMPLICITTHRESHOLD 1.e40
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

//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
//#define PRINTYGC_LEFT
//#define PRINTYGC_RIGHT
#define MINX -50.
#define MAXX 50.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define TNX 512
#define TNY 1
#define TNZ 1
#define NTX 1 //for MPI and OMP
#define NTY 1
#define NTZ 1

#define PERIODIC_XBC
#define PERIODIC_YBC
#define PERIODIC_ZBC


/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0      

#define DOFIXUPS 0
#define DORADFIXUPS 0
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-12
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 50
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
#define DTOUT1 100. //res
#define DTOUT2 1.e40 //avg
#define TMAX 1.e100 //time to stop

/************************************/
//test specific
/***********************************/
#define GAMMA (5./3.)


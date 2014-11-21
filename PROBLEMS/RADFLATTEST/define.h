/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1
#define MODYFIKUJKRZYSIE 0

/************************************/
//define MSTEPS
/************************************/
//#define MSTEP
#define MSTEP_LIMITBC
#define NUMMSTEPLEVELS 20

/************************************/
//radiation
/************************************/
#define RADIATION

#define RADCLOSURE VETCLOSURE
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES


#define ALLOWRADCEILINGINIMPLICIT
#define BASICRADIMPLICIT
//#define NCOMPTONIZATION 
#define RADOUTPUTVELS

//#define RADIMPLICITFIXVEL

/************************************/
//coordinates / resolution
/************************************/
#define MKS1R0 0.
#define MYCOORDS MINKCOORDS
#define OUTCOORDS MINKCOORDS
#define RADCLOSURECOORDS MYCOORDS

#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT
#define MINX 0.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define TNX 256
#define TNY 1
#define TNZ 1
#define NTX 1 //for MPI and OMP
#define NTY 1
#define NTZ 1

//#define SUBZONES
#define OUTPUTAFTERSUBZONES
#define NSUBZONES 2
#define SUBZONESOVERLAP 4

#define SPECIFIC_BC

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define SHUFFLELOOPS 0      

#define DOFIXUPS 0
#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-8
#define RADIMPEPS 1.e-8
#define RADIMPMAXITER 15
#define GAMMAMAXRAD 50.

/************************************/
//output
/************************************/
#define SILOOUTPUT 0 //to silo file
#define OUTOUTPUT 1 //to out file
#define AVGOUTPUT 0
#define RADOUTPUT 0
#define SCAOUTPUT 0
#define ALLSTEPSOUTPUT 0 //whether to output every step
//#define NSTEPSTOP 5e4 //stop after this number of steps
#define NOUTSTOP 10000 //stop after this number of outputs
#define DTOUT1 .1 //res
#define DTOUT2 1.e40 //avg
#define TMAX 1.e100 //time to stop

/************************************/
//test specific
/***********************************/
#define GAMMA (5./3.)


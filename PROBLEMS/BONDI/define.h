/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1
#define MODYFIKUJKRZYSIE 0

/************************************/
//radiation
/************************************/
#define RADIATION
#define ALLOWRADCEILINGINIMPLICIT
#define BASICRADIMPLICIT
#define NCOMPTONIZATION 
#define RADOUTPUTVELS

/************************************/
//coordinates / resolution
/************************************/
#define MKS1R0 0.
#define MYCOORDS MKS1COORDS
#define OUTCOORDS BLCOORDS
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT
#define PRINTINSIDEBH
#define RMIN 10.
#define RMAX 1.e5

#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))


#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#define MINZ -1.
#define MAXZ 1.

#define TNX 64
#define TNY 1
#define TNZ 1
#define NTX 16
#define NTY 1
#define NTZ 1

//#define SUBZONES
#define OUTPUTAFTERSUBZONES
#define NSUBZONES 2
#define SUBZONESOVERLAP 8

#define SPECIFIC_BC
#define FIX_TEMPERATURE
//#define FIX_PRESSURERHO
//#define FIX_VELBONDI
#define FIX_VELOUTBONDI


/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define SHUFFLELOOPS 0      

#define DOFIXUPS 1
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
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1e10 //stop after this number of steps
#define NOUTSTOP 10000 //stop after this number of outputs
#define DTOUT1 1.e5 //res
#define DTOUT2 1.e40 //avg
#define TMAX 1.e100 //time to stop

/************************************/
//test specific
/***********************************/
#define GAMMA (5./3.)
#define MDOT 1.e2
#define TAMB 1.e9*(1.e4/RMAX)
//#define UURHORATIOMIN (calc_PEQ_ufromTrho(TAMB,1.))


#define PRADGASINIT 1.e-10
#define MASS 10.
#define BHSPIN 0.
#define MDOTEDD 2.23/16.*1.e18*MASS //cm/s
#define RHOAMB 1.e-25
#define MUGAS 1.
//#define GAMMA (long double)(1.+1./3.*((1.+PRADGAS)/(.5+PRADGAS)))


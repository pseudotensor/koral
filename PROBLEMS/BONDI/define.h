/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM 1
//#define MODYFIKUJKRZYSIE 1
//#define FLUXMETHOD HLL_FLUX
//#define TEST124
//#define GDETIN 1
//#define WAVESPEEDSATFACES


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
//#define EXPLICIT_LAB_RAD_SOURCE

//#define OVERWRITERADWAVESPEEDSWITHHD
//#define RESETNPH
//#define SKIPRADSOURCE
//#define PUTNFFLOOR

#define RADIMPLICITTHRESHOLD 1.e-2
#define RADIMPCONV 1.e-12
#define RADIMPEPS 1.e-8
#define U2PCONV 1.e-12
#define ALLOWRADCEILINGINIMPLICIT
#define BASICRADIMPLICIT
#ifdef RADIATION
#define NCOMPTONIZATION
//#define DAMPCOMPTONIZATIONATBH
#endif

//#define SKIPFANCYOPACITIES
//#define OPACSIMPLE
#define RADOUTPUTVELS

/************************************/
//coordinates / resolution
/************************************/
#define myMKSCOORDS
//#define FLAT

#define TKST0 0.1
#define MKSR0 0.
#define MKSH0 0.1
#define MKSMY1 0.001
#define MKSMY2 0.2
#define MKSMP0 1.5
#define METRICAXISYMMETRIC

#define RMIN 10.
#define RBONDI 1.e6 //(TAMB=3.267e12/RBONDI)
#define RMAX (RBONDI*10.)
#define RMAXOUT RBONDI 

#ifdef myTKSCOORDS //modified Kerr-Shild further from axis
#define MYCOORDS TKS3COORDS
#define METRICNUMERIC
#define METRICTIMEDEPENDENT
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#ifdef myMKSCOORDS
#define MYCOORDS MKS1COORDS
//#define METRICNUMERIC
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#endif

#ifdef myMSPHCOORDS
#define MYCOORDS MSPH1COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#endif


#define MINZ -1.
#define MAXZ 1.

#define TNX 128
#define TNY 1
#define TNZ 1
#define NTX 4//for MPI and OMP
#define NTY 1
#define NTZ 1

#define OUTCOORDS BLCOORDS
#define PRINTXGC_LEFT
#define PRINTGC_RIGHT
#define PRINTINSIDEBH

#define SELFTIMESTEP
#define SELFTIMESTEP_POWRADIUS 1.25
//#define SHORTERTIMESTEP

//#define SUBZONES
#define SUBZONES_NSTEPSTEP 10
//#define OUTPUTAFTERSUBZONES
#define NSUBZONES (4.)
#define SUBZONESOVERLAP 0

#define SPECIFIC_BC
//#define FIX_TEMPERATURE
//#define FIX_PRESSURERHO
//#define FIX_VELBONDI
#define FIX_VELOUTBONDI
#define FULLBONDI
//#define INFLOW


/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2
#define TSTEPLIM (get_tsteplimiter())
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define SHUFFLELOOPS 0      

#define DOFIXUPS 1


#define RADIMPMAXITER 15
#define GAMMAMAXRAD 3.
/************************************/
//output
/************************************/
#define SILOOUTPUT 0 //to silo file
#define OUTOUTPUT 1 //to out file
#define AVGOUTPUT 0
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1e40 //stop after this number of steps
#define NOUTSTOP 10000 //stop after this number of outputs
#ifdef SELFTIMESTEP
#define DTOUT1 (RMAX*1000.*TSTEPLIM) //res
#else
#define DTOUT1 (10000000.)//(RMAX*100.) //res
#endif
#define DTOUT2 (DTOUT1*100000.) //avg
#define TMAX 1.e100 //time to stop

/************************************/
//test specific
/***********************************/
#define GAMMA (5./3.)
#define MDOT 1.e3
//#define TAMB (1.e8*(1.e5/RTEMPOUT))

#ifndef FULLBONDI
#define TAMB (3.267e12/RBONDI)
#else
#define TAMB (3.267e12/RBONDI)
#endif


//#define UURHORATIOMIN (calc_PEQ_ufromTrho(TAMB,1.))


#define PRADGASINIT 1.e-5
#define MASS 10.
#define BHSPIN 0.
#define MDOTEDD 2.23/16.*1.e18*MASS //cm/s
#define RHOAMB 1.e-25
#define MUGAS 1.
//#define GAMMA (long double)(1.+1./3.*((1.+PRADGAS)/(.5+PRADGAS)))

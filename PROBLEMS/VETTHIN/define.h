//#define RESTART
#define RESTARTNUM -1


#define RADIATION

#define myVET

#ifdef myVET
#define RADCLOSURE VETCLOSURE
#define VETFLUXCOSACCEPT 1.0
#define VETFEACCEPT 0.0
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES
#define SOCCERBALL 0
#if (SOCCERBALL==0) 
#define USEDUALNEIGHBOR
#define NUMANGLES 80
#endif
#if(SOCCERBALL==1)
#define NUMANGLES 80
#endif
#if(SOCCERBALL==2)
#define NUMANGLES 160
#endif
#if(SOCCERBALL==3)
#define NUMANGLES 48
#endif
#endif

#define BEAM1
//#define DISK
//#define HOURGLASS


#define DOFIXUPS 0
#define SKIPRADSOURCE

#ifndef myVET
#define RADVISCOSITY SHEARVISCOSITY
#endif

#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define ALPHARADVISC 0.1
#define RADVISCMAXVELDAMP
#define MAXRADVISCVEL 0.5


#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 20.


#define OMSCALE 1.

//#define myMKS1COORDS
#define myMKER1COORDS
//#define myMSPH1COORDS
//#define mySPHCOORDS

#ifdef myMSPH1COORDS
#define MYCOORDS MSPH1COORDS
#define RADCLOSURECOORDS SPHCOORDS
#endif

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#define RADCLOSURECOORDS BLCOORDS
#endif

#ifdef mySPHCOORDS
#define MYCOORDS SPHCOORDS//KSCOORDS//SPHCOORDS//KERRCOORDS
#define RADCLOSURECOORDS SPHCOORDS
#endif

#ifdef myMKER1COORDS
#define MYCOORDS MKER1COORDS
#define RADCLOSURECOORDS BLCOORDS
#endif


#define OUTCOORDS SPHCOORDS

#define IMAGETYPE "gif"

#define OUTVEL VEL4

#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 80
#define RADOUTPUTINZAMO
#define PRINTGC_RIGHT
#define OUTOUTPUT 0
#define SILOOUTPUT 1
#define SILO2D_XZPLANE
#define SIMOUTPUT 0

#ifdef DISK
#define DTOUT1 1.
#define RMIN 4.
#define RMAX 30.
#endif

#ifdef HOURGLASS
#define DTOUT1 1.
#define RMIN 2.2
#define RMAX 5.
#define RBEAM1 2.8
#define RBEAM2 3.2
#endif

#ifdef BEAM1
#define DTOUT1 .1
#define RMIN 2.2//4.//2.2
#define RMAX 4.//7.//5.
#define RBEAM1 2.8//4.8//2.8
#define RBEAM2 3.2//5.2//3.2
#endif

#ifdef myMSPH1COORDS
#define MKS1R0 0.
#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))
#define TNX 60
#endif

#ifdef myMKS1COORDS
#define MKS1R0 0.
#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))
#define TNX 60
#endif

#ifdef myMKER1COORDS
#define MKS1R0 0.
#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))
#define TNX 60
#endif

#ifdef mySPHCOORDS
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#define TNX 60
#endif

#define TNY 60
#define TNZ 1

#define NTX 2
#define NTY 2
#define NTZ 1


#define MINY 0.02
#define MAXY Pi/2
#define HALFTHETA
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC
//#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5
#define UTPOT .98
//#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define TSTEPLIM .4
#define TIMESTEPPING RK2IMEX
#define FLUXLIMITER 0
#define MINMOD_THETA 1.


#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

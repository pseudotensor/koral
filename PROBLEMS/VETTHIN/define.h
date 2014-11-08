//#define RESTART
#define RESTARTNUM 0


#define RADIATION

//#define myVET

#ifdef myVET
#define RADCLOSURE VETCLOSURE
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES
#endif

//#define BEAM1
#define DISK
#define DOFIXUPS 0
#define SKIPRADSOURCE

#define RADVISCOSITY SHEARVISCOSITY
//#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define ALPHARADVISC 0.1
#define RADVISCMAXVELDAMP
#define MAXRADVISCVEL 0.1


#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 20.


#define OMSCALE 1.
//#define HOURGLASS

//#define myMKS1COORDS
//#define myMKER1COORDS
#define myMSPH1COORDS
//#define mySPHCOORDS

#ifdef myMSPH1COORDS
#define MYCOORDS MSPH1COORDS
#endif

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#endif

#ifdef mySPHCOORDS
#define MYCOORDS SPHCOORDS//KSCOORDS//SPHCOORDS//KERRCOORDS
#endif

#ifdef myMKER1COORDS
#define MYCOORDS MKER1COORDS
#endif

#define RADCLOSURECOORDS SPHCOORDS
#define OUTCOORDS SPHCOORDS

#define IMAGETYPE "gif"

#define OUTVEL VEL4
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 2250.
#define RADOUTPUTINZAMO
#define PRINTGC_RIGHT
#define OUTOUTPUT 0
#define SILOOUTPUT 1
#define SILO2D_XZPLANE
#define SIMOUTPUT 0

#ifdef DISK
#define RMIN 4.
#define RMAX 30.
#endif

#ifdef HOURGLASS
#define RMIN 2.2
#define RMAX 5.
#define RBEAM1 2.8
#define RBEAM2 3.2
#endif

#ifdef BEAM1
#define RMIN 2.2
#define RMAX 5.
#define RBEAM1 2.8
#define RBEAM2 3.2
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
#define TNX 40
#endif

#ifdef myMKER1COORDS
#define MKS1R0 0.
#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))
#define TNX 40
#endif

#ifdef mySPHCOORDS
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#define TNX 40
#endif

#define TNY 40
#define TNZ 1

#define NTX 2
#define NTY 2
#define NTZ 1


#define MINY 0.001
#define MAXY Pi/2.
#define HALFTHETA
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC
#define CORRECT_POLARAXIS
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
#define FLUXLIMITER 0
#define MINMOD_THETA 1.


#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

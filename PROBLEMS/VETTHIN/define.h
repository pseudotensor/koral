#define RESTART
#define RESTARTNUM -1


#define RADIATION

#define myVET

#ifdef myVET
#define RADCLOSURE VETCLOSURE
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES
#endif


#define DOFIXUPS 1
#define SKIPRADSOURCE

//#define RADVISCOSITY SHEARVISCOSITY
#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL .5


#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 20.

//#define myMKS1COORDS
#define myMSPH1COORDS
#define OMSCALE 1.
//#define HOURGLASS

#ifdef myMSPH1COORDS
#define MYCOORDS MSPH1COORDS
#else
#define MYCOORDS SPHCOORDS//KSCOORDS//SPHCOORDS//KERRCOORDS
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

#ifdef myMSPH1COORDS
#define MKS1R0 0.
#define MINX (log(2.-MKS1R0))
#define MAXX (log(50.-MKS1R0))
#define TNX 50
#else
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 40.//27.8
#define TNX 40
#endif

#define TNY 40
#define TNZ 1

#define NTX 4
#define NTY 2
#define NTZ 1


#define MINY (0.02*Pi/4.)
#define MAXY Pi/2.
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
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5


#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

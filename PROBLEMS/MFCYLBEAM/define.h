#define RADIATION
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define myMCYL1COORDS

#define MULTIRADFLUID
#define MFCORRECTPHI
#define MFREDISTRIBUTEMETHOD 2
#define MORERADFLUIDS
#define NRF 6

#define OMSCALE 1.

#ifdef myMCYL1COORDS
#define MYCOORDS MCYL1COORDS
#else
#define MYCOORDS CYLCOORDS
#endif

#define OUTCOORDS CYLCOORDS

#define IMAGETYPE "gif"
#define OUTVEL VELPRIMRAD
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000
#define RADOUTPUTINZAMO
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT


#ifdef myMCYL1COORDS
#define MKS1R0 -1.
#define MINX (log(0.001-MKS1R0))
#define MAXX (log(10.-MKS1R0))
#define NX 50
#else
#define MINX  0.001
#define MAXX 10.
#define NX 50
#endif


#define NY 1
#define NZ 5
#define YZXDUMP



#define MINY -1.
#define MAXY 1.

#define MINZ 0.

//#define FULLPHI

#ifdef FULLPHI
#define MAXZ 2.*Pi
#define PRINTZONEMORE
#else
#define MAXZ Pi/2.5
#endif

#define SPECIFIC_BC

#define GAMMA (4./3.)


#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define RK2_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

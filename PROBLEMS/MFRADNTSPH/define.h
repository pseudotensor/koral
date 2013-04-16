#define RADIATION
//#define RADSOURCEOFF
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

//#define myMKS1COORDS

#define OMSCALE 1.

#define MULTIRADFLUID
#define MFCORRECTPHI
#define MFSKEW 15.
#define MFMINVEL 1.e-4
#define MFFRACMETHOD 1
#define MFFRACSKEW 1.
#define MFFRAC 0.05

#define MFREDISTRIBUTEMETHOD 2
#define MORERADFLUIDS
#define NRF 6

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#define OUTCOORDS KERRCOORDS
#else
#define MYCOORDS SPHCOORDS//KSCOORDS//SPHCOORDS//KERRCOORDS
#define OUTCOORDS MYCOORDS
#endif


#define IMAGETYPE "gif"
#define OUTVEL VEL4
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
//#define CGSOUTPUT
#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF
//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT
//#define PRINTYGC_LEFT

#ifdef myMKS1COORDS
#define MKS1R0 0.
#define MINX (log(1.9-MKS1R0))
#define MAXX (log(15.-MKS1R0))
#define NX 30
#else
#define MINX (.5*r_horizon_BL(BHSPIN))
#define MAXX 15.//27.8
#define NX 50
#endif

#define NY 30
#define NZ 1


#define MINY (0.1*Pi/4.)
#define MAXY Pi/2.

#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5
#define UTPOT .98
//#define RHOATMMIN  rhoCGS2GU(1.e-4)
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

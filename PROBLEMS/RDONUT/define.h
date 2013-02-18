#define MASS 10.
#define BHSPIN 0.
#define MYCOORDS KSCOORDS
#define OUTCOORDS KERRCOORDS
#define OUTVEL VELR
#define DTOUT1 10.
#define ALLSTEPSOUTPUT 0

#define NX 150
#define NY 50
#define NZ 1
#define MINX (.8*r_horizon_BL(BHSPIN))
#define MAXX 27.8
#define MINY 0.*Pi/4.
#define MAXY Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define SPECIFIC_BC

#define GAMMA (4./3.)
#define KKK 1.
#define ELL 4.5
#define UTPOT 1.
#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e9,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e6))

#define INT_ORDER 1
#define RK3_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define NODONUT 0
#define INFLOWING 0
#define NSTEPSTOP 50e10
#define NOUTSTOP 1000.
#define CGSOUTPUT

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40


#define RADIATION
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.


#define DISCRETESRC
#define FULLPHI

#define MYCOORDS MINKCOORDS

#define IMAGETYPE "gif"
#define OUTVEL VELPRIMRAD
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000
#define RADOUTPUTINZAMO
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT


#define MINX  -20.
#define MAXX 20.
#define MINY  -20.
#define MAXY 20.
#define NX 80
#define NY 80
#define NZ 1

#define MINZ 0.
#define MAXZ 1.

#define SPECIFIC_BC

#define GAMMA (4./3.)
#define OMSCALE 1.

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

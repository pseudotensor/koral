#define RADIATION
#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define MYCOORDS CYLCOORDS

#define IMAGETYPE "gif"
#define OUTVEL VELPRIMRAD
#define DTOUT1 10.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 100
#define RADOUTPUTINZAMO
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT


#define MINX  0.
#define MAXX 20.

#define NX 200
#define NY 1
#define NZ 1


#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.
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

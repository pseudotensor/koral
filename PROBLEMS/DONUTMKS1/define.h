#define BHSPIN 0.
#define MYCOORDS MKS1COORDS
//#define MYCOORDS KSCOORDS
#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define DTOUT1 1.e1
#define ALLSTEPSOUTPUT 0

#define NX 50
#define NY 50
#define NZ 1
#define MINX 0.25
#define MAXX 3.5
//#define MINX 4.
//#define MAXX 35.
#define MKS1R0 .5
#define MINY 0.*Pi/4.
#define MAXY Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define SPECIFIC_BC
#define PRINTGC_LEFT
#define PRINTGC_RIGHT

#define GAMMA (4./3.)
#define KKK 0.03
#define ELL 4.5
#define UTPOT 1.
#define RHOATMMIN  1.e-4
#define UINTATMMIN  1.e-6

#define INT_ORDER 1
#define RK3_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define NODONUT 0
#define INFLOWING 0
#define NSTEPSTOP 50e10
#define NOUTSTOP 1000.

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40


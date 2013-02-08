#define BHSPIN 0.0
#define MYCOORDS MKS1COORDS
//#define MYCOORDS KSCOORDS
#define OUTCOORDS KERRCOORDS

#define DTOUT1 1.e-1
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 51
#define PRINTGC_LEFT
#define PRINTGC_RIGHT
#define OUTVEL VEL3

#define NX 100
#define NY 1
#define NZ 1
//#define MINX (.8*r_horizon_BL(BHSPIN))
//#define MINX 4.
//#define MAXX 20.
#define MINX .2
#define MAXX 3.
#define MKS1R0 .5

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.
#define MINZ -Pi/20.
#define MAXZ Pi/20.

#define RK3STEPPING
#define INT_ORDER 1
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define SPECIFIC_BC
//#define COPY_XBC
//#define COPY_YBC
//#define COPY_ZBC

#define GAMMA (5./3.)
#define RHOATMMIN 1.e0
#define UINTATMMIN 1.e-2

#define RHOFLOOR 1.e-20
#define UFLOOR 1.e-15
#define EFLOOR 1.e-40

#define BHSPIN 0.
#define MYCOORDS SCHWCOORDS
#define OUTCOORDS SCHWCOORDS

#define DTOUT1 1.e2
#define ALLSTEPSOUTPUT 0
#define NOUTSTOP 51
#define PRINTGC_LEFT
#define PRINTGC_RIGHT

#define NX 50
#define NY 1
#define NZ 1
#define MINX 100.
#define MAXX 200.
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
#define PAR_D 1.e0
#define PAR_U 1.e-3

#define RHOFLOOR 1.e-20
#define UFLOOR 1.e-15
#define EFLOOR 1.e-40

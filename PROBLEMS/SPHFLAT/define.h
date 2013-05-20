#define VERBOSE0 0
#define BHSPIN 0.
#define MYCOORDS SPHCOORDS

#define TIME_STEPPING RK2
#define INT_ORDER 1
#define NX 20
#define NY 20
#define NZ 1

#define TSTEPLIM .5
#define NOUTSTOP 30

#define SPECIFIC_BC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 5.e0
#define ALLSTEPSOUTPUT 1

#define GAMMA (ldouble)(5./3.)
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define MINX 6.
#define MAXX 40.
#define MINY 0.1*Pi/2.
#define MAXY Pi/2.
#define MINZ 0.
#define MAXZ 1.

#define RHO_AMB 1.e-3
#define U_AMB 1.e-7

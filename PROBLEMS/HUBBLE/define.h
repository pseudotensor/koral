#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0

#define MYCOORDS MINKCOORDS

#define INT_ORDER 1
#define NY 1
#define NZ 1

#define ENFORCEENTROPY
#define TSTEPLIM .5
#define TIMESTEPPING RK2
#define NX 512


#define INITTSTEPLIM (TSTEPLIM/10.)

#define SPECIFIC_BC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define CALCL1_HUBBLE
#define TMAX 1.e4

#define ALLSTEPSOUTPUT 0

#define GAMMA (1.4)
#define MINX -1.
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.


#define RHO0 1.
#define VPRIME 1.e-5
#define UINT0 1.e-10
#define DTOUT1 1.e3
#define NOUTSTOP 101

#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT

#define EFLOOR 1.e-50
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

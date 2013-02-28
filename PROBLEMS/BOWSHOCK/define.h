#define U2PPREC 1.e-5
#define U2PRADPREC 1.e-5
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define BHSPIN 0.

//efine CYLINDRICAL
#define SPHERICAL

#ifdef CYLINDRICAL
#define MYCOORDS CYLCOORDS
#define MINY -1.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 2.*Pi
#define YZXDUMP
#define PRINTZONEMORE
#endif

#ifdef SPHERICAL
#define MYCOORDS SPHCOORDS
#define YSLICE NY-1
#define MINY 0*Pi/2.
#define MAXY 1.*Pi/2.
#define MINZ 0.
#define MAXZ Pi
#endif


#define VELINX -0.1
#define UINTFRAC .001

#define RK3STEPPING
#define INT_ORDER 1
#define NX 20
#define NY 20
#define NZ 20
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define NOUTSTOP 1000
//#define U2P_TEMP

#define SPECIFIC_BC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define DTOUT1 30.
#define ALLSTEPSOUTPUT 0
#define GAMMA (ldouble)(5./3.)

#define MINX 6.
#define MAXX 40.


#define RHOFLOOR 1.e-20
#define UFLOOR 1.e-15
#define PAR_D 1.e0
#define PAR_E 1.e-4

#define RHO_AMB 1.e-3
#define U_AMB 1.e-7
#define EFLOOR 1.e-40
//#define LOGXGRIDREF
//#define LOGPAR1 2.5
//#define LOGPAR2 2.

//#define PRINTGC_LEFT
//#define PRINTGC_RIGHT


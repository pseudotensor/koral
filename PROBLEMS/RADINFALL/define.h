#define BHSPIN 0.
#define MYCOORDS SCHWCOORDS
#define ANAL_PROFILE
#define INT_ORDER 1
#define NX 100
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define NOUTSTOP 30
//#define U2P_TEMP

#define SPECIFIC_BC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 5.e0
#define ALLSTEPSOUTPUT 0
#define GAMMA (ldouble)(5./3.)
#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
#define MINX 6.
#define MAXX 40.
#define MINY 0*Pi/2.
#define MAXY 1.*Pi/2.
#define MINZ 0.
#define MAXZ 1.
#define RHOFLOOR 1.e-20
#define UFLOOR 1.e-15
#define PAR_D 1.e0
#define PAR_E 1.e-2

#define RHO_AMB 1.e-3
#define U_AMB 1.e-4
#define EFLOOR 1.e-40
//#define LOGXGRIDREF
//#define LOGPAR1 2.5
//#define LOGPAR2 2.
//#define PRINTGC_LEFT

#define BLOB
#define BLOBX 30.
#define BLOBMAG 1.
#define BLOBSIG 1.

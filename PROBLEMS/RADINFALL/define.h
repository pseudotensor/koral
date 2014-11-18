#define BHSPIN 0.

#define myMKS

#ifdef myKERR
#define MYCOORDS BLCOORDS
#define RMIN 5.
#define RMAX 40.
#define MINX RMIN
#define MAXX RMAX
#endif

#ifdef myMKS
#define MYCOORDS MKER1COORDS
//#define METRICNUMERIC
#define MKSR0 0
#define RMIN 2.01
#define RMAX 30.
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(RMAX-MKSR0))
#endif

#define OUTCOORDS BLCOORDS 

#define MINY .99*Pi/2.
#define MAXY 1.01*Pi/2.

#define ANAL_PROFILE
#define INT_ORDER 1
#define TNX 512
#define TNY 1
#define TNZ 1
#define NTX 4
#define NTY 1
#define NTZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)
#define NOUTSTOP 30
#define OUTOUTPUT 1 //to out file
//#define U2P_TEMP

#define SPECIFIC_BC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define DTOUT1 5.e0
#define ALLSTEPSOUTPUT 0
#define GAMMA (ldouble)(5./3.)
//#define COPY_XBC
#define COPY_YBC
#define COPY_ZBC
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
#define BLOBX 20.
#define BLOBMAG 10.
#define BLOBSIG .1

//#define LABRADFLUXES

#define TMAX 200.
#define RADIATION

#define MYCOORDS MINKCOORDS

#define NX 150
#define NY 50
#define NZ 1
#define ANGLE .3


#define MINX -6.
#define MAXX 3.//3.
#define MINY 0.
#define MAXY 3.//1.5

#define BEAMY (MAXY/4.)//.3


#define TSTEPLIM 1.//kind of courant limiter
#define INT_ORDER 1
#define RK3STEPPING
#define INITTSTEPLIM (TSTEPLIM/10.)//for the 1st time step
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define ALLSTEPSOUTPUT 0
#define GAMMA (1.4)
#define DTOUT1 1.

//#define IMPLICIT_FF_RAD_SOURCE
//#define EXPLICIT_RAD_SOURCE
#define MASS 1./MSUNCM


#define MINZ -1.
#define MAXZ 1.

#define RHOFLOOR 1.e-50
#define UFLOOR 1.e-65
#define EFLOOR 1.e-40

#define TAMB 1.e7
#define TLEFT 1.e9

#define RHOAMB 1.e-4

//#define BLOB
#define RHOBLOB 1.e3

#define BLOBW 5.e-2
#define KAPPA 0.//.1

#define NLEFT 0.999

#define SPECIFIC_BC

#define RADOUTPUTINZAMO
//#define PRINTGC_LEFT

#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-6
#define RADFORCEPREC 1.e-5

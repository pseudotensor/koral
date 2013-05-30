#define RADIATION
#define SKIPRADSOURCE
#define SIMPLERADVISCOSITY
#define ALPHARADVISC 1.
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.

#define myMCYL1COORDS

//#define DISCRETESRC
//#define FULLPHI

#ifdef myMCYL1COORDS
#define MYCOORDS MCYL1COORDS
#else
#define MYCOORDS CYLCOORDS
#endif

#define OUTCOORDS CYLCOORDS

#define IMAGETYPE "gif"
#define OUTVEL VELPRIMRAD
#define DTOUT1 10.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000
#define RADOUTPUTINZAMO
//#define RADOUTPUTVELS
#define PRINTXGC_LEFT
#define PRINTXGC_RIGHT

#define GAMMAMAXRAD 10.

#ifdef myMCYL1COORDS
#define MKS1R0 -1.
#define MINX (log(0.01-MKS1R0))
#define MAXX (log(10.-MKS1R0))
#define NX 100
#else
#define MINX  .0001
#define MAXX 10.
#define NX 200
#endif


#define NY 1
#ifndef FULLPHI
#define NZ 1
#else
#define NZ 100
#endif
#define YZXDUMP



#define MINY -1.
#define MAXY 1.

#define MINZ 0.



#ifdef FULLPHI
#define MAXZ 2.*Pi
#define PRINTZONEMORE
#else
#define MAXZ Pi/2.5
#endif

#define SPECIFIC_BC

#define GAMMA (4./3.)
#define OMSCALE 1.

#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define RK2_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

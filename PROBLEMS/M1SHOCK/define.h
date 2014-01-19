#define INT_ORDER 1
#define TSTEPLIM .2
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
#define TIMESTEPPING RK2IMEX

//#define RESTART
#define RESTARTNUM 1

#define RADIATION 
//#define SKIPRADSOURCE

#define RADPERTM1CONV
#define RADPERTM1DELTA 1.

#define RADVISCOSITY SHEARVISCOSITY
//#define NUMRADWAVESPEEDS
#define RADVISCNUDAMP
#define ALPHARADVISC 1.
#define MAXRADVISCVEL 1.
#define RADVISCMFPCONST (.5/NX)
#define ZEROTIMEINSHEAR

#define MYCOORDS MINKCOORDS //metric

/************************************/
//rmhd floors
/************************************/
#define UURHORATIOMIN 1.e-50
#define UURHORATIOMAX 1.e50
#define EERHORATIOMIN 1.e-50
#define EERHORATIOMAX 1.e50
#define EEUURATIOMIN 1.e-50
#define EEUURATIOMAX 1.e50
#define GAMMAMAXRAD 5000.
#define GAMMAMAXHD 5000.

#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF

#define NX 200 //x-resolution
#define NY 1 //y-resolution
#define NZ 1 //z=rezolution

#define MINX 0. //size of the domain
#define MAXX 1.
#define MINY -1.
#define MAXY 1.
#define MINZ -1.
#define MAXZ 1.

#define SPECIFIC_BC //whether the boundary conditions are more complicated and given in PROBLEMS/XXX/bc.c

#define DTOUT1 0.1 
#define NOUTSTOP 50
#define ALLSTEPSOUTPUT 0 
#define PRINTXGC_LEFT //if x-left ghost cells are to be printed out, 
#define PRINTXGC_RIGHT //if x-right ghost cells are to be printed out 

#undef SIGMA_RAD
#define SIGMA_RAD 1.e-45
#define GAMMA (5./3.)
#define TEMPIN 1.e7
#define TEMPOUT 2.e8
#define FBEAM 0.9

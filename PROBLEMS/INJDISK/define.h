/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 8

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define SKIPRADSOURCE

/************************************/
//magnetic choices
/************************************/
#define MAGNFIELD
#define GDETIN 1
//#define VECPOTGIVEN
#define MAGBETA .01 //pmag/pgas at the boundary

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX 
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define DOFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
#define ZEROTIMEINSHEAR
#define ALPHARADVISC 1.
#define MAXRADVISCVEL 1./3.

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 1
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define MKS1R0 0.
#define ROUT 50.
#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MINX (log(1.575-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define NX 150
#define NY 100
#define NZ 1
#endif


#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)

//#define HALFY
#ifdef HALFY
#undef NY
#undef MAXY
#define NY 20
#define MAXY (Pi/2.)
#endif

#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define AVGOUTPUT 1
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 10.
#define DTOUT2 250.
#define PRINTXGC_RIGHT

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)


#define DISKRHO rhoCGS2GU(1.e-2)
#define DISKRCIR 20.
#define DISKVR (-1./sqrt(ROUT)*0.1) //v_kepl * alpha
#define VSZERO sqrt(1./ROUT/ROUT/ROUT*DISKH*DISKH*3.3) //determines H/R
#define NPOLI 3.
#define ELLA 0.2
#define MAGNOMEGA 1.e-2

//outdated
#define DISKH (.2*ROUT)
#define DISKTEMP 1.e9

#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

//************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
//#define RADIATION
//#define COMPTONIZATION

/************************************/
//magnetic choices
/************************************/
//if we want a magnetic field, uncomment MAGNFIELD
//#define MAGNFIELD
#define GDETIN 1

#define MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#ifdef RADIATION
#define TIMESTEPPING RK2IMEX
#else
#define TIMESTEPPING RK2HEUN
#endif
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define DOFIXUPS 1

/************************************/
//viscosity choices
/************************************/
//#define RADVISCOSITY SHEARVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define ALPHARADVISC .1
#define RADVISCMAXVELDAMP
#define MAXRADVISCVEL 1.

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 10
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
#define myMKS2COORDS
#define ROUT 100.
#define RMIN 1.8

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MKSR0 0.
#define MYCOORDS MKS1COORDS
#define MINX (log(3.6-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild with more cells towards the eq.plane
#define MKSR0 0.
#define MKSH0 0.6 //makes cells smaller towards equatorial plane
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(ROUT-MKSR0))
#define MINY (0.0025)
#define MAXY (1.-0.0025)
#endif

//total resolution
#define TNX 240
#define TNY 160
#define TNZ 96
//number of tiles
#define NTX 2
#define NTY 2
#define NTZ 1

//#define HALFTHETA //symmetry wrt eq. plane?
#ifdef HALFTHETA
#undef TNY
#undef MAXY
#define TNY 50
#ifdef myMKS1COORDS 
#define MAXY (Pi/2.)
#endif
#ifdef myMKS2COORDS 
#define MAXY (.5)
#endif
#endif



#define PHIWEDGE (2.*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define PERTMAGN 1.e-2
#define SPECIFIC_BC
#define PERIODIC_ZBC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000

#define SCAOUTPUT 1
#define SILOOUTPUT 1
#define RADOUTPUT 1
#define AVGOUTPUT 1
#define COORDOUTPUT 2
#define SILO2D_XZPLANE
#define PRINTXGC_RIGHT

#define DTOUT1 1.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)


//parameters for inflow
#define DISKVR (-sqrt(2./ROUT))
#define DISKHR 0.1
#define DISKHCGS (DISKHR*ROUT*MASS*1477.) 
#define MDOTEDD 10.
#define MDOT (MDOTEDD*2.23e18)
//TODO:correct!!!
#define DISKSIGMA surfdensCGS2GU(MDOT/((-3.e10)*DISKVR*2.*DISKHCGS))
#define DISKRCIR 25.
#define DISKTEMP 1.e8
#define VERTBTIME 10.

#define MAXRADIUS4DYNAMO DISKRCIR

#define DISKH (DISKHR*ROUT)
#define DISKRHO (DISKSIGMA/2./DISKH)

#define MAGNOMEGA 0.//(2.*M_PI/1000.)//0.
//#define MAGNOMEGA 0.*(1.5e-2*pow(200./50.,-1.5)) //omega should follow the free fall time?
#define MAGBETA 0.1

//atmosphere
#define RHOATMMIN  1.e-20
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)


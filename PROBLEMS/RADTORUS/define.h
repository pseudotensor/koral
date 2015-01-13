/************************************/
//restart
/************************************/
#define RESTART 
#define RESTARTNUM -1
#define RESTARTGENERALINDICES
#define BHDISK_PROBLEMTYPE
#define DOFIXUPS 0
#define UNPERTURBED
#define SKIPHDEVOLUTION

/************************************/
//radiation
/************************************/
#define RADIATION
#define U2PCONV 1.e-8
#define RADIMPCONV 1.e-8
#define RADIMPEPS 1.e-6
#define RADIMPMAXITER 15

#define myVET

#ifdef myVET
#define RADCLOSURE VETCLOSURE
#define EVOLVEINTENSITIES
#define RADSTARTWITHM1INTENSITIES
#endif

//#define BALANCEENTROPYWITHRADIATION

/************************************/
//magnetic fields
/************************************/
//#define MAGNFIELD
//#define VECPOTGIVEN //we provide vector potential
//#define MAXBETA .01 //target max pgas/pgas
//#define BETANORMFULL //normalize everywhere

/************************************/
//coordinates / resolution
/************************************/

#define RMIN 2.5
#define RMAX 14.

#define myMKER1COORDS

#ifdef myMKER1COORDS
#define MYCOORDS MKER1COORDS
#define RADCLOSURECOORDS BLCOORDS
#endif

#ifdef myMKER1COORDS
#define MKS1R0 0.
#define MINX (log(RMIN-MKS1R0))
#define MAXX (log(RMAX-MKS1R0))
#endif

#define DTH (0.01)
#define MINY (DTH)
#define MAXY (M_PI-DTH)
#define MINZ -1.
#define MAXZ 1.
#define TNX 50 // Total number of cells in X 
#define TNY 100
#define TNZ 1
#define NTX 2 //number of tiles in X 
#define NTY 2
#define NTZ 1
#define SPECIFIC_BC
#define PERIODIC_ZBC

#define CORRECT_POLARAXIS
#define U2P_EQS U2P_EQS_NOBLE
#define U2P_SOLVER U2P_SOLVER_W
#define NCCORRECTPOLAR 2


/************************************/
//Output
/************************************/
#define SILOOUTPUT 1 //to silo file
#define RADOUTPUT 1
#define SCAOUTPUT 1
#if (TNZ==1)
#define SILO2D_XZPLANE
#endif
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1.e10 //stop after this number of steps
#define NOUTSTOP 9000 //stop after this number of outputs
#define DTOUT1 1. //res
#define DTOUT2 1.e20 //avg

/************************************/
//reconstruction / stepping
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//problem params
/************************************/
#define MASS 10.
#define TORUSENTR 5.e2
#define RHOAMB 1.e-30
#define UUAMB 1.e-2*RHOAMB //temp ~ uu/rho
#define ERADATMMIN 1.e-30
#define RADIMPLICITTHRESHOLD 1.e10

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e6
#define EEUURATIOMIN 1.e-60
#define EEUURATIOMAX 1.e60
#define ERADLIMIT 1.e-50
#define RHOFLOOR (RHOAMB/10.)
#define GAMMAMAXRAD 20.

/************************************/
//physics
/************************************/
#define GAMMA (5./3.)
//when constructing rad-pressure supported torus we want to have pressure like for a gamma=4/3 gas, because radiation pressure has effective gamma = 4/3 
#define EFFGAMMA (4./3.) //opt.thick
#define FRACMICHEL 0.06
#define PERTMAGN 0.01

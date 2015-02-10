/************************************/
//restart
/************************************/
#define RESTART 
#define RESTARTNUM -1
#define RESTARTGENERALINDICES
#define BHDISK_PROBLEMTYPE
#define DOFIXUPS 1
#define UNPERTURBED

//#define SKIPRADSOURCE
//#define SKIPHDEVOLUTION
#define FIXEDALLBUTTEMP
//#define SKIPHDBUTENERGY

//#define HEATINGRATEPERMASS 1.e-4
//#define HEATINGRATEPERMASSSQ (1.e-4/1.e-13)
//#define HEATINGLIMITINGRHO 1.e-14
/************************************/
//radiation
/************************************/
#define RADIATION
#define OMSCALE 1.


/************************************/
//viscosity choices
/************************************/
#ifdef RADIATION
#define RADVISCOSITY SHEARVISCOSITY
//#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC .5
#define MAXRADVISCVEL .5
#endif

#define ALLOWRADCEILINGINIMPLICIT
#define BASICRADIMPLICIT

#define U2PCONV 1.e-10
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-6
#define RADIMPMAXITER 15

//#define myVET

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

#define RMIN 3.5
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

#define DTH (0.04)
#define MINY (DTH)
#define MAXY (M_PI-DTH)
#define MINZ -M_PI
#define MAXZ M_PI
#define TNX 50//150 // Total number of cells in X 
#define TNY 150//300
#define TNZ 1
#define NTX 2 //number of tiles in X 
#define NTY 2
#define NTZ 1
#define SPECIFIC_BC
#define PERIODIC_ZBC

//#define CORRECT_POLARAXIS
#define U2P_EQS U2P_EQS_NOBLE
#define U2P_SOLVER U2P_SOLVER_W
#define NCCORRECTPOLAR 2


/************************************/
//Output
/************************************/
#define SILOOUTPUT 1 //to silo file
#define RADOUTPUT 1
#define SIMOUTPUT 2
#define SCAOUTPUT 1
#define AVGOUTPUT 1
#if (TNZ==1)
#define SILO2D_XZPLANE
#endif
#define ALLSTEPSOUTPUT 0 //whether to output every step
#define NSTEPSTOP 1.e10 //stop after this number of steps
#define NOUTSTOP 9000 //stop after this number of outputs
#define DTOUT1 1. //res
#define DTOUT2 10. //avg

/************************************/
//reconstruction / stepping
/************************************/
#define INT_ORDER 2
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .001
#define FLUXLIMITER 0
#define MINMOD_THETA 2.

/************************************/
//problem params
/************************************/
#define MASS 10.
#define TORUSENTR 1.2e3//5.e2 - higher density, 1.2e3 - order lower density 
//#define TEMPTORUS 5.e6 //overwrites the temperature
#define RHOAMB 1.e-26
#define UUAMB 1.e-4*RHOAMB //temp ~ uu/rho, but see belo
#define TEMPAMB 5.e6 //overwrites the ambient temperature
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
#define GAMMAMAXRAD 5.

/************************************/
//physics
/************************************/
#define GAMMA (5./3.)
//when constructing rad-pressure supported torus we want to have pressure like for a gamma=4/3 gas, because radiation pressure has effective gamma = 4/3 
#define EFFGAMMA (4./3.) //opt.thick
#define FRACMICHEL 0.06
#define PERTMAGN 0.01

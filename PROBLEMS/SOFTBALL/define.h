/************************************/
//restart
/************************************/
//#define RESTART 
#define RESTARTNUM -1
#define RESTARTGENERALINDICES
#define BHDISK_PROBLEMTYPE

#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1

//#define UNPERTURBED
//#define SKIPHDEVOLUTION

//tests
//#define COOLINGTOWARDSENTROPY
//#define TARGETLOGENTROPY 1.5
//#define THETANOCOOL 0.2

//this controls the heating rate per unit mass per time
//this more or less should equal the luminosity of the rad-supported torus divided by the total mass
//#define HEATINGRATEPERMASS 1.

/************************************/
//radiation
/************************************/
//#define RADIATION
#define RADIMPLICITTHRESHOLD 1.e0
#define MAXRADIMPDAMPING 1.e-6
#define BALANCEENTROPYWITHRADIATION
#define COMPTONIZATION
#define ALLOWRADCEILINGINIMPLICIT
#define RADIMPCONV 1.e-8
#define RADIMPEPS 1.e-6
#define RADIMPMAXITER 50
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR 1.e-1
#define RADIMPCONVRELENTR 1.e-4
#define RADIMPCONVRELENTRERR .999

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

#define MYCOORDS KERRCOORDS
#define RMIN 5.
#define RMAX 14.
#define MINX RMIN
#define MAXX RMAX 
#define DTH .3
#define MINY (M_PI/2.-DTH)
#define MAXY (M_PI/2.+DTH)
#define MINZ 0.
#define MAXZ (2.*M_PI) 
#define TNX 80 // Total number of cells in X 
#define TNY 80
#define TNZ 1
#define NTX 2 //number of tiles in X 
#define NTY 2
#define NTZ 1
#define SPECIFIC_BC
#define PERIODIC_ZBC

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
#define DTOUT1 10. //res
#define DTOUT2 1.e20 //avg

/************************************/
//reconstruction / stepping
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

/************************************/
//problem params
/************************************/
#define MASS 10.
#define TORUSENTR 5.e2
#define RHOAMB 1.e-30
#define UUAMB 1.e-2*RHOAMB //temp ~ uu/rho
#define ERADATMMIN 1.e-40

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e6
#define EEUURATIOMIN 1.e-15
#define EEUURATIOMAX 1.e6
#define ERADLIMIT 1.e-50
#define RHOFLOOR (RHOAMB/10.)
#define GAMMAMAXRAD 20.

/************************************/
//physics
/************************************/
#define GAMMA (5./3.)
//when constructing rad-pressure supported torus we want to have pressure like for a gamma=4/3 gas, because radiation pressure has effective gamma = 4/3 
#ifdef RADIATION
#define EFFGAMMA (4./3.) //opt.thick
#else
#define EFFGAMMA (5./3.) //opt.thin
#endif

#define FRACMICHEL 0.06
#define PERTMAGN 0.01

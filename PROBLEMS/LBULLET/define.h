/************************************/
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
#define COMPTONIZATION

/************************************/
//magnetic choices
/************************************/
#define MIMICDYNAMO
#define ALPHAFLIPSSIGN                                                        
#define ALPHADYNAMO 0.03
//#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
#define MAXBETA .01 //target pmag/pgas int the midplane

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX //test IMEX with radiation etc!!!
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 1.

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 50.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#define MKS1R0 0.

#define FULLPHI

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS KERRCOORDS
#define MINX 400.
#define MAXX 4000.
//total resolution
#define TNX 100
#define TNY 128
#define TNZ 200
//number of tiles
#define NTX 4
#define NTY 4
#define NTZ 4
#endif

#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)
//#define MAXY (Pi/2.) //change in postinit.c
#define MINZ 0.
#define MAXZ 2.*Pi
#define SPECIFIC_BC

/************************************/
//output
/************************************/
//#define OUTPUTPERCORE
#define OUTCOORDS KERRCOORDS                                                                    

#define PRINTXGC_LEFT 1

#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define AVGOUTPUT 1
#define SIMOUTPUT 1
//#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 2000.
#define DTOUT2 5000.


/***********************************/
// Bullet choice
/**********************************/
#define PI 3.14159
#define BULLET_PHI
//#define BULLET_THETA


#define BULLETPOS 2500.
#define BULLETRAD 200.
#define BULLETTH PI/5.
#define BULLETRHO 1.e0




/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#define NTORUS 8

#if(NTORUS==8) //flat sigma
#define LT_KAPPA 1.e-3
#define EXPECTEDHR 0.4
#define LT_XI 0.7
#define LT_R1 35.
#define LT_R2 200.
#define LT_GAMMA 5./3.
#define LT_RIN 7.5
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./50.)
#endif

#if(NTORUS==7) //flat sigma
#define LT_KAPPA 3.e2
#define EXPECTEDHR 0.4
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./50.)
#endif

#if(NTORUS==6) //for not-so-hyper
#define LT_KAPPA 1.5e3
#define EXPECTEDHR 0.4
#define LT_XI 0.95
#define LT_R1 16.
#define LT_R2 200.
#define LT_GAMMA 4./3.
//#define LT_RIN 10.25
#define LT_RIN 10.6
#undef MAXBETA
#define MAXBETA (1./25.) //eq.plane
#endif

#if(NTORUS==5) //single toroidal loop

#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e-2
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./30.) //target pmag/pgas inside torus
#define BETANORMFULL
//#define BETANORMFACTOR 2.e-10
#endif

#if(NTORUS==4) //a=0 SANE, no rad, denser loops
#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e-2
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./30.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

#if(NTORUS==3) //a=0 SANE, no rad!
#define EXPECTEDHR 0.4
#define LT_KAPPA 1.e-2
#define LT_XI 0.708
#define LT_R1 42.
#define LT_R2 1000.
#define LT_GAMMA 5./3.
#define LT_RIN 10.
#undef MAXBETA
#define MAXBETA (1./30.) //target pmag/pgas inside torus
#define BETANORMFULL
#endif

#if(NTORUS==1) //original (2nd koral paper)
#define LT_KAPPA 1.5e3
#define LT_XI 0.9
#define LT_R1 31.75
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 15.
#endif

#if(NTORUS==2) //for Yucong?
#define LT_KAPPA 2.e3
#define LT_XI 0.95
#define LT_R1 16.
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 10.
#endif

#define RHOATMMIN  1.e-24
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

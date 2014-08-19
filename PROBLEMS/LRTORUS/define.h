/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM 86

/************************************/
//radiation choices
/************************************/
#define RADIATION
#define BALANCEENTROPYWITHRADIATION

#define RADIMPCONV 1.e-8
#define U2PCONV 1.e-12

//#define OPACBELLLIN
//#define RADCLOSURE MINERBOCLOSURE
//#define RADCLOSUREMIXATBH
//#define RADCLOSURECOORDS BLCOORDS
//#define COMPTONIZATION
//#define DAMPCOMPTONIZATIONATBH


/************************************/
//magnetic choices
/************************************/
//#define MIMICDYNAMO
//#define SKIPHDEVOLUTION //temp!
#define ALPHAFLIPSSIGN                                                        
#define CALCHRONTHEGO
#define ALPHADYNAMO 0.314
#define DAMPBETA
#define BETASATURATED 0.1
#define ALPHABETA 6.28
#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN
#define MAXBETA .01 //target pmag/pgas int the midplane

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 0
#define REDUCEORDERWHENNEEDED
#define REDUCEORDERTEMP 1.e7

/************************************/
//viscosity choices
/************************************/
//#define RADVISCOSITY SHEARVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.5

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
#define MASS 3.e5
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define MKS2R0 1.
#define MKS2H0 0.7

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(1.5-MKS2R0))
#define MAXX (log(100.-MKS2R0))
#define MINY (0.001)
#define MAXY (1.-0.001)
//total resolution
#define TNX 64
#define TNY 64
#define TNZ 1
//number of tiles
#define NTX 1
#define NTY 1
#define NTZ 1
#endif

//efine MINY (0.0025*Pi/2.)
//efine MAXY (Pi-0.0025*Pi/2.)
//#define MAXY (Pi/2.) //change in postinit.c
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
//#define OUTPUTPERCORE
#define OUTCOORDS KERRCOORDS                                                                    
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
#define AVGAVGOUTPUT 1
#define ANARELRADOUTPUT 1
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 1.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#define NTORUS 8

#if(NTORUS==8) //flat sigma
#define LT_KAPPA 1.e1
#define EXPECTEDHR 0.4
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#define LT_GAMMA 4./3.
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (.1) //eq.plane
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
#define ERADATMMIN  (calc_LTE_EfromT(1.e3)/10)

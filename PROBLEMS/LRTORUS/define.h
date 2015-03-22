/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE
//#define PERTURBVEL (-.5)

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1
#define PERTMAGN 1.e-2
//#define ENFORCEENTROPY

/************************************/
//radiation choices
/************************************/
#define BALANCEENTROPYWITHRADIATION
#define RADIATION
#define COMPTONIZATION

/************************************/
//magnetic choices
/************************************/
//#define MIMICDYNAMO
//#define CALCHRONTHEGO
//#define THETAANGLE 0.25
//#define ALPHAFLIPSSIGN                                                        
//#define ALPHADYNAMO 0.314
//#define DAMPBETA
//#define BETASATURATED 0.1
//#define ALPHABETA 6.28

#define MAGNFIELD
#define GDETIN 1
#define VECPOTGIVEN

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2IMEX
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 1
#define DORADFIXUPS 1

/************************************/
//viscosity choices
/************************************/
#define RADVISCOSITY SHEARVISCOSITY
//#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.3

/************************************/
//rmhd floors
/************************************/
//#define CORRECT_POLARAXIS_3D
#define CORRECT_POLARAXIS
//#define POLARAXISAVGIN3D
#define U2P_EQS U2P_EQS_NOBLE
#define U2P_SOLVER U2P_SOLVER_W
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 10000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 10.
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
#define myMKS3COORDS
#define MKSR0 0.
#define MKSH0 0.6
#define MKSMY1 0.001
#define MKSMY2 0.025
#define MKSMP0 1.2
#define METRICAXISYMMETRIC

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(1.85-MKSR0))
#define MAXX (log(1000.-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(1.85-MKSR0))
#define MAXX (log(1000.-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#define PHIWEDGE (2.*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

//total resolution                                                                                                                               
#define TNX 120//264 //12*22                                                                                                                          
#define TNY 80//192 //12*16                                                                                                                          
#define TNZ 4//132 //12*11                                                                                                                          

//number of tiles                                                                                                                                
#define NTX 12//22
#define NTY 4//16
#define NTZ 1//11

#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

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

#define BOXOUTPUT 1
#define BOXR1 15.
#define BOXR2 20.
#define BOXITH 20 
#define VAROUTPUT 1
#define VARRADIUS 100.
#define NVARCUTS 20

#define DTOUT3 1.
#define DTOUT4 1.

#define SILOOUTPUT 1
#define OUTOUTPUT 0
#define RADOUTPUT 1
#define SCAOUTPUT 1
#define AVGOUTPUT 1
#define SILO2D_XZPLANE
#define CBAUTOSCALE
#define DTOUT1 10.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#define NTORUS 7

#if(NTORUS==81) //
#define LT_KAPPA 2.e2
#define LT_XI 0.705
#define LT_R1 40.
#define LT_R2 1000.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 10.
#define BETANORMEQPLANE
//#define BETANORMFACTOR 1.e-04
#undef MAXBETA
#define MAXBETA (.1) 
#endif


#if(NTORUS==79 || NTORUS==80) //
#define LT_KAPPA 5.e2
#define LT_XI 0.705
#define LT_R1 40.
#define LT_R2 1000.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 10.
#define BETANORMEQPLANE
//#define BETANORMFACTOR 1.e-04
#undef MAXBETA
#define MAXBETA (.005) 
#endif

#if(NTORUS==78) //
#define LT_KAPPA 5.e2
#define LT_XI 0.96
#define LT_R1 14.
#define LT_R2 400.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 10.
#define BETANORMEQPLANE
#undef MAXBETA
#define MAXBETA (.1) 
#endif

#if(NTORUS==77) //flat sigma, single poloidal loop
#define LT_KAPPA 5.e2
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (.05) 
#endif

#if(NTORUS==7) //flat sigma
#define LT_KAPPA 2.e2
#define EXPECTEDHR 0.3
#define LT_XI 0.975
#define LT_R1 30.
#define LT_R2 200.
#ifdef RADIATION
#define LT_GAMMA 4./3.
#else
#define LT_GAMMA 5./3.
#endif
#define LT_RIN 22.
#define BETANORMFULL
#undef MAXBETA
#define MAXBETA (1./10.)
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

#define RHOATMMIN  1.e-20
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

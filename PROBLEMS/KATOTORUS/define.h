#define NONRELMHD
#define RADIMPLICITTHRESHOLD 1.e0
#define MAXRADIMPDAMPING 1.e-6
#define NONRELMHDENTROPYCUT 1.e-10 // Tcut = 3e12*this number

/************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/ 
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
#define RADIATION
//#define SKIPRADSOURCE
#define BALANCEENTROPYWITHRADIATION
#define COMPTONIZATION
#define ALLOWRADCEILINGINIMPLICIT
#define RADIMPLICITFIXVEL
#define RADIMPCONVRELERR 1.e-4
//#define BASICRADIMPLICIT
//#define RADIMPSTARTWITHEXP
//#define ALLOWFORENTRINF4DPRIM

//#define U2P_EQS U2P_EQS_JONS
//#define U2P_SOLVER U2P_SOLVER_WP


#define U2PCONV 1.e-12
#define RADIMPCONV 1.e-10
#define RADIMPEPS 1.e-5
#define RADIMPMAXITER 50

/************************************/
//magnetic choices
/************************************/
#define MIMICDYNAMO
#define CALCHRONTHEGO
#define THETAANGLE 0.25
#define ALPHAFLIPSSIGN                                                        

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
#define TSTEPLIM 0.5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define SHUFFLELOOPS 0
#define DOFIXUPS 0
#define DORADFIXUPS 0

/************************************/
//viscosity choices
/************************************/
#ifdef RADIATION
#define RADVISCOSITY SHEARVISCOSITY
#endif
#define ACCELRADVISCOSITY
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define RADVISCMAXVELDAMP
#define ALPHARADVISC 0.1
#define MAXRADVISCVEL 0.1

/************************************/
//blackhole
/************************************/
#define MASS 6.62
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
//#define myMKS2COORDS
#define mySPHCOORDS
//#define myCYLCOORDS
#define RMIN 4.
#define RMAX 100.
#define MKSR0 -300.
#define MKSH0 0.8
#define MKSMY1 0.001
#define MKSMY2 0.2
#define MKSMP0 1.5
#define METRICAXISYMMETRIC

#ifdef myMSPH1COORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS MSPH1COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY (0.45)
#define MAXY (M_PI-MINY)
#endif

#ifdef mySPHCOORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS SPHCOORDS
#define MINX RMIN
#define MAXX 100.
#define MINY (0.05)
#define MAXY (M_PI-MINY)
#endif

#ifdef myCYLCOORDS //modified Kerr-Shild
#define PWPOTENTIAL
#define MYCOORDS CYLCOORDS
#define MINX RMIN
#define MAXX 100.
#define MINY (-60.)
#define MAXY 60.
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY (0.001)
#define MAXY (1.-0.001)
#endif

#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(1.85-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

#define PHIWEDGE (M_PI/2.)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

//total resolution
#define TNX 300//128 //28*9
#define TNY 500//192 //26*9
#define TNZ 1 //2*8
//number of tiles
#define NTX 4
#define NTY 4
#define NTZ 1

#define SPECIFIC_BC
#define PERIODIC_ZBC
//#define PERIODIC_XBC
//#define PERIODIC_YBC

/************************************/
//output
/************************************/

#define DTOUT3 1.
#define BOXOUTPUT 1
#define BOXR1 10.
#define BOXR2 15.
#define BOXITH 30 //distance from eq.plane in cells                                                                                                             
#define OUTCOORDS MYCOORDS//KERRCOORDS                                                                    
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
#if(TNZ==1)
#define SILO2D_XZPLANE
#endif
#define CBAUTOSCALE
#define DTOUT1 1.
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

#ifndef RADIATION
#undef GAMMA
#define GAMMA (4./3.)
#endif

#define NTORUS 1

#if(NTORUS==1) //Jiang+14
#define KT_A 0.4
#define KT_R0 (25.*2.)
#define KT_RHO0 (10.*rhoCGS2GU(1.e-2))
#define KT_T0 (100.*1.e7)
#undef MAXBETA
#define MAXBETA (1.25*3.294/10.) //eq.plane
#endif

#define RHOFLOOR 1.e-50
#define RHOATMMIN  1.e-20
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)


/************************************/
//rmhd floors
/************************************/
#ifndef myCYLCOORDS
#define CORRECT_POLARAXIS
#endif
//#define POLARAXISAVGIN3D
#define NCCORRECTPOLAR 2
#define UURHORATIOMIN 1.e-20
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
#define GAMMAMAXHD 2.

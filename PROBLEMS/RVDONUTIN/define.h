/************************************/
//radiation choices
/************************************/
#define ALLOW_EXPLICIT_RAD_SOURCE 1
#define GAMMAMAXRAD 1000.

/************************************/
//hydro choices
/************************************/
#define ALLOWENTROPYU2P 1
#define FIXUPAFTERENTROPY 0
#define DOFIXUPS 0

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.8
#define FLUXMETHOD HLL_FLUX
//#define WAVESPEEDSATFACES 
#define GDETIN 0
#define SKIPRADSOURCE
//#define SKIPRADWAVESPEEDLIMITER
 
/************************************/
//hd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e3
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e6

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.

/************************************/
//simple viscosity
/************************************/
//#define VISCOSITY
//#define ENFORCERADWAVESPEEDS
#define SIMPLEVISCOSITY
#define ALPHATOTALPRESSURE
#define RMINVISC 4.

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MKS1R0 -2.
#define MINX (log(1.-MKS1R0))
#define NX 40
#define NY 20
#define NZ 1
#else //Schwarzschild
#define MYCOORDS SCHWCOORDS
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 16.
#define NX 48
#define NY 32
#define NZ 1
#endif
#define MINY (0.005*Pi/2.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
//#define OUTCOORDS KERRCOORD                                                                        S
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
///#define RADOUTPUTINZAMO
#define PRINTINSIDEBH
//#define RADOUTPUTINFF
#define CGSOUTPUT
#define PRINTXGC_LEFT
//                                                                                                                                                                                                                                                                                                                                                                                               #define PRINTGC_RIGHT

/************************************/
//common physics / atmosphere
/************************************/
#define GAMMA (4./3.)
#define NODONUT 0
#define INFLOWING 0

/************************************/
//model choice
/************************************/
#define NDONUT 3

/************************************/
#if (NDONUT==12) //mdot = 0.5, alpha = 0.1, r=30
/************************************/
#define MAXX (log(30.-MKS1R0))
#undef NY
#define NY 80
#define RADIATION
#define ELL 5.868
#define ALPHAVISC .1
#define URIN (4.e6/CCC)
#define KKK 131.859
#define UTPOT 0.98424
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)

/************************************/
#elif (NDONUT==11) //mdot = 0.5, alpha = 0.1, r=50
/************************************/
#define MAXX (log(50.-MKS1R0))
#define RADIATION
#define ELL 7.3657
#define ALPHAVISC .1
#define URIN (1.6e6/CCC)
#define KKK 23.41
#define UTPOT 0.99025
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)

/************************************/
#elif (NDONUT==10) //mdot = 100, alpha = 0.1, r=30
/************************************/
#define MAXX (log(30.-MKS1R0))
#define RADIATION
#define ELL 5.868
#define ALPHAVISC .1
#define URIN (7.29e8/CCC)
#define KKK 18455.
#define UTPOT 1.015//.9974
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)

/************************************/
#elif (NDONUT==9) //mdot = 100, alpha = 0.01
/************************************/
#define MAXX (log(50.-MKS1R0))
#define RADIATION
#define ELL 7.3657
#define ALPHAVISC .01
#define URIN (4.00e7/CCC)
#define KKK 5548.
#define UTPOT 1.008//.9980
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-20
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)

/************************************/
#elif (NDONUT==8) //mdot = 100, alpha = 0.1
/************************************/
#define MAXX (log(50.-MKS1R0))
#define RADIATION
#define ELL 7.3657
#define ALPHAVISC .1
#define URIN (4.00e8/CCC)
#define KKK 11954.//4676.
#define UTPOT 1.008//.9980
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)

/************************************/
#elif (NDONUT==7) //mdot = 10, alpha = 0.1
/************************************/
#define MAXX (log(50.-MKS1R0))
#define RADIATION
#define ELL 7.3657
#define ALPHAVISC .1
#define URIN (8.00e7/CCC)
#define KKK 2381.
#define UTPOT .9937
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)

/************************************/
#elif (NDONUT==6) //mdot = 1, alpha = 0.1
/************************************/
#define MAXX (log(50.-MKS1R0))
#define RADIATION
#define ELL 7.3657
#define ALPHAVISC .1
#define URIN (4.80e6/CCC)
#define KKK 490.
#define UTPOT .9912
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-21
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)


/************************************/
#elif (NDONUT==5) //mdot = 100, alpha = 0.1, r=15 (and below)
/************************************/
#define MAXX (log(15.3-MKS1R0))
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (1.57e9/CCC)
#define KKK 9713.
#define UTPOT .9925
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)


/************************************/
#elif (NDONUT==4) //mdot = 10, alpha = 0.01
/************************************/
#define MAXX (log(15.3-MKS1R0))
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .01
#define URIN (3.92e7/CCC)
#define KKK 3300.
#define UTPOT .983
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))

/************************************/
#elif (NDONUT==3) //mdot = 10, alpha = 0.1
/************************************/
//333
//#undef NX
//#define NX 80
//#undef NY
//#define NY 50
#define MAXX (log(15.3-MKS1R0))
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (3.92e8/CCC)
#define KKK 7127.
#define UTPOT .983
#define DTOUT1 5.e-1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10)

/************************************/
#elif (NDONUT==2) //mdot = 1, alpha = 0.1
/************************************/
#define MAXX (log(15.3-MKS1R0))
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (5.23e7/CCC)
#define KKK 1894.
#define UTPOT .9734
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-23
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))

/************************************/
#elif (NDONUT==1) //mdot = 0.5, alpha = 0.1
/************************************/
#define MAXX (log(15.3-MKS1R0))
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (1.57e7/CCC)
#define KKK 291.
#define UTPOT .9704
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))




/************************************/
#elif (NDONUT==0) //pure hd
/************************************/
#define MAXX (log(15.3-MKS1R0))
#define ALPHAVISC .1
#define URIN 0.5
#define KKK 9.e-4//1.e-4
#define UTPOT .99
#define RHOATMMIN  3.e-3
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e9,RHOATMMIN))
#define DTOUT1 10.e-1
#endif






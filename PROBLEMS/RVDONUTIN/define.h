/************************************/
//radiation choices
/************************************/
#define ALLOW_EXPLICIT_RAD_SOURCE 0
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
#define RK2_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

/************************************/
//hd floors
/************************************/
#define UURHORATIOMIN 1.e-7
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-7
#define EERHORATIOMAX 1.e3

/************************************/
//blackhole
/************************************/
#define MASS 10.
#define BHSPIN 0.

/************************************/
//simple viscosity
/************************************/
#define VISCOSITY
#define SIMPLEVISCOSITY
#define ALPHATOTALPRESSURE
#define RMINVISC 2.

/************************************/
//coordinates / resolution
/************************************/
#define myMKS1COORDS
#ifdef myMKS1COORDS //modified Kerr-Shild
#define MYCOORDS MKS1COORDS
#define MKS1R0 -2.
#define MINX (log(1.5-MKS1R0))
#define MAXX (log(15.3-MKS1R0))
#define NX 64
#define NY 32
#define NZ 1
#else //Scharzschild
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
#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define RADOUTPUTINZAMO
#define CGSOUTPUT
//#define PRINTGC_LEFT
#define PRINTGC_RIGHT

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
#if (NDONUT==1) //mdot = 0.5, alpha = 0.1
/************************************/
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (1.57e7/CCC)
#define KKK 291.
#define UTPOT .9704
#define DTOUT1 1.e5
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))

/************************************/
#elif (NDONUT==2) //mdot = 1, alpha = 0.1
/************************************/
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (5.23e7/CCC)
#define KKK 1894.
#define UTPOT .9734
#define DTOUT1 1.e5
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))

/************************************/
#elif (NDONUT==3) //mdot = 10, alpha = 0.1
/************************************/
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .1
#define URIN (3.92e8/CCC)
#define KKK 7127.
#define UTPOT .983
#define DTOUT1 1.e5
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))

/************************************/
#elif (NDONUT==4) //mdot = 10, alpha = 0.01
/************************************/
#define RADIATION
#define ELL 4.5
#define ALPHAVISC .01
#define URIN (3.92e7/CCC)
#define KKK 3300.
#define UTPOT .983
#define DTOUT1 1.e5
#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))

/************************************/
#elif (NDONUT==0) //pure hd
/************************************/
#define ALPHAVISC .1
#define URIN 0.5
#define KKK 9.e-4//1.e-4
#define UTPOT .99
#define RHOATMMIN  3.e-3
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e9,RHOATMMIN))
#define DTOUT1 10.e-1
#endif






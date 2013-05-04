//#define RADIATION

//#define RADSOURCEOFF
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE
#define ALLOW_EXPLICIT_RAD_SOURCE 0

#define ALLOWENTROPYU2P 1
#define ALLOWCOLDU2P 0
#define FIXUPAFTERENTROPY 0
#define DOFIXUPS 0
#define WAVESPEEDSATFACES
//#define AVERAGEONLYUINT
//#define FULLDISSIPATION
#define UURHORATIOMIN 1.e-10
#define UURHORATIOMAX 1.e3
#define EERHORATIOMIN 1.e-7
#define EERHORATIOMAX 1.e3
#define RHORHOMAXRATIOMIN 1.e-20

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define myMKS1COORDS

#ifdef myMKS1COORDS
#define MYCOORDS MKS1COORDS
#else
#define MYCOORDS SCHWCOORDS
#endif

#define VELPRIM VELR
//#define BLOB

//#define VISCOSITY
#define SIMPLEVISCOSITY
#define ALPHATOTALPRESSURE
#define RMINVISC 2.

#define OUTCOORDS KERRCOORDS
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 5000.
#define RADOUTPUTINZAMO
//#define RADOUTPUTINFF
//#define PRINTGC_LEFT
//efine PRINTGC_RIGHT

#define PAR_D 1.e0
#define PAR_E 1.e-8

#ifdef myMKS1COORDS
#define MKS1R0 -2.
#define MINX (log(1.5-MKS1R0))
#define MAXX (log(15.3-MKS1R0))//(log(16.-MKS1R0))
#define NX 64
#else
#define MINX (1.5*r_horizon_BL(BHSPIN))
#define MAXX 16.
#define NX 16//
#endif

#define NY 32
#define NZ 1


#define MINY (0.005*Pi/2.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC
//#define PUREAXISOUTFLOW

#define GAMMA (4./3.)
#define ELL 4.5


#ifdef RADIATION

//mdot = 0.5
#define ALPHAVISC .1
#define URIN (1.57e7/CCC)
#define KKK 291.
#define UTPOT .9704

//mdot = 1
/*
#define ALPHAVISC .1
#define URIN (5.23e7/CCC)
#define KKK 1894.
#define UTPOT .9734
*/

//mdot = 10
/*
#define ALPHAVISC .1
#define URIN (3.92e8/CCC)
#define KKK 7127.
#define UTPOT .983
*/

//mdot = 10
/*
#define ALPHAVISC .01
#define URIN (3.92e7/CCC)
#define KKK 3300.
#define UTPOT .983
*/

#define RHOATMMIN  1.e-22
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6))
#define DTOUT1 1.e5
//#define CGSOUTPUT

#else //purehd

#define ALPHAVISC .1
#define URIN 0.
#define KKK 9.e-4//1.e-4
#define UTPOT .99
#define RHOATMMIN  1.e-6
#define UINTATMMIN 1.e-8

#define DTOUT1 1.e1

#endif


#define INT_ORDER 1
#define RK2_STEPPING
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 2.

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

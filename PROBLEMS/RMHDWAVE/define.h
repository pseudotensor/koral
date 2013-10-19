#define MYCOORDS MINKCOORDS

#define NX 64
#define NY 1
#define NZ 1

#define TSTEPLIM .3

#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 2.
#define INT_ORDER 1

#define GAMMA (ldouble)(5./3.)
#define MINX 0
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

//#define EDDINGTON_APR
#define KAPPAES 0.
#define KK (2.*M_PI)
#define MAGNFIELD

#define RADOUTPUTINZAMO
#define CALCL1_RMHDWAVE
#define TIMESTEPPING RK2


#define NUMERO 10

#if (NUMERO==1) //sonic wave
#define KAPPA 0.
#define RHOFAC 1.e-1
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (1.50023e-7*RHOFAC)
#define DUIM 0.
#define DV1RE (0.00001*RHOFAC)
#define DV1IM 0.
#define DV2RE 0.
#define DV2IM 0.
#define B1ZERO 0.
#define B2ZERO 0.
#define DB2RE 0.
#define DB2IM 0.
#define OMRE 0.0628319
#define OMIM 0.
#define NOUTSTOP 11
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2.*M_PI/OMRE/10.
#endif

#if (NUMERO==10) //fast-magnetosonic wave, no-rad
#define KAPPA 0.
#define RHOFAC 1.e-1
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (1.50023e-7*RHOFAC)
#define DUIM 0.
#define DV1RE (0.0000161788*RHOFAC)
#define DV1IM 0.
#define DV2RE (-9.99788e-6*RHOFAC)
#define DV2IM 0.
#define B1ZERO 0.0100008
#define B2ZERO 0.0100008
#define DB2RE (0.0000161808*RHOFAC)
#define DB2IM 0.
#define OMRE 0.101654
#define OMIM 0.
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==11)
#define RADIATION
#undef SIGMA_RAD
#define SIGMA_RAD 10.
#define PP 0.01
#define CC 1.e2
#define KAPPA 0.01
#define RHOFAC 10.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define DVRE (9.99998e-6*RHOFAC)
#define DVIM (8.48878e-9*RHOFAC)
#define DURE (1.66666e-3*RHOFAC)
#define DUIM (2.82938e-6*RHOFAC)
#define DEERE (-4.52724e-5*RHOFAC)
#define DEEIM (2.78566e-5*RHOFAC)
#define DF1RE (-5.83678e-6*RHOFAC)
#define DF1IM (-9.48194e-6*RHOFAC)
#define DF2RE (-5.83678e-6*RHOFAC)
#define DF2IM (-9.48194e-6*RHOFAC)
#define OMRE 6.28317e-2
#define OMIM 5.33366e-5
#define DTOUT1 1.e-0
#endif

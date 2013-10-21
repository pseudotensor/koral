#define MYCOORDS MINKCOORDS

#define NY 1
#define NZ 1

#define TSTEPLIM .6

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

#define RADIMPCONV 1.e-12
#define U2PCONV 1.e-12

#define NX 512

#define NUMERO 102

#if (NUMERO==1) //sonic wave
#define KAPPA 0.
#define RHOFAC 1.e-2
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
#define RHOFAC 1.e-2
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

#if (NUMERO==11) //slow-magnetosonic wave, no-rad
#define KAPPA 0.
#define RHOFAC 1.e-2
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (1.50023e-7*RHOFAC)
#define DUIM 0.
#define DV1RE (6.18031e-6*RHOFAC)
#define DV1IM 0.
#define DV2RE (0.0000100001*RHOFAC)
#define DV2IM 0.
#define B1ZERO 0.0100008
#define B2ZERO 0.0100008
#define DB2RE (-6.18108e-6*RHOFAC)
#define DB2IM 0.
#define OMRE 0.038832
#define OMIM 0.
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==101) //rad modified, sonic wave, tau=0.01, P=0.01
#define RADIATION
#define KAPPA 0.01
#undef SIGMA_RAD 
#define SIGMA_RAD 2.46889e-42
#define RHOFAC 1.e-1
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (1.5001477607e-7*RHOFAC)
#define DUIM (7.6399157051e-10*RHOFAC)
#define DV1RE (9.999775486e-6*RHOFAC)
#define DV1IM (2.5465054939e-8*RHOFAC)
#define DV2RE (0.)
#define DV2IM (0.)
#define B1ZERO 0.0
#define B2ZERO 0.0
#define DB2RE 0.
#define DB2IM 0.
#define EEZERO 1.8002700405e-6
#define DEERE (3.272468477e-14*RHOFAC)
#define DEEIM (3.44343726e-13*RHOFAC)
#define DF1RE (-2.4099469e-11*RHOFAC)
#define DF1IM (7.58186629475e-12*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.0628304424
#define OMIM 0.000160001659
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11

#endif

#if (NUMERO==102) //rad modified, sonic wave, tau=0.1, P=0.01
#define RADIATION
#define KAPPA 0.1
#undef SIGMA_RAD 
#define SIGMA_RAD 2.46889e-42
#define RHOFAC 1.e-1
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (1.49254e-7*RHOFAC)
#define DUIM (7.547e-9*RHOFAC)
#define DV1RE (9.97763e-6*RHOFAC)
#define DV1IM (2.53022e-7*RHOFAC)
#define DV2RE (0.)
#define DV2IM (0.)
#define B1ZERO 0.0
#define B2ZERO 0.0
#define DB2RE 0.
#define DB2IM 0.
#define EEZERO 1.8002700405e-6
#define DEERE (3.22869e-12*RHOFAC)
#define DEEIM (3.8497e-12*RHOFAC)
#define DF1RE (-3.35009e-11*RHOFAC)
#define DF1IM (7.48081e-11*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.0626913
#define OMIM 0.00158978
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11

#endif

#if (NUMERO==110) //rad modified, fast-msonic wave, tau=0.01, P=0.01
#define RADIATION
#define KAPPA 0.01
#undef SIGMA_RAD 
#define SIGMA_RAD 2.46889e-42
#define RHOFAC 1.e-2
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (1.50019e-7*RHOFAC)
#define DUIM (4.72231e-10*RHOFAC)
#define DV1RE (0.0000161788*RHOFAC)
#define DV1IM (7.03971e-9*RHOFAC)
#define DV2RE (-9.99794e-6*RHOFAC)
#define DV2IM (9.73070e-9*RHOFAC)
#define B1ZERO 0.0100008
#define B2ZERO 0.0100008
#define DB2RE (0.0000161809*RHOFAC)
#define DB2IM (-8.70405e-9*RHOFAC)
#define EEZERO 1.80027e-6
#define DEERE (3.34303e-14*RHOFAC)
#define DEEIM (5.56961e-13*RHOFAC)
#define DF1RE (-3.88935e-11*RHOFAC)
#define DF1IM (7.6322e-12*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.101654
#define OMIM 0.0000442318
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11

#endif

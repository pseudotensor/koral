#define MYCOORDS MINKCOORDS

#define NY 1
#define NZ 1

#define TSTEPLIM .6

#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.
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
//#define EXPLICIT_RAD_SOURCE

#define RADOUTPUTINFF
#define CALCL1_RMHDWAVE
#define TIMESTEPPING RK2

#define RADIMPCONV 1.e-12
#define U2PCONV 1.e-14

#define NX 128

#define NUMERO 1003

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


#if (NUMERO==104) //rad modified, sonic wave, tau=0.1, P=0.1, CC=1e1
#define RADIATION
#define KAPPA 0.1
#undef SIGMA_RAD 
#define SIGMA_RAD 2.36051e-47
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.0000151557*RHOFAC)
#define DUIM (7.69693e-7*RHOFAC)
#define DV1RE (0.0000997992*RHOFAC)
#define DV1IM (2.55207e-6*RHOFAC)
#define DV2RE (0.)
#define DV2IM (0.)
#define B1ZERO 0.0
#define B2ZERO 0.0
#define DB2RE 0.
#define DB2IM 0.
#define EEZERO 0.00182741
#define DEERE (1.33148e-10*RHOFAC)
#define DEEIM (3.60017e-8*RHOFAC)
#define DF1RE (-2.52471e-7*RHOFAC)
#define DF1IM (7.40041e-8*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.627057
#define OMIM 0.0160351
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==105) //rad modified, sonic wave, tau=10, P=10, CC=1e1
#define RADIATION
#define KAPPA 10.
#undef SIGMA_RAD 
#define SIGMA_RAD 2.36051e-45
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.000011707*RHOFAC)
#define DUIM (1.88153e-6*RHOFAC)
#define DV1RE (0.000266251*RHOFAC)
#define DV1IM (0.0000633514*RHOFAC)
#define DV2RE (0.)
#define DV2IM (0.)
#define B1ZERO 0.0
#define B2ZERO 0.0
#define DB2RE 0.
#define DB2IM 0.
#define EEZERO 0.182741
#define DEERE (0.000205419*RHOFAC)
#define DEEIM (0.000149859*RHOFAC)
#define DF1RE (-0.0000207307*RHOFAC)
#define DF1IM (0.0000377556*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 1.6729
#define OMIM 0.398049
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif


#if (NUMERO==1001) //rad modified, fast msonic wave, tau=0.1, P=0.1, CC=1e1
#define RADIATION
#define KAPPA 0.1
#undef SIGMA_RAD 
#define SIGMA_RAD 2.36051e-47
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.0000151983*RHOFAC)
#define DUIM (4.81572e-7*RHOFAC)
#define DV1RE (0.000160252*RHOFAC)
#define DV1IM (7.18765e-7*RHOFAC)
#define DV2RE (-0.0000979957*RHOFAC)
#define DV2IM (1.01386e-6*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (0.000162345*RHOFAC)
#define DB2IM (-9.13696e-7*RHOFAC)
#define EEZERO 0.00182741
#define DEERE (1.48645e-9*RHOFAC)
#define DEEIM (6.06323e-8*RHOFAC)
#define DF1RE (-3.95433e-7*RHOFAC)
#define DF1IM (8.51173e-8*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 1.00689
#define OMIM 0.00451613
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==1101) //rad modified, slow msonic wave, tau=0.1, P=0.1, CC=1e1
#define RADIATION
#define KAPPA 0.1
#undef SIGMA_RAD 
#define SIGMA_RAD 2.36051e-47
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.0000150172*RHOFAC)
#define DUIM (1.22284e-6*RHOFAC)
#define DV1RE (0.0000615366*RHOFAC)
#define DV1IM (1.81827e-6*RHOFAC)
#define DV2RE (0.0000989793*RHOFAC)
#define DV2IM (6.53003e-6*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (-0.0000614822*RHOFAC)
#define DB2IM (-5.89827e-6*RHOFAC)
#define EEZERO 0.00182741
#define DEERE (1.96430e-10*RHOFAC)
#define DEEIM (2.18733e-8*RHOFAC)
#define DF1RE (-1.65186e-7*RHOFAC)
#define DF1IM (7.17813e-8*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.386646
#define OMIM 0.01142450
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==1002) //rad modified, fast msonic wave, tau=10., P=10., CC=1e1
#define RADIATION
#define KAPPA 10.
#undef SIGMA_RAD 
#define SIGMA_RAD 2.36051e-45
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.0000117290*RHOFAC)
#define DUIM (1.70920e-6*RHOFAC)
#define DV1RE (0.0002786280*RHOFAC)
#define DV1IM (0.0000519949*RHOFAC)
#define DV2RE (-0.0000352003*RHOFAC)
#define DV2IM (9.58975e-6*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (0.00011243400*RHOFAC)
#define DB2IM (-5.64666e-6*RHOFAC)
#define EEZERO 0.18274100
#define DEERE (0.000207174*RHOFAC)
#define DEEIM (0.000136067*RHOFAC)
#define DF1RE (-0.0000183034*RHOFAC)
#define DF1IM (0.0000363096*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 1.750670
#define OMIM 0.32669300
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==1102) //rad modified, slow msonic wave, tau=10., P=10., CC=1e1
#define RADIATION
#define KAPPA 10.
#undef SIGMA_RAD 
#define SIGMA_RAD 2.360507338595936e-45
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.009137055837563452
#define DURE (9.57383336055e-6*RHOFAC)
#define DUIM (1.3116722304e-6*RHOFAC)
#define DV1RE (0.00009336446865*RHOFAC)
#define DV1IM (0.000012053581584*RHOFAC)
#define DV2RE (0.0001400622766700*RHOFAC)
#define DV2IM (0.000324581754*RHOFAC)
#define B1ZERO 0.10075854437197568
#define B2ZERO 0.10075854437197568
#define DB2RE (-0.0000923995722968773*RHOFAC)
#define DB2IM (-0.000325350086418*RHOFAC)
#define EEZERO 0.18274111675126906
#define DEERE (0.000034908069853*RHOFAC)
#define DEEIM (0.000104592131222 *RHOFAC)
#define DF1RE (-0.0000349080806*RHOFAC)
#define DF1IM (7.303348358574e-6*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.586626257659714
#define OMIM 0.07573488671206455
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==1003) //rad modified, slow msonic wave, tau=10., P=10., CC=1e2
#define RADIATION
#define KAPPA 10.
#undef SIGMA_RAD 
#define SIGMA_RAD 2.46889e-39
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0000900135
#define DURE (8.99241e-8*RHOFAC)
#define DUIM (7.51995e-10*RHOFAC)
#define DV1RE (5.07163e-6*RHOFAC)
#define DV1IM (7.2926e-7*RHOFAC)
#define DV2RE (6.61826e-6*RHOFAC)
#define DV2IM (1.62925e-6*RHOFAC)
#define B1ZERO 0.0100008
#define B2ZERO 0.0100008
#define DB2RE (-3.23804e-6*RHOFAC)
#define DB2IM (-1.3090e-6*RHOFAC)
#define EEZERO 0.00180027
#define DEERE (-7.12812e-9*RHOFAC)
#define DEEIM (5.99678e-8*RHOFAC)
#define DF1RE (-1.2559e-8*RHOFAC)
#define DF1IM (-1.49156e-9*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.031866
#define OMIM 0.00458207 
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif



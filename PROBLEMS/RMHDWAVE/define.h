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

#define NX 64

#define NUMERO 1101

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



/************** old **************/


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

#if (NUMERO==102) //rad modified, sonic wave, tau=0.1, P=0.01, CC=1e2
#define RADIATION
#define KAPPA 0.1
#undef SIGMA_RAD 
#define SIGMA_RAD 2.46889e-42
#define RHOFAC 1.e-3
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

#if (NUMERO==103) //rad modified, sonic wave, tau=0.01, P=0.01, CC=1e3
#define RADIATION
#define KAPPA 0.01
#undef SIGMA_RAD 
#define SIGMA_RAD 2.46999e-36
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 9.00001e-7
#define DURE (1.49226e-9*RHOFAC)
#define DUIM (7.57753e-11*RHOFAC)
#define DV1RE (9.97736e-7*RHOFAC)
#define DV1IM (2.53158e-8*RHOFAC)
#define DV2RE (0.)
#define DV2IM (0.)
#define B1ZERO 0.0
#define B2ZERO 0.0
#define DB2RE 0.
#define DB2IM 0.
#define EEZERO 1.8e-8
#define DEERE (3.22543e-16*RHOFAC)
#define DEEIM (3.85373e-16*RHOFAC)
#define DF1RE (-3.35928e-14*RHOFAC)
#define DF1IM (7.48005e-14*RHOFAC)
#define DF2RE (0.)
#define DF2IM (0.)
#define OMRE 0.00626896
#define OMIM 0.000159064
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/100.
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

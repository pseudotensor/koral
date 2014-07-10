#define MYCOORDS MINKCOORDS

#define NY 1
#define NZ 1

#define TSTEPLIM .6
#define FULLDISSIPATION

#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define INT_ORDER 1

#define GAMMA (ldouble)(5./3.)
#define MINX 0
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define KAPPAES 0.
#define KK (2.*M_PI)
#define MAGNFIELD
//#define EXPLICIT_RAD_SOURCE

#define RADOUTPUTINFF
#define CALCL1_RMHDWAVE
#define TIMESTEPPING RK2IMEX
#define OUTOUTPUT 1

#define RADIMPCONV 1.e-12
#define U2PCONV 1.e-14

#define NX 128

#define NUMERO 1

#if (NUMERO==1) //sonic wave
#define KAPPA 0.
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.0000152284*RHOFAC)
#define DUIM 0.
#define DV1RE (0.0001*RHOFAC)
#define DV1IM 0.
#define DV2RE 0.
#define DV2IM 0.
#define B1ZERO 0.
#define B2ZERO 0.
#define DB2RE 0.
#define DB2IM 0.
#define OMRE 0.628319
#define OMIM 0.
#define NOUTSTOP 11
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2.*M_PI/OMRE/10.
#endif

#if (NUMERO==10) //fast-magnetosonic wave, no-rad
#define KAPPA 0.
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (0.0000152284*RHOFAC)
#define DUIM 0.
#define DV1RE (0.000160294*RHOFAC)
#define DV1IM 0.
#define DV2RE (-0.0000979087*RHOFAC)
#define DV2IM 0.
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (0.0001623030*RHOFAC)
#define DB2IM 0.
#define OMRE 1.007160
#define OMIM 0.
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==11) //slow-magnetosonic wave, no-rad
#define KAPPA 0.
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.0091370600
#define DURE (0.0000152284*RHOFAC)
#define DUIM 0.
#define DV1RE (0.0000617707*RHOFAC)
#define DV1IM 0.
#define DV2RE (0.000100118*RHOFAC)
#define DV2IM 0.
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (-0.0000625516*RHOFAC)
#define DB2IM 0.
#define OMRE 0.388117
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
#define DURE (0.0000151984*RHOFAC)
#define DUIM (4.81575e-7*RHOFAC)
#define DV1RE (0.000160251*RHOFAC)
#define DV1IM (7.23831e-7*RHOFAC)
#define DV2RE (-0.0000979544*RHOFAC)
#define DV2IM (9.83679e-7*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (0.000162344*RHOFAC)
#define DB2IM (-8.96662e-7*RHOFAC)
#define EEZERO 0.00182741
#define DEERE (1.48421e-9*RHOFAC)
#define DEEIM (6.06322e-8*RHOFAC)
#define DF1RE (-3.95433e-7*RHOFAC)
#define DF1IM (8.51051e-8*RHOFAC)
#define DF2RE (2.3668e-7*RHOFAC)
#define DF2IM (2.11182e-8*RHOFAC)
#define OMRE 1.00689
#define OMIM 0.00454797
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
#define DURE (0.0000150174*RHOFAC)
#define DUIM (1.22299e-6*RHOFAC)
#define DV1RE (0.0000615333*RHOFAC)
#define DV1IM (1.8314e-6*RHOFAC)
#define DV2RE (0.0000989772*RHOFAC)
#define DV2IM (6.54186e-6*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (-0.0000614882*RHOFAC)
#define DB2IM (-5.88315e-6*RHOFAC)
#define EEZERO 0.00182741
#define DEERE (1.96703e-10*RHOFAC)
#define DEEIM (2.18721e-8*RHOFAC)
#define DF1RE (-1.65181e-7*RHOFAC)
#define DF1IM (7.17520e-8*RHOFAC)
#define DF2RE (-2.23679e-7*RHOFAC)
#define DF2IM (-7.43141e-8*RHOFAC)
#define OMRE 0.386625
#define OMIM 0.01150700
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
#define DURE (0.0000117305*RHOFAC)
#define DUIM (1.71290e-6*RHOFAC)
#define DV1RE (0.0002784990*RHOFAC)
#define DV1IM (0.0000523804*RHOFAC)
#define DV2RE (-0.0000281093*RHOFAC)
#define DV2IM (6.25588e-6*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (0.00011017000*RHOFAC)
#define DB2IM (-4.03337e-6*RHOFAC)
#define EEZERO 0.18274100
#define DEERE (0.000207294*RHOFAC)
#define DEEIM (0.000136364*RHOFAC)
#define DF1RE (-0.0000183331*RHOFAC)
#define DF1IM (0.0000363664*RHOFAC)
#define DF2RE (2.67581e-7*RHOFAC)
#define DF2IM (1.24272e-6*RHOFAC)
#define OMRE 1.749860
#define OMIM 0.32911600
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

#if (NUMERO==1102) //rad modified, slow msonic wave, tau=10., P=10., CC=1e1
#define RADIATION
#define KAPPA 10.
#undef SIGMA_RAD 
#define SIGMA_RAD 2.36051e-45
#define RHOFAC 1.e-3
#define RHOZERO 1.
#define DRRE (1.e-3*RHOFAC)
#define DRIM 0.
#define UZERO 0.00913706
#define DURE (9.46189e-6*RHOFAC)
#define DUIM (1.21376e-6*RHOFAC)
#define DV1RE (0.0000834269*RHOFAC)
#define DV1IM (0.0000120829*RHOFAC)
#define DV2RE (0.000113633*RHOFAC)
#define DV2IM (0.000272697*RHOFAC)
#define B1ZERO 0.100759
#define B2ZERO 0.100759
#define DB2RE (-0.0000803823*RHOFAC)
#define DB2IM (-0.000303114*RHOFAC)
#define EEZERO 0.18274111675126906
#define DEERE (0.0000259666*RHOFAC)
#define DEEIM (0.0000967891*RHOFAC)
#define DF1RE (-0.0000198263*RHOFAC)
#define DF1IM (5.4761e-6*RHOFAC)
#define DF2RE (3.66075e-6*RHOFAC)
#define DF2IM (-1.1475e-6*RHOFAC)
#define OMRE 0.524187
#define OMIM 0.075919
#define TMAX 2.*M_PI/OMRE*2.
#define DTOUT1 2*M_PI/OMRE/10.
#define NOUTSTOP 11
#endif

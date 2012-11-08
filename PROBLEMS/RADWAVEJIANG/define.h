#define U2PPREC 1.e-7
#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define MINKOWSKI

#define RK3STEPPING
#define INT_ORDER 1
#define NX 100
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)

#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5

#define DTOUT1 1.0
#define ALLSTEPSOUTPUT 0

#define GAMMA (ldouble)(5./3.)
#define MINX 0
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define EDDINGTON_APR
#define RADIATION

#define KAPPAES 0.

#define NWAVE 1

#undef SIGMA_RAD

#if (NWAVE==1)
#define PP 0.01
#define CC 1.e4
#define KAPPA 0.01
#define DRRE 1.e-3
#define DRIM 0.
#define DVRE 1.27177e-3
#define DVIM 8.15409e-5
#define DPRE 1.61075e-3
#define DPIM 2.07402e-4
#define DERE 1.79137e-8
#define DEIM 8.56498e-9
#define DFRE -1.32035e-6
#define DFIM 3.88814e-6
#define OMRE 7.99077
#define OMIM 0.512336

#define RHO 1.
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)

#elif (NWAVE==4)
#define PP 100.
#define CC 1.e4
#define KAPPA 10.
#define DRRE 1.e-3
#define DRIM 0.
#define DVRE 9.99947e-4
#define DVIM 1.07703e-5
#define DPRE 9.99998e-4
#define DPIM 1.60499e-7
#define DERE -6.60096e-9
#define DEIM 6.41367e-7
#define DFRE -1.00130e-9
#define DFIM 5.35966e-11
#define OMRE 6.28285
#define OMIM 6.76716e-2

#define RHO 1.
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UNIT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)

#endif


#define EFLOOR 1.e-50
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

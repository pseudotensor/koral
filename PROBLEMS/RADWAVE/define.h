#define U2PPREC 1.e-6
#define U2PRADPREC 1.e-7
#define RADFORCEPREC 1.e-5
#define VERBOSE0 0
#define MINKOWSKI

#define RK3STEPPING
#define INT_ORDER 1
#define NY 1
#define NZ 1
#define TSTEPLIM .5
#define INITTSTEPLIM (TSTEPLIM/10.)

#define PERIODIC_XBC
#define COPY_YBC
#define COPY_ZBC
#define FLUXLIMITER 0
#define MINMOD_THETA 2.

#define ALLSTEPSOUTPUT 0

#define GAMMA (ldouble)(5./3.)
#define MINX 0
#define MAXX 1.
#define MINY 0.
#define MAXY 1.
#define MINZ 0.
#define MAXZ 1.

#define EDDINGTON_APR
#define KAPPAES 0.
#define KAPPA 0.
#define NX 100

//#define RADIATION

#undef SIGMA_RAD

#define NWAVE 4


#if (NWAVE==1) //density wave advected with the gas
//#define FLUXDISSIPATIONOFF
#define PP 0.1
#define CC 1.e6
#define VX 1.e-3
#define DTOUT1 (.05/VX)
#define RHO 1.
#define AAA 1.e-5
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#endif

#if (NWAVE==2) //sound wave
//#define FLUXDISSIPATIONOFF
#define PP 0.01
#define CC 1.e6
#define DTOUT1 (.05*CC)
#define VX 0.
#define RHO 1.
#define AAA 1.e-5
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#endif

#if (NWAVE==3) //radiative density wave advected with the gas
#define FLUXDISSIPATIONOFF
#define PP 10.
#define CC 1.e6
#define VX 1.e-2
#define DTOUT1 (.0005/VX)
#define RHO 1.
#define AAA 1.e-5
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#define ERAD calc_LTE_EfromT(TEMP)
#define RADIATION
#undef KAPPAES
#define KAPPAES 10.
#endif


#if (NWAVE==4) //sound wave with radiation
#define FLUXDISSIPATIONOFF
#define PP 1.
#define CC 1.e2
#define DTOUT1 (.005*CC)
#define VX 0.
#define RHO 1.
#define AAA 1.e-1
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHO)
#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#define ERAD calc_LTE_EfromT(TEMP)
#define RADIATION
#undef KAPPA
#define KAPPA 100.
#define ERADFACTOR .5
#define GASFACTOR .5
#endif




#define EFLOOR 1.e-50
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

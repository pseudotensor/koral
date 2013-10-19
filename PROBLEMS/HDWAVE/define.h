#define MYCOORDS MINKCOORDS

#define NY 1
#define NZ 1
#define TSTEPLIM .6
#define TIMESTEPPING RK2
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

#define NX 256

#define NWAVE 2

//#ifdef NWAVE1 //density wave advected with the gas
#if(NWAVE==1)
//#define FLUXDISSIPATIONOFF //switches of the dissipative term in LAXF 
#define CC 1.e6
#define VELX 1.e-3
#define DTOUT1 (.05/VELX)
#define RHOZERO 1.
#define AAA 1.e-5
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHOZERO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHOZERO)
//#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)

#endif

//#ifdef NWAVE2 //hydro sound wave
#if(NWAVE==2)
#define CALCL1_HDWAVE
#define NOUTSTOP 41
#define CC 1.e3 //1/cs
#define DTOUT1 (.05*CC)
#define VELX 0.
#define RHOZERO 1.
#define AAA 1.e-4
#define ERAD 1.
#define KK 2.*Pi
#define UINT (1./CC/CC)*RHOZERO/GAMMA/(GAMMA-1.-1./CC/CC)
#define TEMP calc_PEQ_Tfromurho(UINT,RHOZERO)
//#define SIGMA_RAD (1./4.*PP*(GAMMA-1.)*UINT/TEMP/TEMP/TEMP/TEMP)
#endif


#define EFLOOR 1.e-50
#define UFLOOR 1.e-40
#define RHOFLOOR 1.e-40

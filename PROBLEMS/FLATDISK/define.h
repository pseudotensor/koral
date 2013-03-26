#define RADIATION
//#define EXPLICIT_RAD_SOURCE
//#define IMPLICIT_FF_RAD_SOURCE

#define MASS 10.
#define BHSPIN 0.
#define GAMMAMAXRAD 1000.

#define MYCOORDS SPHCOORDS
//#define LABRADFLUXES

//#define WIDENPRESSURE
#define WIDENPRESSUREPOWER 1.5

#define IMAGETYPE "jpg"
#define OUTVEL VEL4
#define DTOUT1 1.
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 100e10
#define NOUTSTOP 1000.
#define RADOUTPUTINZAMO
#define PRINTYGC_RIGHT

#define MINX 4.
#define MAXX 100.

#define NX 40
#define NY 40
#define NZ 1


#define MINY (0.01*Pi/4.)
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC
//#define COPY_XBC
//#define COPY_YBC
//#define COPY_ZBC

#define GAMMA (4./3.)
#define KKK 1.e-4
#define ELL 4.5

//#define RHOATMMIN  rhoCGS2GU(1.e-4)
#define RHOATMMIN  1.e-2
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e11,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(1.e9))

#define INT_ORDER 1
#define RK2_STEPPING
#define TSTEPLIM 1.
#define FLUXLIMITER 0
#define MINMOD_THETA 1.

#define NODONUT 0
#define INFLOWING 0

#define RHOFLOOR 1.e-40
#define UFLOOR 1.e-40
#define EFLOOR 1.e-40

int find_globalindex(double r, double th, double ph, int gi[3]);

//************************************/
//SPH input parameters
/************************************/

#define SPHRHOCUT 1.e-40
#define SPHRHOCUTBFIELD 1.e-12
#define SPHPRECUTBFIELD 1.e-21
#define SPHTHBCUT 0.002 //max theta from eq. plane for initial B field
#define SPHBFIELDRHO

#define SPHSMEARX 0
#define SPHSMEARY 1
#define SPHSMEARZ 0

#define SPHMASSNORM 1.9891e32 //mass of the star in cgs




#define FULLPHI

//************************************/
//general
/************************************/
#define BHDISK_PROBLEMTYPE

/************************************/
//restart
/************************************/
#define RESTART
#define RESTARTGENERALINDICES
#define RESTARTNUM -1

/************************************/
//radiation choices
/************************************/
#define RADIATION
#define RADIMPLICITTHRESHOLD 1.e0
#define MAXRADIMPDAMPING 1.e-6
//#define SKIPRADSOURCE
#define BALANCEENTROPYWITHRADIATION
#define COMPTONIZATION
#define ALLOWRADCEILINGINIMPLICIT
//#define RADIMPLICITFIXVEL
                                                                                                                                                
#define RADIMPCONV 1.e-8
#define RADIMPEPS 1.e-6
#define RADIMPMAXITER 50
#define RADIMPCONVREL 1.e-6
#define RADIMPCONVRELERR 1.e-1
#define RADIMPCONVRELENTR 1.e-4
#define RADIMPCONVRELENTRERR .999

//#define BASICRADIMPLICIT
//#define RADIMPSTARTWITHEXP
//#define ALLOWFORENTRINF4DPRIM

//#define U2P_EQS U2P_EQS_JONS
//#define U2P_SOLVER U2P_SOLVER_WP


#define U2PCONV 1.e-12
#define RADIMPMAXITER 50


/************************************/
//magnetic choices
/************************************/
//if we want a magnetic field, uncomment MAGNFIELD
//#define MAGNFIELD
#define MAXBETA 0.1
#define VECPOTGIVEN
#define MPI4CORNERS
#define GDETIN 1

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#ifdef RADIATION
#define TIMESTEPPING RK2IMEX
#else
#define TIMESTEPPING RK2HEUN
#endif
#define TSTEPLIM .5
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
#define DOFIXUPS 1
#define DOU2PRADFIXUPS 0
#define DOU2PMHDFIXUPS 1
#define DORADIMPFIXUPS 1


/************************************/
//viscosity choices
/************************************/
#ifdef RADIATION
#define RADVISCOSITY SHEARVISCOSITY
#endif
#define RADVISCMFPSPH
#define RADVISCNUDAMP
#define ALPHARADVISC .1
#define RADVISCMAXVELDAMP
#define MAXRADVISCVEL 1.

/************************************/
//rmhd floors
/************************************/
#define CORRECT_POLARAXIS
#define NCCORRECTPOLAR 1
#define UURHORATIOMIN 1.e-9 // 1K: u/rho = 7.259162e+12
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-20
#define EERHORATIOMAX 1.e20
#define EEUURATIOMIN 1.e-20
#define EEUURATIOMAX 1.e20
#define B2UURATIOMIN 0.
#define B2UURATIOMAX 100000.
#define B2RHORATIOMIN 0.
#define B2RHORATIOMAX 100.
#define GAMMAMAXRAD 50.
#define GAMMAMAXHD 50.

/************************************/
//blackhole
/************************************/
#define MASS 1.e5//1.e1
#define BHSPIN 0.0

/************************************/
//coordinates / resolution
/************************************/
#define myMKS2COORDS
#define METRICAXISYMMETRIC
#define ROUT 1000.
#define RMIN 1.7

#ifdef myMKS1COORDS //modified Kerr-Shild
#define MKSR0 0.
#define MYCOORDS MKS1COORDS
#define MINX (log(3.6-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.0025*Pi/2.)
#define MAXY (Pi-0.0025*Pi/2.)
#endif

#ifdef myMKS2COORDS //modified Kerr-Shild with more cells towards the eq.plane
#define MKSR0 0.
#define MKSH0 0.9 //makes cells smaller towards equatorial plane
#define MYCOORDS MKS2COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(ROUT-MKSR0))
#define MINY (0.01)
#define MAXY (1.-MINY)
#endif


#ifdef myMKS3COORDS //modified Kerr-Shild further from axis
#define MKSMY1 0.001
#define MKSMY2 0.4
#define MKSH0 0.85
#define MKSR0 0.
#define MKSMP0 1.5
#define METRICNUMERIC
#define MYCOORDS MKS3COORDS
#define MINX (log(RMIN-MKSR0))
#define MAXX (log(100.-MKSR0))
#define MINY 0.
#define MAXY 1.
#endif

//total resolution
#define TNX 128// 224//128//64
#define TNY 64//128//64//32
#define TNZ 32//96//64//32
//number of tiles
#define NTX 2
#define NTY 4
#define NTZ 4

//#define HALFTHETA //symmetry wrt eq. plane?
#ifdef HALFTHETA
#undef TNY
#undef MAXY
#define TNY 50
#ifdef myMKS1COORDS 
#define MAXY (Pi/2.)
#endif
#ifdef myMKS2COORDS 
#define MAXY (.5)
#endif
#endif



#define PHIWEDGE (2.*M_PI)
#define MINZ (-PHIWEDGE/2.)
#define MAXZ (PHIWEDGE/2.)

#define PERTMAGN 1.e-2
#define SPECIFIC_BC
#define PERIODIC_ZBC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define OUTVEL VEL4
#define ALLSTEPSOUTPUT 0
#define RADOUTPUTINZAMO
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000


#define PRINTXGC_RIGHT
#define PRINTZONEMORE
#define SCAOUTPUT 1
#define SILOOUTPUT 1
#define RADOUTPUT 0
#define AVGOUTPUT 1
#define COORDOUTPUT 2

#define DTOUT1 .01
#define DTOUT2 1000.

/************************************/
//common physics / torus / atmosphere
/************************************/
#define GAMMA (5./3.)

//parameters for the magnetic field
#define VERTBTIME 1000.
#define MAGNOMEGA 0.//(2.*M_PI/1000.)//0.
#define MAGBETA 0.1

//atmosphere
#define RHOATMMIN  1.e-20
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  1.e-35//(calc_LTE_EfromT(3.e6)/10/100000.)


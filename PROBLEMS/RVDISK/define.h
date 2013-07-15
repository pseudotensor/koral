/************************************/
//restart
/************************************/
//#define RESTART
#define RESTARTNUM 1

/************************************/
//radiation choices
/************************************/
#define RADIATION
//#define SKIPRADSOURCE
//#define SKIPRADWAVESPEEDLIMITER
#define ALLOW_EXPLICIT_RAD_SOURCE 0

/************************************/
//hydro choices
/************************************/
#define ALLOWENTROPYU2P 1
#define DOFIXUPS 0

/************************************/
//reconstruction / Courant
/************************************/
#define INT_ORDER 1
#define TIMESTEPPING RK2K1K2
#define TSTEPLIM .6
#define FLUXLIMITER 0
#define MINMOD_THETA 1.5
//#define FLUXMETHOD HLL_FLUX
//#define WAVESPEEDSATFACES 
#define GDETIN 0

/************************************/
//viscosity choices
/************************************/
#define HDVISCOSITY SIMPLEVISCOSITY
#define ALPHATOTALPRESSURE
#define RMINVISC 2.
#define RADVISCOSITY SHEARVISCOSITY
#define TAUSUPPRESSPARAM 100. //the larger the less prad
#define ALPHARADVISC 1.
//#define ENFORCERADWAVESPEEDS

/************************************/
//rhd floors
/************************************/
#define UURHORATIOMIN 1.e-15
#define UURHORATIOMAX 1.e2
#define EERHORATIOMIN 1.e-15
#define EERHORATIOMAX 1.e6
#define EEUURATIOMIN 1.e-15
#define EEUURATIOMAX 1.e6
#define ERADLIMIT 1.e-50
#define RHOFLOOR 1.e-50
#define GAMMAMAXRAD 50.

/************************************/
//blackhole
/************************************/
#define MASS 10.

/************************************/
//coordinates / resolution
/************************************/
#define MYCOORDS MKS1COORDS
#define MAXY Pi/2.
#define MINZ -1.
#define MAXZ 1.
#define SPECIFIC_BC

/************************************/
//output
/************************************/
#define OUTCOORDS KERRCOORDS                                                                    
#define RADOUTPUTINZAMO
//#define PRINTINSIDEBH
//#define PRINTXGC_LEFT
//#define PRINTXGC_RIGHT
#define ALLSTEPSOUTPUT 0
#define NSTEPSTOP 1.e10
#define NOUTSTOP 5000
#define CGSOUTPUT

/************************************/
//common physics 
/************************************/
#define GAMMA (4./3.)

/************************************/
//model choice
/************************************/
#define NDISK 990

/************************************/
#if (NDISK==990) //mdot = 10, r=15/50 alpha=0.1
/************************************/

#undef RMINVISC
#define RMINVISC 0.
#undef RADIATION
#undef ALPHATOTALPRESSURE
#undef HDVISCOSITY 
#define HDVISCOSITY SHEARVISCOSITY
#undef RADVISCOSITY 
#define RADVISCOSITY NOVISCOSITY
#define ALPHAHDVISC 0.1
#define BHSPIN 0.0
#define MKS1R0 -10.
#define MDOTOUT 50.
#define RKEP 30.
#define ROUT 50.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))

#undef MYCOORDS
#define MYCOORDS KERRCOORDS
#undef MINX
#define MINX 3.
#undef MAXX
#define MAXX ROUT

#define MINY (0.05*Pi/2.)
#define NX 60
#define NY 60
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif


/************************************/
#if (NDISK==100) //mdot = 10, r=15/50 alpha=0.1
/************************************/
#define ALPHAHDVISC .1
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 50.
#define RKEP 15.
#define ROUT 50.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 50
#define NY 30
#define NZ 1
#define DTOUT1 5.e-1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==101) //mdot = 10, r=50/100, alpha=0.1
/************************************/
#define ALPHAHDVISC .1
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 50.
#define RKEP 50.
#define ROUT 100.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 70
#define NY 40
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==102) //mdot = 10, r=50/100, alpha=0.02
/************************************/
#define ALPHAHDVISC .02
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 50.
#define RKEP 50.
#define ROUT 100.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 70
#define NY 40
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==103) //mdot = 1, r=50/100, alpha=0.1
/************************************/
#define ALPHAHDVISC .1
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 10.
#define RKEP 50.
#define ROUT 100.
#define INJTHETA ((1.-0.2)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 70
#define NY 40
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==104) //mdot = 50, r=50/100, alpha=0.1
/************************************/
#define ALPHAHDVISC .1
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 200.
#define RKEP 50.
#define ROUT 100.
#define INJTHETA ((1.-0.4)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 70
#define NY 40
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==105) //mdot = 10, r=50/100, alpha=0.1, spin=0.9
/************************************/
#define ALPHAHDVISC .1
#define BHSPIN 0.9
#define MKS1R0 -3.
#define MDOTOUT 50.
#define RKEP 50.
#define ROUT 100.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.6-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 70
#define NY 40
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==106) //mdot = 10, r=100/200, alpha=0.1
/************************************/
#define ALPHAHDVISC .1
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 50.
#define RKEP 100.
#define ROUT 200.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 80
#define NY 50
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif

/************************************/
#if (NDISK==107) //mdot = 10, r=50/100, alpha=0.1, TAUSUP=10
/************************************/
#undef TAUSUPPRESSPARAM
#define TAUSUPPRESSPARAM 10. //the larger the less prad
#define ALPHAHDVISC .1
#define BHSPIN 0.0
#define MKS1R0 -4.
#define MDOTOUT 50.
#define RKEP 50.
#define ROUT 100.
#define INJTHETA ((1.-0.3)*M_PI/2.)
#define MINX (log(.8-MKS1R0))
#define MAXX (log(ROUT-MKS1R0))
#define MINY (0.01*Pi/2.)
#define NX 70
#define NY 40
#define NZ 1
#define DTOUT1 5.e1
#define RHOATMMIN  1.e-25
#define UINTATMMIN  (calc_PEQ_ufromTrho(1.e10,RHOATMMIN))
#define ERADATMMIN  (calc_LTE_EfromT(3.e6)/10.)
#endif


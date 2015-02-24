
/*********************/
/*********************/
/*********************/
/*********************/
//important choices **/
/*********************/
/*********************/
/*********************/
/*********************/

#ifndef VELPRIM
#define VELPRIM VELR
#endif

#ifndef VELPRIMRAD
#define VELPRIMRAD VELR
#endif

#ifndef DORADFIXUPS
#define DORADFIXUPS 1
#endif

#ifndef DOFIXUPS
#define DOFIXUPS 1
#endif

#ifndef ALLOWENTROPYU2P
#define ALLOWENTROPYU2P 1
#endif

#ifndef ALLOWCOLDU2P
#define ALLOWCOLDU2P 0
#endif

#define SMALL 1.e-50 //small number 
#define BIG (1./SMALL) //big number

//min uint over rho 
#ifndef UURHORATIOMIN
#define UURHORATIOMIN 1.e-50
#endif

//uint over rho for u2p_cold
#ifndef UURHORATIOU2PCOLD
#define UURHORATIOU2PCOLD 1.e-10
#endif

//max uint over rho
#ifndef UURHORATIOMAX 
#define UURHORATIOMAX 1.e2
#endif

//min Erad over rho
#ifndef EERHORATIOMIN
#define EERHORATIOMIN 1.e-30
#endif

//max Erad over rho
#ifndef EERHORATIOMAX 
#define EERHORATIOMAX 1.e30
#endif

//min Erad over uint
#ifndef EEUURATIOMIN
#define EEUURATIOMIN 1.e-30
#endif

//max Erad over uint
#ifndef EEUURATIOMAX 
#define EEUURATIOMAX 1.e30
#endif

//min absolute Erad
#ifndef ERADFLOOR
#define ERADFLOOR (10.*SMALL)
#endif

//min B^2 over uint
#ifndef B2UURATIOMIN 
#define B2UURATIOMIN 0.
#endif

//max B^2 over uint
#ifndef B2UURATIOMAX 
#define B2UURATIOMAX 100.
#endif

//max B^2 over Ehat
#ifndef B2EERATIOMAX
#define B2EERATIOMAX 100.
#endif

//min B^2 over rho 
#ifndef B2RHORATIOMIN 
#define B2RHORATIOMIN 0.
#endif

//max B^2 over rho 
#ifndef B2RHORATIOMAX 
#define B2RHORATIOMAX 100.
#endif

/*********************/
/*********************/
/*********************/
/*********************/
//passive definitions
/*********************/
/*********************/
/*********************/
/*********************/

//default number of radiative fluids
//redefined later if single fluid
#ifndef NRF
#define NRF 6
#endif

#ifndef MFSKEW
#define MFSKEW 20.
#endif

#ifndef MFMINVEL
#define MFMINVEL 1.e-3
#endif

#ifndef MFFRACSCALE
#define MFFRACSCALE 1.
#endif

#ifndef MFREDISTRIBUTEMETHOD
#define MFREDISTRIBUTEMETHOD 2
#endif

#ifndef MYCOORDS2
#define MYCOORDS2 MYCOORDS
#endif

#ifndef OUTCOORDS
#define OUTCOORDS MYCOORDS
#endif

#ifndef BHSPIN
#define BHSPIN 0.
#endif

#ifndef NOUTSTOP
#define NOUTSTOP 1e50 //max n of outputs
#endif

#ifndef NSTEPSTOP
#define NSTEPSTOP 1e50 //max n of steps
#endif

#ifndef TMAX
#define TMAX 1.e50  //max time
#endif

#ifndef VERBOSE0
#define VERBOSE0 0 //verbose level
#endif

#ifndef EXPLICIT_LAB_RAD_SOURCE
#ifndef EXPLICIT_SUBSTEP_RAD_SOURCE
#ifndef IMPLICIT_LAB_RAD_SOURCE
#define IMPLICIT_LAB_RAD_SOURCE
#endif
#endif
#endif

#ifndef ALLOW_EXPLICIT_RAD_SOURCE
#define ALLOW_EXPLICIT_RAD_SOURCE 0 //whether to allow reducing implicit_lab to explicit
#endif

#ifndef INT_ORDER
#define INT_ORDER 1 //reconstruction order
#endif

#ifndef TIMESTEPPING
#define TIMESTEPPING RK2IMEX //time stepping
#endif

#ifndef NG
#if (INT_ORDER==0)
#define NG 2 //number of ghost cells
#endif
#if (INT_ORDER==1)
#define NG 2 
#endif
#if (INT_ORDER==2)
#define NG 3
#endif
#if (INT_ORDER==4)
#define NG 4
#endif
#endif



#ifndef NUM_INPUTARG
#define NUM_INPUTARG 0 //number of input arguments in the command line
#endif

//number of hydro variables
#ifndef TRACER
#define NVHD (6)
#else
#define NVHD (6+1)
#endif

//number of magneto-hydro variables
#ifdef MAGNFIELD
#define NVMHD (NVHD+3)
#else
#define NVMHD (NVHD)
#endif

//number of total variables
#ifdef RADIATION

//number of radiative quantities per fluid
#ifndef NCOMPTONIZATION
#define NRADVAR 4
#else
#define NRADVAR 5
#endif

#ifndef MULTIRADFLUID
#undef NRF
#define NRF 1
#endif

#define NV (NVMHD+NRADVAR*NRF)

#else //no RADIATION

#define NV (NVMHD)
#define NRADVAR 4 //not used
#endif

#ifndef GAMMA
#define GAMMA (5./3.) //gamma
#endif

#ifndef MASS
#define MASS 1./MSUNCM //default mass of the BH used to calibrate radiation constant, Solar mass units
#endif

#ifndef U2PRADPREC
#define U2PRADPREC 1.e-5
#endif

#ifndef NSCALARS
#define NSCALARS 13
#endif

#ifndef NAVGVARS
//#ifdef BHDISK_PROBLEMTYPE
#ifdef RADIATION
#define NAVGVARS (151+3*NV) //added to existing NV 
#else
#define NAVGVARS (110+3*NV)
#endif
//#else
//#define NAVGVARS (0)
//#endif
#endif

#ifndef NRADPROFILES
#define NRADPROFILES (47-1)
#endif

#ifndef NTHPROFILES
#define NTHPROFILES 6
#endif

#ifndef NANARELRADPROFILES
#define NANARELRADPROFILES 4
#endif

#ifndef OUTVEL
#define OUTVEL VEL3
#endif

#ifndef RHOATMMIN
#define RHOATMMIN 1.
#endif

#ifndef ERADATMMIN
#define ERADATMMIN 1.
#endif

#ifndef UINTATMMIN
#define UINTATMMIN 1.e-2
#endif

#ifndef EEFLOOR
#define EEFLOOR 1.e-50
#endif

#ifndef RHOFLOOR
#define RHOFLOOR 1.e-50
#endif

#ifndef UUFLOOR
#define UUFLOOR 1.e-50
#endif

#ifndef GAMMAMAXHD
#define GAMMAMAXHD 100.
#endif

#ifndef GAMMAMAXRAD
#define GAMMAMAXRAD 100.
#endif

#ifndef MASSCM
#define MASSCM (MASS*MSUNCM) //mass in cm
#endif

#ifndef MAXEXPLICITSUBSTEPCHANGE
#define MAXEXPLICITSUBSTEPCHANGE 1.e-2
#endif

#ifndef IMAGETYPE
#define IMAGETYPE "gif"
#endif

#ifndef GAMMASMALLLIMIT
#define GAMMASMALLLIMIT (1.0-1E-10) // at what point above which assume gamma^2=1.0
#endif

#ifndef FLUXMETHOD
#define FLUXMETHOD LAXF_FLUX
#endif

#ifndef MFWEDGESTYPE
#define MFWEDGESTYPE 1
#endif

#ifndef GDETIN
#define GDETIN 1 //whether to include metric determinant into the fluxes; must be on for magnetic fields 'cos flux_ct assumes that
#endif

#ifndef MODYFIKUJKRZYSIE
#if (GDETIN==1)
#define MODYFIKUJKRZYSIE 1
#else
#define MODYFIKUJKRZYSIE 0
#endif
#endif

#ifndef HDVISCOSITY
#define HDVISCOSITY NOVISCOSITY
#endif

#ifndef RADVISCOSITY
#define RADVISCOSITY NOVISCOSITY
#endif

#ifndef MAXRADVISCVEL
#define MAXRADVISCVEL 1.
#endif

#ifndef SHUFFLELOOPS
#define SHUFFLELOOPS 0
#endif

#ifndef MKSR0
#define MKSR0 0.
#endif

//whether to check if the advection operator keeps entropy increasing,
//if not invert with the independently evolved entropy
#ifndef VERIFYENTROPYAFTERADVECTION
#define VERIFYENTROPYAFTERADVECTION 0
#endif

#ifndef OUTOUTPUT
#define OUTOUTPUT 0
#endif

#ifndef SCAOUTPUT
#define SCAOUTPUT 0
#endif

#ifndef RADOUTPUT
#define RADOUTPUT 0
#endif

#ifndef AVGOUTPUT
#define AVGOUTPUT 0
#endif

#ifndef SILOOUTPUT
#define SILOOUTPUT 0
#endif

#ifndef GRIDOUTPUT
#define GRIDOUTPUT 0
#endif

#ifndef SIMOUTPUT
#define SIMOUTPUT 0
#endif

#ifndef DTOUT2
#define DTOUT2 DTOUT1
#endif

#ifndef B2RHOFLOORFRAME
#define B2RHOFLOORFRAME ZAMOFRAME
#endif

#ifndef NCCORRECTPOLAR
#define NCCORRECTPOLAR 2
#endif

#ifndef NCELLSINSIDEHORIZON
#define NCELLSINSIDEHORIZON 6
#endif

#ifndef ALLSTEPSOUTPUT
#define ALLSTEPSOUTPUT 0
#endif

#ifndef U2PCONV
#define U2PCONV 1.e-10
#endif

#ifndef RADIMPCONV
#define RADIMPCONV 1.e-10
#endif

#ifndef RADIMPEPS
#define RADIMPEPS 1.e-6
#endif

#ifndef RADIMPMAXITER
#define RADIMPMAXITER 50
#endif

#define NUMEPSILON DBL_EPSILON

/*********************/
/*********************/
/*********************/
/*********************/
/***** wrappers ******/
/*********************/
/*********************/
/*********************/
/*********************/

#define GAMMAM1 (GAMMA-1.) //gamma - 1

#define IGAMMAR (GAMMAM1/GAMMA)

#define LCM (MASSCM) //unit of length in cm

#define TSEC (MASSCM/CCC) //unit of time in seconds

#define GMC2CM (MASSCM) //gravitational radius in cm

/*********************/
/*********************/
/*********************/
/*********************/
/***** mpi-spec ******/
/*********************/
/*********************/
/*********************/
/*********************/

#if !defined(MPI) && !defined(OMP)
#undef NTX
#undef NTY
#undef NTZ
#define NTX 1 //number of tiles in X
#define NTY 1
#define NTZ 1
#endif

#ifndef NTX
#define NTX 1 //number of tiles in X
#endif

#ifndef NTY
#define NTY 1
#endif

#ifndef NTZ
#define NTZ 1
#endif



#ifdef MPI
#ifndef NX
#define NX (TNX/NTX)
#endif 

#ifndef NY
#define NY (TNY/NTY)
#endif
 
#ifndef NZ
#define NZ (TNZ/NTZ)
#endif 

#else //OMP or single core

#ifndef NX
#define NX (TNX)
#endif 

#ifndef NY
#define NY (TNY)
#endif
 
#ifndef NZ
#define NZ (TNZ)
#endif 
#endif


#ifndef RADCLOSURE
#define RADCLOSURE M1CLOSURE
#endif

#ifndef RADCLOSURECOORDS
#define RADCLOSURECOORDS MYCOORDS
#endif


#ifndef MPIMSGBUFSIZE
#define MPIMSGBUFSIZE 12
#endif

#ifndef MPI4CORNERS //required by MAGNFIELD and VETCLOSURE
#ifdef MAGNFIELD
#define MPI4CORNERS
#endif
#endif

#ifndef MPI4CORNERS //required by MAGNFIELD and VETCLOSURE
#ifdef RADIATION //required by MAGNFIELD and VETCLOSURE
#if (RADCLOSURE==VETCLOSURE)
#define MPI4CORNERS
#endif
#endif
#endif

#ifdef MPI4CORNERS
#undef MPIMSGBUFSIZE
#if (TNX>1 && TNY>1 && TNY>1) //3d
#define MPIMSGBUFSIZE 52
#else //2d
#define MPIMSGBUFSIZE 20
#endif
#endif


#ifndef MPI
#ifndef OUTPUTPERCORE
#define OUTPUTPERCORE
#endif
#endif

//ZERO short solver

#ifndef NUMANGLES
#ifdef EVOLVEINTENSITIES
#define NUMANGLES 80 		//Must be equal to number of lines in angle grid file
#else
#define NUMANGLES 1
#endif
#endif

//to limit memory footprint

#if (RADCLOSURE==VETCLOSURE)
#define SXVET SX
#define SYVET SY
#define SZVET SZ
#else
#define SXVET 1
#define SYVET 1
#define SZVET 1
#endif

#ifndef NUMDUALANGLES
#define NUMDUALANGLES 156 	//Must be equal to number of lines in dual angle grid file
#endif

#ifndef NUMMSTEPLEVELS
#define NUMMSTEPLEVELS 1
#endif

#ifndef RADIMPLICITTHRESHOLD
#define RADIMPLICITTHRESHOLD 1.e-2
#endif

#ifndef SUBZONES_NSTEPSTEP
#define SUBZONES_NSTEPSTEP 100
#endif

#ifndef SILOCOORDS
#define SILOCOORDS MINKCOORDS
#endif

#ifndef VETFLUXCOSACCEPT
#define VETFLUXCOSACCEPT 1.
#endif

#ifndef VETFEACCEPT
#define VETFEACCEPT 0.
#endif

#ifndef U2P_EQS
#define U2P_EQS U2P_EQS_NOBLE
#endif


#ifndef U2P_SOLVER
#define U2P_SOLVER U2P_SOLVER_W
#endif


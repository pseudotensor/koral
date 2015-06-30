//restart - related 
//#define RESTART
#define RESTARTNUM -1 //-1 means from the last dump

#define MYCOORDS MINKCOORDS //coordinate system, Cartesian
//SPHCOORDS
//CYLCOORDS
//BLCOORDS (Kerr)
//KSCOORDs (Kerr-Shild)
//MKS2COORDS (KS with more cells at the equatorial plane and logarithmic cells radius)
//...

#define MASS 1./MSUNCM //mass of the BH in Msun, this one gives 1Rg=1cm

#define MINX -1. //size of the box units of Rg 
#define MAXX 1.
#define MINY -1. //in radians for sph-like coordinates, still in Rg for cylindrical
#define MAXY 1.
#define MINZ -1. //in rad. for sph and cyl
#define MAXZ 1.

#define TNX 60 //total number of cells in x
#define TNY 60
#define TNZ 1

#define NTX 2 //number of tiles in x
#define NTY 2
#define NTZ 1

//don't go below 3 in one dimension per tile unless 2D

#define TIMESTEPPING RK2HEUN //time stepping algorithm
//RK2IMEX - better for radiation

#define MAGNFIELD //whether we evolve magn. fields
#define VECPOTGIVEN //whether you give vector potenital initially
#define MAXBETA .1

#define INT_ORDER 1 //integration order 
//1 - linear interpolation, most stable
//2 - PPM (parabolic piece wise), less stable 
//4 - MP5

#define TSTEPLIM .5 //Courant-limiter, not to exceed .7, the more, the less stable

#define SPECIFIC_BC //whether you specify some boundary condition manually
//#define PERIODIC_XBC //periodic, handled automatically
//#define PERIODIC_YBC
#define PERIODIC_ZBC

#define FLUXLIMITER 0 //reconstruction MINMOD parameters 
#define MINMOD_THETA 1. //1. - most diffusive,  and stable, 2. - least diffusive

#define DTOUT1 .5 //dt for dumping snapshots in units GM/c3
#define NOUTSTOP 100 //max n of dumps
#define NSTEPSTOP 1e50 //max n of time-steps
#define TMAX 1.e50  //max time


//output choices 
#define SILOOUTPUT 1 //SILO file, .silo, for Visit
#define SCAOUTPUT 0 //scalars vs t
#define RADOUTPUT 0 //radial profiles vs t
#define AVGOUTPUT 0 //avg quantities output
#define COORDOUTPUT 0 //to dump ascii file with coordinates
#define SIMOUTPUT 0 //simple output to ascii file, to use in python

//problem parameters:
#define PL1X1 .4
#define PL1X2 .6
#define PL1Y -.2
#define PL1VX -.05
#define PL1RHO 1.
#define PL1TEMP 1.e8

#define PL2X1 -.6
#define PL2X2 -.4
#define PL2Y .2
#define PL2VX .05
#define PL2RHO 1.
#define PL2TEMP 1.e8

#define AMBRHO 1.e-5
#define AMBTEMP 1.e10


#define U_AMB 1.e-7
#define VELWIND 0.1


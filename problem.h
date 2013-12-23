//KORAL - problem.h
//choice of the problem plus some definitions

//available problems:

//1 RADBEAM2D - beam of light 
//2 RADINFALL - radial inflow
//3 DONUT - 2d Polish donut
//4 GEODESICINFALL - geodesic infall with blobs or not
//5 EDDINFALL - infall with flux from inside
//6 RADTUBE- radiative shock tubes as in Farris et al 09
//7 BONDI - like in Fragile's paper
//8 HDTUBE - relativistic shock tube
//9 HDTUBE2D - in 2d
//10 RADPULSE - radiative blob spreading around
//11 RADSHADOW - radiative shadow
//12 RADATM - atmosphere enlighted
//13 DONUTOSC - 2d Polish donut oscillating
//14 RADWAVEBC - 1d linear rad wave imposed on boundary
//15 RADWAVE - 1d linear rad wave with periodic BC
//16 RADPULSE3D - radiative blob spreading around
//17 RADDBLSHADOW - radiative shadow with two beams inclined
//18 ATMSTATIC - hydro atmosphere 
//19 RADBEAM2DKS - beam of light in KS coordinates
//20 ATMKS - radial atmosphere infalling in KS
//21 DONUTKS - 2d Polish donut in KS
//22 DONUTMKS1 - 2d Polish donut in MKS1
//23 ATMMKS1 - radial atmosphere infalling in MKS1
//24 RADBEAMFLAT - beam of light in Cartesian 
//25 RDONUT - 2d radiative Polish donut in KS
//26 RADBEAM2DKSVERT - 2d radiative beam in r,theta plane
//27 RADFLATNESS - flat but with non-zero four-force
//28 BOWSHOCK - bow shock hydro test
//29 RADWALL - flat with wall
//30 RADNT - emission from midplane
//31 FLATDISK - emission from flat disk
//32 CYLBEAM - beam towards the axis in cylindrical
//33 RADDOT - radiating dots
//34 MFPULSE - multi fluid pulse
//35 MFBEAMS - multi fluid colliding beams
//36 MFDOTS - multi fluid radiating dots
//37 MFCYLBEAM - beam towards the axis in cylindrical with multifluids
//40 CYLBEAMCART - similar to discrete CYLBEAM but in cartesian 
//41 FLATDOT - dot which may be bigger than a dot
//42 RVDONUT - radiative and viscous dougnut
//43 RVDONUTIN - radiative and viscous dougnut inflowing
//44 RADNTCYL - emission from midplane in cylindrical
//45 MFRADNTCYL - multi fluid emission from midplane in cylindrical
//46 MFRADNTSPH - multi fluid emission from midplane in spherical
//47 MFDONUT - inflowing donut with radiation/viscosity and multi-fluids
//48 RVDTEST - testing the entropy near the axis
//49 RVRING - 1d donut ring
//50 HDWAVE - hydro wave for testing entropy
//51 HUBBLE - hubble-type test as in the WHAM paper
//52 SPHFLAT - spherical flat to test Christoffels
//53 1DDONUT - 1d donut equatorial plane structure
//54 RVDISK - viscous rad. disk damped from larger radius
//55 G2STATIC - cloud through static atmosphere
//56 EDDSPH - 1d Edd Sphere
//57 LIMOTORUS - torus from Bob's/Akshay's paper
//58 LUKE - beam hitting a net gas flow
//59 BBBLOB - rad blobs
//60 JONTESTS - testing jon's fails.dat
//61 MAGTESTS - simple mag tests
//62 MAGTUBES - magnetic tubes
//63 ORSZAG - Orszag-Tang vortex
//64 MAGDONUT - MHD donut with poloidal magnetic fields
//65 MRTORUS - RMHD torus with radiation
//66 KTORUS - RMHD Newtonian torus with radiation
//67 LRTORUS - RMHD limo torus
//68 RMHDWAVE - radiation modified linear magnetosonic waves
//69 INJDISK - magnetized gas damped from outer boundary
//70 M1SHOCK - light beam penetrating absorbing medium, Monreal & Frank 08

#define PROBLEM 59

#if(PROBLEM==70)

#define PR_DEFINE "PROBLEMS/M1SHOCK/define.h"
#define PR_BC "PROBLEMS/M1SHOCK/bc.c"
#define PR_INIT "PROBLEMS/M1SHOCK/init.c"
#define PR_KAPPA "PROBLEMS/M1SHOCK/kappa.c"
#define PR_KAPPAES "PROBLEMS/M1SHOCK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/M1SHOCK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/M1SHOCK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/M1SHOCK/dump.c"
#define PR_TOOLS "PROBLEMS/M1SHOCK/tools.c"

#endif

#if(PROBLEM==69)

#define PR_DEFINE "PROBLEMS/INJDISK/define.h"
#define PR_BC "PROBLEMS/INJDISK/bc.c"
#define PR_INIT "PROBLEMS/INJDISK/init.c"
#define PR_KAPPA "PROBLEMS/INJDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/INJDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/INJDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/INJDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/INJDISK/dump.c"
#define PR_TOOLS "PROBLEMS/INJDISK/tools.c"
#define PR_POSTINIT "PROBLEMS/INJDISK/postinit.c"

#endif

#if(PROBLEM==68)

#define PR_DEFINE "PROBLEMS/RMHDWAVE/define.h"
#define PR_BC "PROBLEMS/RMHDWAVE/bc.c"
#define PR_INIT "PROBLEMS/RMHDWAVE/init.c"
#define PR_KAPPA "PROBLEMS/RMHDWAVE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RMHDWAVE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RMHDWAVE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RMHDWAVE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RMHDWAVE/dump.c"
#define PR_TOOLS "PROBLEMS/RMHDWAVE/tools.c"

#endif

#if(PROBLEM==67)

#define PR_DEFINE "PROBLEMS/LRTORUS/define.h"
#define PR_BC "PROBLEMS/LRTORUS/bc.c"
#define PR_INIT "PROBLEMS/LRTORUS/init.c"
#define PR_KAPPA "PROBLEMS/LRTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/LRTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LRTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LRTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LRTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/LRTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/LRTORUS/postinit.c"

#endif

#if(PROBLEM==66)

#define PR_DEFINE "PROBLEMS/KTORUS/define.h"
#define PR_BC "PROBLEMS/KTORUS/bc.c"
#define PR_INIT "PROBLEMS/KTORUS/init.c"
#define PR_KAPPA "PROBLEMS/KTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/KTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/KTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/KTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/KTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/KTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/KTORUS/postinit.c"

#endif

#if(PROBLEM==65)

#define PR_DEFINE "PROBLEMS/MRTORUS/define.h"
#define PR_BC "PROBLEMS/MRTORUS/bc.c"
#define PR_INIT "PROBLEMS/MRTORUS/init.c"
#define PR_KAPPA "PROBLEMS/MRTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MRTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MRTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MRTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MRTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/MRTORUS/tools.c"
#define PR_POSTINIT "PROBLEMS/MRTORUS/postinit.c"

#endif

#if(PROBLEM==64)

#define PR_DEFINE "PROBLEMS/MAGDONUT/define.h"
#define PR_BC "PROBLEMS/MAGDONUT/bc.c"
#define PR_INIT "PROBLEMS/MAGDONUT/init.c"
#define PR_KAPPA "PROBLEMS/MAGDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/MAGDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MAGDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MAGDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MAGDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/MAGDONUT/tools.c"
//#define PR_POSTINIT "PROBLEMS/MAGDONUT/postinit.c"

#endif

#if(PROBLEM==63)

#define PR_DEFINE "PROBLEMS/ORSZAG/define.h"
#define PR_BC "PROBLEMS/ORSZAG/bc.c"
#define PR_INIT "PROBLEMS/ORSZAG/init.c"
#define PR_KAPPA "PROBLEMS/ORSZAG/kappa.c"
#define PR_KAPPAES "PROBLEMS/ORSZAG/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ORSZAG/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ORSZAG/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ORSZAG/dump.c"
#define PR_TOOLS "PROBLEMS/ORSZAG/tools.c"

#endif

#if(PROBLEM==62)

#define PR_DEFINE "PROBLEMS/MAGTUBES/define.h"
#define PR_BC "PROBLEMS/MAGTUBES/bc.c"
#define PR_INIT "PROBLEMS/MAGTUBES/init.c"
#define PR_KAPPA "PROBLEMS/MAGTUBES/kappa.c"
#define PR_KAPPAES "PROBLEMS/MAGTUBES/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MAGTUBES/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MAGTUBES/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MAGTUBES/dump.c"
#define PR_TOOLS "PROBLEMS/MAGTUBES/tools.c"

#endif

#if(PROBLEM==61)

#define PR_DEFINE "PROBLEMS/MAGTESTS/define.h"
#define PR_BC "PROBLEMS/MAGTESTS/bc.c"
#define PR_INIT "PROBLEMS/MAGTESTS/init.c"
#define PR_KAPPA "PROBLEMS/MAGTESTS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MAGTESTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MAGTESTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MAGTESTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MAGTESTS/dump.c"
#define PR_TOOLS "PROBLEMS/MAGTESTS/tools.c"

#endif

#if(PROBLEM==60)

#define PR_DEFINE "PROBLEMS/JONTESTS/define.h"
#define PR_BC "PROBLEMS/JONTESTS/bc.c"
#define PR_INIT "PROBLEMS/JONTESTS/init.c"
#define PR_KAPPA "PROBLEMS/JONTESTS/kappa.c"
#define PR_KAPPAES "PROBLEMS/JONTESTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/JONTESTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/JONTESTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/JONTESTS/dump.c"
#define PR_TOOLS "PROBLEMS/JONTESTS/tools.c"

#endif

#if(PROBLEM==59)

#define PR_DEFINE "PROBLEMS/BBBLOB/define.h"
#define PR_BC "PROBLEMS/BBBLOB/bc.c"
#define PR_INIT "PROBLEMS/BBBLOB/init.c"
#define PR_KAPPA "PROBLEMS/BBBLOB/kappa.c"
#define PR_KAPPAES "PROBLEMS/BBBLOB/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BBBLOB/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BBBLOB/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BBBLOB/dump.c"
#define PR_TOOLS "PROBLEMS/BBBLOB/tools.c"

#endif

#if(PROBLEM==58)

#define PR_DEFINE "PROBLEMS/LUKE/define.h"
#define PR_BC "PROBLEMS/LUKE/bc.c"
#define PR_INIT "PROBLEMS/LUKE/init.c"
#define PR_KAPPA "PROBLEMS/LUKE/kappa.c"
#define PR_KAPPAES "PROBLEMS/LUKE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LUKE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LUKE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LUKE/dump.c"
#define PR_TOOLS "PROBLEMS/LUKE/tools.c"

#endif

#if(PROBLEM==57)

#define PR_DEFINE "PROBLEMS/LIMOTORUS/define.h"
#define PR_BC "PROBLEMS/LIMOTORUS/bc.c"
#define PR_INIT "PROBLEMS/LIMOTORUS/init.c"
#define PR_KAPPA "PROBLEMS/LIMOTORUS/kappa.c"
#define PR_KAPPAES "PROBLEMS/LIMOTORUS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/LIMOTORUS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/LIMOTORUS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/LIMOTORUS/dump.c"
#define PR_TOOLS "PROBLEMS/LIMOTORUS/tools.c"

#endif

#if(PROBLEM==56)

#define PR_DEFINE "PROBLEMS/EDDSPH/define.h"
#define PR_BC "PROBLEMS/EDDSPH/bc.c"
#define PR_INIT "PROBLEMS/EDDSPH/init.c"
#define PR_KAPPA "PROBLEMS/EDDSPH/kappa.c"
#define PR_KAPPAES "PROBLEMS/EDDSPH/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/EDDSPH/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/EDDSPH/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/EDDSPH/dump.c"
#define PR_TOOLS "PROBLEMS/EDDSPH/tools.c"

#endif

#if(PROBLEM==55)

#define PR_PREPINIT "PROBLEMS/G2STATIC/prepinit.c"
#define PR_DEFINE "PROBLEMS/G2STATIC/define.h"
#define PR_BC "PROBLEMS/G2STATIC/bc.c"
#define PR_INIT "PROBLEMS/G2STATIC/init.c"
#define PR_KAPPA "PROBLEMS/G2STATIC/kappa.c"
#define PR_KAPPAES "PROBLEMS/G2STATIC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/G2STATIC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/G2STATIC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/G2STATIC/dump.c"
#define PR_TOOLS "PROBLEMS/G2STATIC/tools.c"
#define PR_DEFS "PROBLEMS/G2STATIC/defs.h"
#define PR_FINGER "PROBLEMS/G2STATIC/finger.c"
#endif

#if(PROBLEM==54)

#define PR_DEFINE "PROBLEMS/RVDISK/define.h"
#define PR_BC "PROBLEMS/RVDISK/bc.c"
#define PR_INIT "PROBLEMS/RVDISK/init.c"
#define PR_KAPPA "PROBLEMS/RVDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDISK/dump.c"
#define PR_TOOLS "PROBLEMS/RVDISK/tools.c"

#endif

#if(PROBLEM==53)

#define PR_DEFINE "PROBLEMS/1DDONUT/define.h"
#define PR_BC "PROBLEMS/1DDONUT/bc.c"
#define PR_INIT "PROBLEMS/1DDONUT/init.c"
#define PR_KAPPA "PROBLEMS/1DDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/1DDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/1DDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/1DDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/1DDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/1DDONUT/tools.c"

#endif

#if(PROBLEM==52)

#define PR_DEFINE "PROBLEMS/SPHFLAT/define.h"
#define PR_BC "PROBLEMS/SPHFLAT/bc.c"
#define PR_INIT "PROBLEMS/SPHFLAT/init.c"
#define PR_KAPPA "PROBLEMS/SPHFLAT/kappa.c"
#define PR_KAPPAES "PROBLEMS/SPHFLAT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/SPHFLAT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/SPHFLAT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/SPHFLAT/dump.c"
#define PR_TOOLS "PROBLEMS/SPHFLAT/tools.c"

#endif

#if(PROBLEM==51)

#define PR_DEFINE "PROBLEMS/HUBBLE/define.h"
#define PR_BC "PROBLEMS/HUBBLE/bc.c"
#define PR_INIT "PROBLEMS/HUBBLE/init.c"
#define PR_KAPPA "PROBLEMS/HUBBLE/kappa.c"
#define PR_KAPPAES "PROBLEMS/HUBBLE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HUBBLE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HUBBLE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HUBBLE/dump.c"
#define PR_TOOLS "PROBLEMS/HUBBLE/tools.c"

#endif

#if(PROBLEM==50)

#define PR_DEFINE "PROBLEMS/HDWAVE/define.h"
#define PR_BC "PROBLEMS/HDWAVE/bc.c"
#define PR_INIT "PROBLEMS/HDWAVE/init.c"
#define PR_KAPPA "PROBLEMS/HDWAVE/kappa.c"
#define PR_KAPPAES "PROBLEMS/HDWAVE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HDWAVE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HDWAVE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HDWAVE/dump.c"
#define PR_TOOLS "PROBLEMS/HDWAVE/tools.c"

#endif


#if(PROBLEM==49)

#define PR_DEFINE "PROBLEMS/RVRING/define.h"
#define PR_BC "PROBLEMS/RVRING/bc.c"
#define PR_INIT "PROBLEMS/RVRING/init.c"
#define PR_KAPPA "PROBLEMS/RVRING/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVRING/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVRING/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVRING/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVRING/dump.c"
#define PR_TOOLS "PROBLEMS/RVRING/tools.c"

#endif

#if(PROBLEM==48)

#define PR_DEFINE "PROBLEMS/RVDTEST/define.h"
#define PR_BC "PROBLEMS/RVDTEST/bc.c"
#define PR_INIT "PROBLEMS/RVDTEST/init.c"
#define PR_KAPPA "PROBLEMS/RVDTEST/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDTEST/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDTEST/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDTEST/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDTEST/dump.c"
#define PR_TOOLS "PROBLEMS/RVDTEST/tools.c"

#endif

#if(PROBLEM==47)

#define PR_DEFINE "PROBLEMS/MFDONUT/define.h"
#define PR_BC "PROBLEMS/MFDONUT/bc.c"
#define PR_INIT "PROBLEMS/MFDONUT/init.c"
#define PR_KAPPA "PROBLEMS/MFDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/MFDONUT/tools.c"

#endif

#if(PROBLEM==46)

#define PR_DEFINE "PROBLEMS/MFRADNTSPH/define.h"
#define PR_BC "PROBLEMS/MFRADNTSPH/bc.c"
#define PR_INIT "PROBLEMS/MFRADNTSPH/init.c"
#define PR_KAPPA "PROBLEMS/MFRADNTSPH/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFRADNTSPH/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFRADNTSPH/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFRADNTSPH/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFRADNTSPH/dump.c"
#define PR_TOOLS "PROBLEMS/MFRADNTSPH/tools.c"

#endif

#if(PROBLEM==45)

#define PR_DEFINE "PROBLEMS/MFRADNTCYL/define.h"
#define PR_BC "PROBLEMS/MFRADNTCYL/bc.c"
#define PR_INIT "PROBLEMS/MFRADNTCYL/init.c"
#define PR_KAPPA "PROBLEMS/MFRADNTCYL/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFRADNTCYL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFRADNTCYL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFRADNTCYL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFRADNTCYL/dump.c"
#define PR_TOOLS "PROBLEMS/MFRADNTCYL/tools.c"

#endif

#if(PROBLEM==44)

#define PR_DEFINE "PROBLEMS/RADNTCYL/define.h"
#define PR_BC "PROBLEMS/RADNTCYL/bc.c"
#define PR_INIT "PROBLEMS/RADNTCYL/init.c"
#define PR_KAPPA "PROBLEMS/RADNTCYL/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADNTCYL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADNTCYL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADNTCYL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADNTCYL/dump.c"
#define PR_TOOLS "PROBLEMS/RADNTCYL/tools.c"

#endif

#if(PROBLEM==43)

#define PR_DEFINE "PROBLEMS/RVDONUTIN/define.h"
#define PR_BC "PROBLEMS/RVDONUTIN/bc.c"
#define PR_INIT "PROBLEMS/RVDONUTIN/init.c"
#define PR_KAPPA "PROBLEMS/RVDONUTIN/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDONUTIN/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDONUTIN/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDONUTIN/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDONUTIN/dump.c"
#define PR_TOOLS "PROBLEMS/RVDONUTIN/tools.c"

#endif

#if(PROBLEM==42)

#define PR_DEFINE "PROBLEMS/RVDONUT/define.h"
#define PR_BC "PROBLEMS/RVDONUT/bc.c"
#define PR_INIT "PROBLEMS/RVDONUT/init.c"
#define PR_KAPPA "PROBLEMS/RVDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RVDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RVDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RVDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RVDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/RVDONUT/tools.c"

#endif

#if(PROBLEM==41)

#define PR_DEFINE "PROBLEMS/FLATDOT/define.h"
#define PR_BC "PROBLEMS/FLATDOT/bc.c"
#define PR_INIT "PROBLEMS/FLATDOT/init.c"
#define PR_KAPPA "PROBLEMS/FLATDOT/kappa.c"
#define PR_KAPPAES "PROBLEMS/FLATDOT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FLATDOT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FLATDOT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FLATDOT/dump.c"
#define PR_TOOLS "PROBLEMS/FLATDOT/tools.c"
#define PR_FINGER "PROBLEMS/FLATDOT/finger.c"

#endif

#if(PROBLEM==40)

#define PR_DEFINE "PROBLEMS/CYLBEAMCART/define.h"
#define PR_BC "PROBLEMS/CYLBEAMCART/bc.c"
#define PR_INIT "PROBLEMS/CYLBEAMCART/init.c"
#define PR_KAPPA "PROBLEMS/CYLBEAMCART/kappa.c"
#define PR_KAPPAES "PROBLEMS/CYLBEAMCART/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/CYLBEAMCART/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/CYLBEAMCART/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/CYLBEAMCART/dump.c"
#define PR_TOOLS "PROBLEMS/CYLBEAMCART/tools.c"

#endif

#if(PROBLEM==37)

#define PR_DEFINE "PROBLEMS/MFCYLBEAM/define.h"
#define PR_BC "PROBLEMS/MFCYLBEAM/bc.c"
#define PR_INIT "PROBLEMS/MFCYLBEAM/init.c"
#define PR_KAPPA "PROBLEMS/MFCYLBEAM/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFCYLBEAM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFCYLBEAM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFCYLBEAM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFCYLBEAM/dump.c"
#define PR_TOOLS "PROBLEMS/MFCYLBEAM/tools.c"

#endif

#if(PROBLEM==36)

#define PR_DEFINE "PROBLEMS/MFDOTS/define.h"
#define PR_BC "PROBLEMS/MFDOTS/bc.c"
#define PR_INIT "PROBLEMS/MFDOTS/init.c"
#define PR_KAPPA "PROBLEMS/MFDOTS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFDOTS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFDOTS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFDOTS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFDOTS/dump.c"
#define PR_TOOLS "PROBLEMS/MFDOTS/tools.c"
#define PR_FINGER "PROBLEMS/MFDOTS/finger.c"

#endif

#if(PROBLEM==35)

#define PR_DEFINE "PROBLEMS/MFBEAMS/define.h"
#define PR_BC "PROBLEMS/MFBEAMS/bc.c"
#define PR_INIT "PROBLEMS/MFBEAMS/init.c"
#define PR_KAPPA "PROBLEMS/MFBEAMS/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFBEAMS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFBEAMS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFBEAMS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFBEAMS/dump.c"
#define PR_TOOLS "PROBLEMS/MFBEAMS/tools.c"

#endif

#if(PROBLEM==34)

#define PR_DEFINE "PROBLEMS/MFPULSE/define.h"
#define PR_BC "PROBLEMS/MFPULSE/bc.c"
#define PR_INIT "PROBLEMS/MFPULSE/init.c"
#define PR_KAPPA "PROBLEMS/MFPULSE/kappa.c"
#define PR_KAPPAES "PROBLEMS/MFPULSE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/MFPULSE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/MFPULSE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/MFPULSE/dump.c"
#define PR_TOOLS "PROBLEMS/MFPULSE/tools.c"

#endif

#if(PROBLEM==33)

#define PR_DEFINE "PROBLEMS/RADDOT/define.h"
#define PR_BC "PROBLEMS/RADDOT/bc.c"
#define PR_INIT "PROBLEMS/RADDOT/init.c"
#define PR_KAPPA "PROBLEMS/RADDOT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADDOT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADDOT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADDOT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADDOT/dump.c"
#define PR_TOOLS "PROBLEMS/RADDOT/tools.c"
#define PR_FINGER "PROBLEMS/RADDOT/finger.c"

#endif

#if(PROBLEM==32)

#define PR_DEFINE "PROBLEMS/CYLBEAM/define.h"
#define PR_BC "PROBLEMS/CYLBEAM/bc.c"
#define PR_INIT "PROBLEMS/CYLBEAM/init.c"
#define PR_KAPPA "PROBLEMS/CYLBEAM/kappa.c"
#define PR_KAPPAES "PROBLEMS/CYLBEAM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/CYLBEAM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/CYLBEAM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/CYLBEAM/dump.c"
#define PR_TOOLS "PROBLEMS/CYLBEAM/tools.c"

#endif

#if(PROBLEM==31)

#define PR_DEFINE "PROBLEMS/FLATDISK/define.h"
#define PR_BC "PROBLEMS/FLATDISK/bc.c"
#define PR_INIT "PROBLEMS/FLATDISK/init.c"
#define PR_KAPPA "PROBLEMS/FLATDISK/kappa.c"
#define PR_KAPPAES "PROBLEMS/FLATDISK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/FLATDISK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/FLATDISK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/FLATDISK/dump.c"
#define PR_TOOLS "PROBLEMS/FLATDISK/tools.c"

#endif

#if(PROBLEM==30)

#define PR_DEFINE "PROBLEMS/RADNT/define.h"
#define PR_BC "PROBLEMS/RADNT/bc.c"
#define PR_INIT "PROBLEMS/RADNT/init.c"
#define PR_KAPPA "PROBLEMS/RADNT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADNT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADNT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADNT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADNT/dump.c"
#define PR_TOOLS "PROBLEMS/RADNT/tools.c"

#endif

#if(PROBLEM==1)

#define PR_DEFINE "PROBLEMS/RADBEAM2D/define.h"
#define PR_BC "PROBLEMS/RADBEAM2D/bc.c"
#define PR_INIT "PROBLEMS/RADBEAM2D/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAM2D/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAM2D/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAM2D/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAM2D/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAM2D/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAM2D/tools.c"

#endif

#if(PROBLEM==2)

#define PR_DEFINE "PROBLEMS/RADINFALL/define.h"
#define PR_BC "PROBLEMS/RADINFALL/bc.c"
#define PR_INIT "PROBLEMS/RADINFALL/init.c"
#define PR_KAPPA "PROBLEMS/RADINFALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADINFALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADINFALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADINFALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADINFALL/dump.c"
#define PR_TOOLS "PROBLEMS/RADINFALL/tools.c"

#endif

#if(PROBLEM==3)

#define PR_DEFINE "PROBLEMS/DONUT/define.h"
#define PR_BC "PROBLEMS/DONUT/bc.c"
#define PR_INIT "PROBLEMS/DONUT/init.c"
#define PR_KAPPA "PROBLEMS/DONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUT/dump.c"
#define PR_TOOLS "PROBLEMS/DONUT/tools.c"

#endif

#if(PROBLEM==4)

#define PR_DEFINE "PROBLEMS/GEODESICINFALL/define.h"
#define PR_BC "PROBLEMS/GEODESICINFALL/bc.c"
#define PR_INIT "PROBLEMS/GEODESICINFALL/init.c"
#define PR_KAPPA "PROBLEMS/GEODESICINFALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/GEODESICINFALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/GEODESICINFALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/GEODESICINFALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/GEODESICINFALL/dump.c"
#define PR_TOOLS "PROBLEMS/GEODESICINFALL/tools.c"

#endif

#if(PROBLEM==5)

#define PR_DEFINE "PROBLEMS/EDDINFALL/define.h"
#define PR_BC "PROBLEMS/EDDINFALL/bc.c"
#define PR_INIT "PROBLEMS/EDDINFALL/init.c"
#define PR_KAPPA "PROBLEMS/EDDINFALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/EDDINFALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/EDDINFALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/EDDINFALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/EDDINFALL/dump.c"
#define PR_TOOLS "PROBLEMS/EDDINFALL/tools.c"

#endif

#if(PROBLEM==6)

#define PR_DEFINE "PROBLEMS/RADTUBE/define.h"
#define PR_BC "PROBLEMS/RADTUBE/bc.c"
#define PR_INIT "PROBLEMS/RADTUBE/init.c"
#define PR_KAPPA "PROBLEMS/RADTUBE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADTUBE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADTUBE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADTUBE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADTUBE/dump.c"
#define PR_TOOLS "PROBLEMS/RADTUBE/tools.c"

#endif

#if(PROBLEM==7)

#define PR_DEFINE "PROBLEMS/BONDI/define.h"
#define PR_BC "PROBLEMS/BONDI/bc.c"
#define PR_INIT "PROBLEMS/BONDI/init.c"
#define PR_KAPPA "PROBLEMS/BONDI/kappa.c"
#define PR_KAPPAES "PROBLEMS/BONDI/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BONDI/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BONDI/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BONDI/dump.c"
#define PR_TOOLS "PROBLEMS/BONDI/tools.c"

#endif

#if(PROBLEM==8)

#define PR_DEFINE "PROBLEMS/HDTUBE/define.h"
#define PR_BC "PROBLEMS/HDTUBE/bc.c"
#define PR_INIT "PROBLEMS/HDTUBE/init.c"
#define PR_KAPPA "PROBLEMS/HDTUBE/kappa.c"
#define PR_KAPPAES "PROBLEMS/HDTUBE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HDTUBE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HDTUBE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HDTUBE/dump.c"
#define PR_TOOLS "PROBLEMS/HDTUBE/tools.c"

#endif

#if(PROBLEM==9)

#define PR_DEFINE "PROBLEMS/HDTUBE2D/define.h"
#define PR_BC "PROBLEMS/HDTUBE2D/bc.c"
#define PR_INIT "PROBLEMS/HDTUBE2D/init.c"
#define PR_KAPPA "PROBLEMS/HDTUBE2D/kappa.c"
#define PR_KAPPAES "PROBLEMS/HDTUBE2D/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/HDTUBE2D/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/HDTUBE2D/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/HDTUBE2D/dump.c"
#define PR_TOOLS "PROBLEMS/HDTUBE2D/tools.c"

#endif

#if(PROBLEM==10)

#define PR_DEFINE "PROBLEMS/RADPULSE/define.h"
#define PR_BC "PROBLEMS/RADPULSE/bc.c"
#define PR_INIT "PROBLEMS/RADPULSE/init.c"
#define PR_KAPPA "PROBLEMS/RADPULSE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADPULSE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADPULSE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADPULSE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADPULSE/dump.c"
#define PR_TOOLS "PROBLEMS/RADPULSE/tools.c"

#endif

#if(PROBLEM==11)

#define PR_DEFINE "PROBLEMS/RADSHADOW/define.h"
#define PR_BC "PROBLEMS/RADSHADOW/bc.c"
#define PR_INIT "PROBLEMS/RADSHADOW/init.c"
#define PR_KAPPA "PROBLEMS/RADSHADOW/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADSHADOW/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADSHADOW/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADSHADOW/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADSHADOW/dump.c"
#define PR_TOOLS "PROBLEMS/RADSHADOW/tools.c"

#endif

#if(PROBLEM==12)

#define PR_DEFINE "PROBLEMS/RADATM/define.h"
#define PR_BC "PROBLEMS/RADATM/bc.c"
#define PR_INIT "PROBLEMS/RADATM/init.c"
#define PR_KAPPA "PROBLEMS/RADATM/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADATM/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADATM/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADATM/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADATM/dump.c"
#define PR_TOOLS "PROBLEMS/RADATM/tools.c"

#endif

#if(PROBLEM==13)

#define PR_DEFINE "PROBLEMS/DONUTOSC/define.h"
#define PR_BC "PROBLEMS/DONUTOSC/bc.c"
#define PR_INIT "PROBLEMS/DONUTOSC/init.c"
#define PR_KAPPA "PROBLEMS/DONUTOSC/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUTOSC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUTOSC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUTOSC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUTOSC/dump.c"
#define PR_TOOLS "PROBLEMS/DONUTOSC/tools.c"

#endif

#if(PROBLEM==14)

#define PR_DEFINE "PROBLEMS/RADWAVEBC/define.h"
#define PR_BC "PROBLEMS/RADWAVEBC/bc.c"
#define PR_INIT "PROBLEMS/RADWAVEBC/init.c"
#define PR_KAPPA "PROBLEMS/RADWAVEBC/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADWAVEBC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADWAVEBC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADWAVEBC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADWAVEBC/dump.c"
#define PR_TOOLS "PROBLEMS/RADWAVEBC/tools.c"

#endif

#if(PROBLEM==15)

#define PR_DEFINE "PROBLEMS/RADWAVE/define.h"
#define PR_BC "PROBLEMS/RADWAVE/bc.c"
#define PR_INIT "PROBLEMS/RADWAVE/init.c"
#define PR_KAPPA "PROBLEMS/RADWAVE/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADWAVE/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADWAVE/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADWAVE/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADWAVE/dump.c"
#define PR_TOOLS "PROBLEMS/RADWAVE/tools.c"

#endif

#if(PROBLEM==16)

#define PR_DEFINE "PROBLEMS/RADPULSE3D/define.h"
#define PR_BC "PROBLEMS/RADPULSE3D/bc.c"
#define PR_INIT "PROBLEMS/RADPULSE3D/init.c"
#define PR_KAPPA "PROBLEMS/RADPULSE3D/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADPULSE3D/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADPULSE3D/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADPULSE3D/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADPULSE3D/dump.c"
#define PR_TOOLS "PROBLEMS/RADPULSE3D/tools.c"

#endif

#if(PROBLEM==17)

#define PR_DEFINE "PROBLEMS/RADDBLSHADOW/define.h"
#define PR_BC "PROBLEMS/RADDBLSHADOW/bc.c"
#define PR_INIT "PROBLEMS/RADDBLSHADOW/init.c"
#define PR_KAPPA "PROBLEMS/RADDBLSHADOW/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADDBLSHADOW/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADDBLSHADOW/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADDBLSHADOW/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADDBLSHADOW/dump.c"
#define PR_TOOLS "PROBLEMS/RADDBLSHADOW/tools.c"

#endif

#if(PROBLEM==18)

#define PR_DEFINE "PROBLEMS/ATMSTATIC/define.h"
#define PR_BC "PROBLEMS/ATMSTATIC/bc.c"
#define PR_INIT "PROBLEMS/ATMSTATIC/init.c"
#define PR_KAPPA "PROBLEMS/ATMSTATIC/kappa.c"
#define PR_KAPPAES "PROBLEMS/ATMSTATIC/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ATMSTATIC/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ATMSTATIC/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ATMSTATIC/dump.c"
#define PR_TOOLS "PROBLEMS/ATMSTATIC/tools.c"

#endif

#if(PROBLEM==19)

#define PR_DEFINE "PROBLEMS/RADBEAM2DKS/define.h"
#define PR_BC "PROBLEMS/RADBEAM2DKS/bc.c"
#define PR_INIT "PROBLEMS/RADBEAM2DKS/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAM2DKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAM2DKS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAM2DKS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAM2DKS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAM2DKS/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAM2DKS/tools.c"

#endif

#if(PROBLEM==20)

#define PR_DEFINE "PROBLEMS/ATMKS/define.h"
#define PR_BC "PROBLEMS/ATMKS/bc.c"
#define PR_INIT "PROBLEMS/ATMKS/init.c"
#define PR_KAPPA "PROBLEMS/ATMKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/ATMKS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ATMKS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ATMKS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ATMKS/dump.c"
#define PR_TOOLS "PROBLEMS/ATMKS/tools.c"

#endif

#if(PROBLEM==21)

#define PR_DEFINE "PROBLEMS/DONUTKS/define.h"
#define PR_BC "PROBLEMS/DONUTKS/bc.c"
#define PR_INIT "PROBLEMS/DONUTKS/init.c"
#define PR_KAPPA "PROBLEMS/DONUTKS/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUTKS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUTKS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUTKS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUTKS/dump.c"
#define PR_TOOLS "PROBLEMS/DONUTKS/tools.c"

#endif

#if(PROBLEM==22)

#define PR_DEFINE "PROBLEMS/DONUTMKS1/define.h"
#define PR_BC "PROBLEMS/DONUTMKS1/bc.c"
#define PR_INIT "PROBLEMS/DONUTMKS1/init.c"
#define PR_KAPPA "PROBLEMS/DONUTMKS1/kappa.c"
#define PR_KAPPAES "PROBLEMS/DONUTMKS1/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/DONUTMKS1/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/DONUTMKS1/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/DONUTMKS1/dump.c"
#define PR_TOOLS "PROBLEMS/DONUTMKS1/tools.c"

#endif



#if(PROBLEM==23)

#define PR_DEFINE "PROBLEMS/ATMMKS1/define.h"
#define PR_BC "PROBLEMS/ATMMKS1/bc.c"
#define PR_INIT "PROBLEMS/ATMMKS1/init.c"
#define PR_KAPPA "PROBLEMS/ATMMKS1/kappa.c"
#define PR_KAPPAES "PROBLEMS/ATMMKS1/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/ATMMKS1/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/ATMMKS1/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/ATMMKS1/dump.c"
#define PR_TOOLS "PROBLEMS/ATMMKS1/tools.c"

#endif

#if(PROBLEM==24)

#define PR_DEFINE "PROBLEMS/RADBEAMFLAT/define.h"
#define PR_BC "PROBLEMS/RADBEAMFLAT/bc.c"
#define PR_INIT "PROBLEMS/RADBEAMFLAT/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAMFLAT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAMFLAT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAMFLAT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAMFLAT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAMFLAT/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAMFLAT/tools.c"

#endif

#if(PROBLEM==25)

#define PR_DEFINE "PROBLEMS/RDONUT/define.h"
#define PR_BC "PROBLEMS/RDONUT/bc.c"
#define PR_INIT "PROBLEMS/RDONUT/init.c"
#define PR_KAPPA "PROBLEMS/RDONUT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RDONUT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RDONUT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RDONUT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RDONUT/dump.c"
#define PR_TOOLS "PROBLEMS/RDONUT/tools.c"

#endif

#if(PROBLEM==26)

#define PR_DEFINE "PROBLEMS/RADBEAM2DKSVERT/define.h"
#define PR_BC "PROBLEMS/RADBEAM2DKSVERT/bc.c"
#define PR_INIT "PROBLEMS/RADBEAM2DKSVERT/init.c"
#define PR_KAPPA "PROBLEMS/RADBEAM2DKSVERT/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADBEAM2DKSVERT/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADBEAM2DKSVERT/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADBEAM2DKSVERT/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADBEAM2DKSVERT/dump.c"
#define PR_TOOLS "PROBLEMS/RADBEAM2DKSVERT/tools.c"

#endif

#if(PROBLEM==27)

#define PR_DEFINE "PROBLEMS/RADFLATNESS/define.h"
#define PR_BC "PROBLEMS/RADFLATNESS/bc.c"
#define PR_INIT "PROBLEMS/RADFLATNESS/init.c"
#define PR_KAPPA "PROBLEMS/RADFLATNESS/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADFLATNESS/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADFLATNESS/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADFLATNESS/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADFLATNESS/dump.c"
#define PR_TOOLS "PROBLEMS/RADFLATNESS/tools.c"

#endif

#if(PROBLEM==28)

#define PR_DEFINE "PROBLEMS/BOWSHOCK/define.h"
#define PR_BC "PROBLEMS/BOWSHOCK/bc.c"
#define PR_INIT "PROBLEMS/BOWSHOCK/init.c"
#define PR_KAPPA "PROBLEMS/BOWSHOCK/kappa.c"
#define PR_KAPPAES "PROBLEMS/BOWSHOCK/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/BOWSHOCK/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/BOWSHOCK/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/BOWSHOCK/dump.c"
#define PR_TOOLS "PROBLEMS/BOWSHOCK/tools.c"

#endif


#if(PROBLEM==29)

#define PR_DEFINE "PROBLEMS/RADWALL/define.h"
#define PR_BC "PROBLEMS/RADWALL/bc.c"
#define PR_INIT "PROBLEMS/RADWALL/init.c"
#define PR_KAPPA "PROBLEMS/RADWALL/kappa.c"
#define PR_KAPPAES "PROBLEMS/RADWALL/kappaes.c"
#define PR_OUT2GIF_2D "PROBLEMS/RADWALL/out2gif_2d.c"
#define PR_OUT2GIF_1D "PROBLEMS/RADWALL/out2gif_1d.c"
#define PR_DUMP "PROBLEMS/RADWALL/dump.c"
#define PR_TOOLS "PROBLEMS/RADWALL/tools.c"

#endif



/*********************/
//including problem specific definitions from PROBLEMS/XXX/define.h
/*********************/

#include PR_DEFINE

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

#ifndef DOFIXUPS
#define DOFIXUPS 0
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
#define EERHORATIOMIN 1.e-20
#endif

//max Erad over rho
#ifndef EERHORATIOMAX 
#define EERHORATIOMAX 1.e20
#endif

//min Erad over uint
#ifndef EEUURATIOMIN
#define EEUURATIOMIN 1.e-20
#endif

//max Erad over uint
#ifndef EEUURATIOMAX 
#define EEUURATIOMAX 1.e20
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
#define NOUTSTOP 1e7 //max n of outputs
#endif

#ifndef NSTEPSTOP
#define NSTEPSTOP 1e7 //max n of steps
#endif

#ifndef TMAX
#define TMAX 1.e50  //max time
#endif

#ifndef VERBOSE0
#define VERBOSE0 0 //verbose level
#endif

#ifndef EXPLICIT_RAD_SOURCE
#ifndef EXPLICIT_SUBSTEP_RAD_SOURCE
#ifndef IMPLICIT_FF_RAD_SOURCE
#ifndef IMPLICIT_LAB_RAD_SOURCE
#define IMPLICIT_LAB_RAD_SOURCE
#endif
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
#if (INT_ORDER==1)
#define NG 3 //number of ghost cells
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

#ifndef MULTIRADFLUID
#undef NRF
#define NRF 1
#endif

#define NV (NVMHD+4*NRF)

#else //no RADIATION

#define NV (NVMHD)
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
#define NSCALARS 5
#endif

#ifndef NAVGVARS
#ifdef BHDISK_PROBLEMTYPE
#ifdef RADIATION
#define NAVGVARS (142)
#else
#define NAVGVARS (109)
#endif
#else
#define NAVGVARS (0)
#endif
#endif

#ifndef NRADPROFILES
#define NRADPROFILES 26
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

#ifndef MKS1R0
#define MKS1R0 0.
#endif

//whether to check if the advection operator keeps entropy increasing,
//if not invert with the independently evolved entropy
#ifndef VERIFYENTROPYAFTERADVECTION
#define VERIFYENTROPYAFTERADVECTION 0
#endif

#ifndef OUTOUTPUT
#define OUTOUTPUT 1
#endif

#ifndef SCAOUTPUT
#define SCAOUTPUT 1
#endif

#ifndef RADOUTPUT
#define RADOUTPUT 1
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
#define U2PCONV 1.e-12
#endif

#ifndef RADIMPCONV
#define RADIMPCONV 1.e-8
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
#ifndef NTX
#define NTX 1 //number of tiles in X
#endif

#ifndef NTY
#define NTY 1
#endif

#ifndef NTZ
#define NTZ 1
#endif

#ifndef TNX
#define TNX (NX*NTX)
#endif

#ifndef TNY
#define TNY (NY*NTY)
#endif

#ifndef TNZ
#define TNZ (NZ*NTZ)
#endif

#ifndef NX
#define NX (TNX/NTX)
#endif 

#ifndef NY
#define NY (TNY/NTY)
#endif
 
#ifndef NZ
#define NZ (TNZ/NTZ)
#endif 

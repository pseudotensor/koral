//KORAL - problem.h
//choice of the problem plus some definitions

//available problems:

//* denotes tested and working problems
//(if you need more, email me :)

//1* RADBEAM2D - beam of light 
//2 RADINFALL - radial inflow
//3* DONUT - 2d Polish donut
//4 GEODESICINFALL - geodesic infall with blobs or not
//5 EDDINFALL - infall with flux from inside
//6* RADTUBE- radiative shock tubes as in Farris et al 09 - assumes Edd.approximation which is currently not handled
//7* BONDI - like in Fragile's paper
//8 HDTUBE - relativistic shock tube
//9 HDTUBE2D - in 2d
//10* RADPULSE - radiative blob spreading around
//11* RADSHADOW - radiative shadow
//12* RADATM - atmosphere enlighted
//13 DONUTOSC - 2d Polish donut oscillating
//14 RADWAVEBC - 1d linear rad wave imposed on boundary
//15* RADWAVE - 1d linear rad wave with periodic BC
//16 RADPULSE3D - radiative blob spreading around
//17* RADDBLSHADOW - radiative shadow with two beams inclined
//18 ATMSTATIC - hydro atmosphere 
//19* RADBEAM2DKS - beam of light in KS coordinates
//20 ATMKS - radial atmosphere infalling in KS
//21 DONUTKS - 2d Polish donut in KS
//22* DONUTMKS1 - 2d Polish donut in MKS1
//23 ATMMKS1 - radial atmosphere infalling in MKS1
//24* RADBEAMFLAT - beam of light in Cartesian 
//25* RDONUT - 2d radiative Polish donut in KS
//26* RADBEAM2DKSVERT - 2d radiative beam in r,theta plane
//27* RADFLATNESS - flat but with non-zero four-force
//28* BOWSHOCK - bow shock hydro test
//29* RADWALL - flat with wall
//30* RADNT - emission from midplane
//31* FLATDISK - emission from flat disk
//32* CYLBEAM - beam towards the axis in cylindrical
//33* RADDOT - radiating dots
//34* MFPULSE - multi fluid pulse
//35* MFBEAMS - multi fluid colliding beams
//36* MFDOTS - multi fluid radiating dots
//37* MFCYLBEAM - beam towards the axis in cylindrical with multifluids
//40* CYLBEAMCART - similar to discrete CYLBEAM but in cartesian 
//41* FLATDOT - dot which may be bigger than a dot
//42* RVDONUT - radiative and viscous dougnut
//43* RVDONUTIN - radiative and viscous dougnut inflowing
//44* RADNTCYL - emission from midplane in cylindrical
//45* MFRADNTCYL - multi fluid emission from midplane in cylindrical

#define PROBLEM 37

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

#define SMALL 1.e-50 //small number 

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
#define NRF 4

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

#ifndef U2PPREC
#define U2PPREC 1.e-7 //precision of the numerical u2p hydro solver
#endif

#ifndef IMPLABPREC
#define IMPLABPREC 1.e-7 //precision for the numerical solver in solve_implicit_lab()
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

#ifndef NX
#define NX 100 //x-resolution
#endif

#ifndef NY
#define NY 1 //y-resolution
#endif 

#ifndef NZ
#define NZ 1 //z-resolution
#endif

#ifndef INT_ORDER
#define INT_ORDER 1 //reconstruction order
#endif

#ifndef RK2STEPPING
#ifndef RK3STEPPING
#ifndef RK4STEPPING
#define RK2STEPPING //time stepping
#endif
#endif
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

#ifdef RADIATION

#ifdef MULTIRADFLUID
#define NV (6+4*NRF)
#else
#define NV 10 //number of variables
#undef NRF
#define NRF 1
#endif

#else
#define NV 6
#endif

#define NVHD 6 //number of hydro variables

#ifndef GAMMA
#define GAMMA (5./3.) //gamma
#endif


#ifndef MASS
#define MASS 1./MSUNCM //default mass of the BH used to calibrate radiation constant, Solar mass units
#endif

#ifndef U2PRADPREC
#define U2PRADPREC 1.e-5
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
#define GAMMAMAXHD 1000.
#endif

#ifndef GAMMAMAXRAD
#define GAMMAMAXRAD 1000.
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


#define LCM (MASSCM) //unit of length in cm

#define TSEC (MASSCM/CCC) //unit of time in seconds

#define GMC2CM (MASSCM) //gravitational radius in cm


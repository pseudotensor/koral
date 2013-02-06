/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
(...)

*/ 
//to transform primitives between coordinates if necessary
ldouble ggks[4][5],GGks[4][5];
calc_g_arb(xxvec,ggks,KSCOORDS);
calc_G_arb(xxvec,GGks,KSCOORDS);
trans_phd_coco(pp, pp, OUTCOORDS,KSCOORDS, xxvec,gg,GG,ggks,GGks);
					  
ldouble mvx=pp[2];
ldouble mvy=pp[3];
ldouble mvz=pp[4];
ldouble mvrel[4]={0,vx,vy,vz};
conv_vels(mvrel,mvrel,VELPRIM,VEL3,ggks,GGks);
mvx=vrel[1];
mvy=vrel[2];
mvz=vrel[3];

v1=mvx;
v2=mvy;
v3=mvz;
 
v4=get_cflag(ENTROPYFLAG,ix,iy,iz);

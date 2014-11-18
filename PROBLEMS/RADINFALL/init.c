
/*
  int
  set_initial_profile()
  {
  int ix,iy,iz;
  for(iz=0;iz<NZ;iz++)
  {
  for(iy=0;iy<NY;iy++)
  {
  for(ix=0;ix<NX;ix++)
  {
*/
struct geometry geom;
fill_geometry_arb(ix,iy,iz,&geom,BLCOORDS);
struct geometry geomint;
fill_geometry(ix,iy,iz,&geomint);

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];

xx=geom.xx;
yy=0.;
zz=0.;
ldouble gg[4][5],eup[4][4],elo[4][4],GG[4][5];

ldouble pp[NV],T;

/************************/

//ambient
mx=my=mz=0.;
pp[0]=RHO_AMB;
pp[1]=U_AMB;
pp[2]=mx;
pp[3]=my;
pp[4]=mz;
ldouble x[4]={0,geom.xx,geom.yy,geom.zz};
	   
#ifdef ANAL_PROFILE	      
//zaczynam jednak od profilu analitycznego:   
ldouble r=geom.xx;
ldouble mD=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
ldouble mE=PAR_E/(powl(r*r*sqrtl(2./r),GAMMA)*powl(1.-2./r,(GAMMA+1.)/4.));
ldouble V=sqrtl(2./r)*(1.-2./r)           ;
ldouble W=1./sqrtl(1.-V*V*geom.gg[1][1]);
ldouble mrho=mD/W;
ldouble muint=mE/W;
	      
//corrected rho:
//muint=PAR_E/(powl(r*r*sqrtl(2./r),GAMMA));
mrho=PAR_D/(r*r*sqrtl(2./r));
	      

pp[0]=mrho; pp[1]=muint; pp[2]=-V; pp[3]=pp[4]=0.;
#endif	      

#ifdef BLOB
pp[0]*=1.+BLOBMAG*exp(-(xx-BLOBX)*(xx-BLOBX)/BLOBSIG/BLOBSIG);
#endif

      //temp
//      pp[0]=1.; pp[1]=1.e-4; pp[2]=0.;

pp[5]=calc_Sfromu(pp[0],pp[1]);	      

//converting from 3vel to VELPRIM
conv_velsinprims(pp,VEL3,VELPRIM,geom.gg,geom.GG);
trans_pall_coco(pp,pp, BLCOORDS, MYCOORDS, geom.xxvec,&geom,&geomint);

	      
if(ix>=-1) //conserved required for ix=-1 only
  p2u(pp,uu,&geomint);	 


/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);

//if(isnan(get_u(p,5,ix,iy,iz))) {printf("pr: %d %d %d S: %Le\n",ix,iy,iz,0.);getchar();}

//mark initialy succesfull u2p_hot step
set_cflag(0,ix,iy,iz,0);



//KORAL - finite.c
//routines related to finite difference and grid

#include "ko.h"

/*********************************/
/****auxiliary flux limiters******/
/*********************************/
ldouble fd_flux_limiter(ldouble r)
{
  if(isnan(r)) return 1.;

  //MC
  if(FLUXLIMITER==1)
    return my_max(0., my_min(my_min(2.*r, .5 + .5*r), 2.));

  //Osher
  if(FLUXLIMITER==2) return my_max(0., my_min(r, MINMOD_THETA));

  //Koren
  if(FLUXLIMITER==3) return my_max(0.,my_min(my_min(2.*r, (1.+2.*r)/3.),2.));

  return -1.;
}

/**********************************************/
/* MINMOD limiter as in Kurganov & Tudmor *****/
/**********************************************/
ldouble minmod_flux_limiter(ldouble a,ldouble b, ldouble c)
{
  if(isnan(a) || isnan(b) || isnan(c)) return 1.;

  ldouble x1,x2;
  if(a < 0. && b < 0. && c < 0.)
    {
      if(a>b) x1=a; else x1=b;
      if(c>x1) x1=c;
      return x1;
     }
  else if(a > 0. && b > 0. && c > 0.)
    {
      if(a<b) x1=a; else x1=b;
      if(c<x1) x1=c;
      return x1;
    }
  else
    return 0.;
}

/************************************************/
/* reconstructs primitives on faces *************/
/************************************************/

ldouble DMM(ldouble X, ldouble Y)
{
  return .5*(my_sign(X)+my_sign(Y))*my_min(fabs(X),fabs(Y));
}

ldouble DM4(ldouble W,ldouble X, ldouble Y, ldouble Z)
{
  return 0.125*(my_sign(W)+my_sign(X))*fabs((my_sign(W)+my_sign(Y))*(my_sign(W)+my_sign(Z)))*my_min(my_min(fabs(W),fabs(X)),my_min(fabs(Y),fabs(Z)));
}

int
avg2point(ldouble *um2,ldouble *um1,ldouble *u0,ldouble *up1,ldouble *up2,ldouble *ul,ldouble *ur,ldouble dxm2,ldouble dxm1,ldouble dx0,ldouble dxp1,ldouble dxp2)
{
      
  ldouble r0[NV],rm1[NV],rp1[NV];
  if(INT_ORDER==1) //linear
    {
      int i;
		     
      for(i=0;i<NV;i++)
	{
	  if(FLUXLIMITER==0)
	    {	  
	      ur[i]=u0[i]+.5*minmod_flux_limiter(MINMOD_THETA*(up1[i]-u0[i]), .5*(up1[i]-um1[i]), MINMOD_THETA*(u0[i]-um1[i]));
	      ul[i]=u0[i]-.5*minmod_flux_limiter(MINMOD_THETA*(up1[i]-u0[i]), .5*(up1[i]-um1[i]), MINMOD_THETA*(u0[i]-um1[i]));
	    }
	  else
	    {
	      r0[i]=(u0[i]-um1[i])/(up1[i]-u0[i]);
	      rm1[i]=(um1[i]-um2[i])/(u0[i]-um1[i]);
	      rp1[i]=(up1[i]-u0[i])/(up2[i]-up1[i]); 
	
	      ur[i]=u0[i]+.5*fd_flux_limiter(r0[i])*(up1[i]-u0[i]);
	      ul[i]=u0[i]-.5*fd_flux_limiter(r0[i])*(up1[i]-u0[i]);      
	    }
	  if(isnan(ur[i]) || isnan (ul[i])) printf("%d %e %e %e %e %e\n",i,um2[i],um1[i],u0[i],up1[i],up2[i]);
	  //u0 remains intact - in linear reconstruction cell averaged equals cell centered
	}
    }

  else if(INT_ORDER==2) //parabolic PPM
    {
      ldouble l,r;
      int iv;      
      ldouble Z1,Z2,DX;

      ldouble dri[NV],drim1[NV],drip1[NV];

      //right face
      Z1=(dxm1+dx0)/(2.*dx0+dxp1);
      Z2=(dxp2+dxp1)/(2.*dxp1+dx0);
      DX=dxm1+dx0+dxp1+dxp2;

      for(iv=0;iv<NV;iv++)
	{
	  dri[iv]=dx0/(dxm1+dx0+dxp1)*
	    ((2.*dxm1+dx0)/(dxp1+dx0)*(up1[iv]-u0[iv])+
	     (dx0+2.*dxp1)/(dxm1+dx0)*(u0[iv]-um1[iv]));
	  drip1[iv]=dxp1/(dx0+dxp1+dxp2)*
	    ((2.*dx0+dxp1)/(dxp2+dxp1)*(up2[iv]-up1[iv])+
	     (dxp1+2.*dxp2)/(dx0+dxp1)*(up1[iv]-u0[iv]));
	  drim1[iv]=dxm1/(dxm2+dxm1+dx0)*
	    ((2.*dxm2+dxm1)/(dx0+dxm1)*(u0[iv]-um1[iv])+
	     (dxm1+2.*dx0)/(dxm2+dxm1)*(um1[iv]-um2[iv]));
	  
	  //limit slope to be monotonic
	  if( (up1[iv]-u0[iv])*(u0[iv]-um1[iv]) > 0.)
	    dri[iv]= my_min(fabs(dri[iv]),my_min(2.*fabs(u0[iv]-um1[iv]), 2.*fabs(u0[iv]-up1[iv]))) * my_sign(dri[iv]);
	  else dri[iv]=0.;
	  if( (up2[iv]-up1[iv])*(up1[iv]-u0[iv]) > 0.)
	    drip1[iv]= my_min(fabs(drip1[iv]),my_min(2.*fabs(up1[iv]-u0[iv]), 2.*fabs(up1[iv]-up2[iv]))) * my_sign(drip1[iv]);
	  else drip1[iv]=0.;
	  if( (u0[iv]-um1[iv])*(um1[iv]-um2[iv]) > 0.)
	    drim1[iv]= my_min(fabs(drim1[iv]),my_min(2.*fabs(um1[iv]-um2[iv]), 2.*fabs(um1[iv]-u0[iv]))) * my_sign(drim1[iv]);
	  else drim1[iv]=0.;
	  
	  ur[iv]=u0[iv]+dx0/(dx0+dxp1)*(up1[iv]-u0[iv])+
	    1./DX*((2.*dxp1*dx0)/(dxp1+dx0)*(Z1-Z2)*(up1[iv]-u0[iv])-
		   dx0*Z1*drip1[iv] + dxp1*Z2*dri[iv]);
	}

      //left face
      Z1=(dxm2+dxm1)/(2.*dxm1+dx0);
      Z2=(dxp1+dx0)/(2.*dx0+dxm1);
      DX=dxm2+dxm1+dx0+dxp1;

      for(iv=0;iv<NV;iv++)
	{
	      ul[iv]=um1[iv]+dxm1/(dxm1+dx0)*(u0[iv]-um1[iv])+
		1./DX*((2.*dx0*dxm1)/(dx0+dxm1)*(Z1-Z2)*(u0[iv]-um1[iv])-
		       dxm1*Z1*dri[iv] + dx0*Z2*drim1[iv]);
	}
       
      //and make sure that the parabola remains monotonic
      
      for(iv=0;iv<NV;iv++)
	{
	  if((ur[iv]-u0[iv])*(u0[iv]-ul[iv])<=0.)
	    {
	      ul[iv]=u0[iv];
	      ur[iv]=u0[iv];
	    }
	  if((ur[iv]-ul[iv])*(ul[iv]-(3.*u0[iv]-2.*ur[iv]))<0.)
	    {
	      ul[iv]=3.*u0[iv]-2.*ur[iv];
	    }
	  if((ur[iv]-ul[iv])*((3.*u0[iv]-2.*ul[iv])-ur[iv])<0.)
	    {
	      ur[iv]=3.*u0[iv]-2.*ul[iv];
	    }
	      
	  /*
	    if((u0[iv]-ul[iv])*(ul[iv]-um1[iv])<0. || (u0[iv]-ur[iv])*(ur[iv]-up1[iv])<0.) 
	    {
	    printf("non-mon parabola: %e | %e || %e || %e | %e\n",um1[iv],ul[iv],u0[iv],ur[iv],up1[iv]);
	    getchar();
	    }
	  */

	  //pass up reconstructed value at center - only if reconstructing average -> center
	  //check consistency!
	  u0[iv]=ul[iv]+.5*(ur[iv]-ul[iv] + 6.*(u0[iv]-.5*(ul[iv]+ur[iv]))*(1.-.5));
	}
    }
  else if(INT_ORDER==4) 
    //MP5 
    //TODO: currently assumes dx=const!
    {
      int iv;
      ldouble alpha=4.;
      ldouble epsm=1.e-10;
      ldouble BB2=1.3333333;

      for(iv=0;iv<NV;iv++)
	{
	  ldouble VL,VR,VOR,VMP,DJM1,DJ,DJP1,DM4JPH,DM4JMH,VUL,VAV,VMD,VLC,VMIN,VMAX;

	  //u_j+1/2,eft -> ur
	  VOR=(3.*um2[iv]-20.*um1[iv]+90.*u0[iv]+60.*up1[iv]-5.*up2[iv])/128.;
	  // VOR=(20.*um2[iv]-13.*um1[iv]+47.*u0[iv]+27.*up1[iv]-3.*up2[iv])/60.;
	  VMP=u0[iv]+DMM(up1[iv]-u0[iv],alpha*(u0[iv]-um1[iv]));
	  if((VOR-u0[iv])*(VOR-VMP)<epsm)
	    VL=VOR;
	  else
	    {
	      DJM1=um2[iv]-2.*um1[iv]+u0[iv];
	      DJ=um1[iv]-2.*u0[iv]+up1[iv];
	      DJP1=u0[iv]-2.*up1[iv]+up2[iv];
	      DM4JPH=DM4(4.*DJ-DJP1,4.*DJP1-DJ,DJ,DJP1);
	      DM4JMH=DM4(4.*DJ-DJM1,4.*DJM1-DJ,DJ,DJM1);
	      VUL=u0[iv]+alpha*(u0[iv]-um1[iv]);
	      VAV=0.5*(u0[iv]+up1[iv]);
	      VMD=VAV-0.5*DM4JPH;
	      VLC=u0[iv]+0.5*(u0[iv]-um1[iv])+BB2*DM4JMH;
	      VMIN=my_max(my_min(u0[iv],my_min(up1[iv],VMD)),my_min(u0[iv],my_min(VUL,VLC)));
	      VMAX=my_min(my_max(u0[iv],my_max(up1[iv],VMD)),my_max(u0[iv],my_max(VUL,VLC)));
	      VL=VOR+DMM(VMIN-VOR,VMAX-VOR);	      	      
	    }
	  ur[iv]=VL;
	  
	  //u_j-1/2,Right -> ul
	  VOR=(3.*up2[iv]-20.*up1[iv]+90.*u0[iv]+60.*um1[iv]-5.*um2[iv])/128.;
	  //VOR=(20.*up2[iv]-13.*up1[iv]+47.*u0[iv]+27.*um1[iv]-3.*um2[iv])/60.;
	  VMP=u0[iv]+DMM(um1[iv]-u0[iv],alpha*(u0[iv]-up1[iv]));
	  if((VOR-u0[iv])*(VOR-VMP)<epsm)
	    VL=VOR;
	  else
	    {
	      DJM1=up2[iv]-2.*up1[iv]+u0[iv];
	      DJ=up1[iv]-2.*u0[iv]+um1[iv];
	      DJP1=u0[iv]-2.*um1[iv]+um2[iv];
	      DM4JPH=DM4(4.*DJ-DJP1,4.*DJP1-DJ,DJ,DJP1);
	      DM4JMH=DM4(4.*DJ-DJM1,4.*DJM1-DJ,DJ,DJM1);
	      VUL=u0[iv]+alpha*(u0[iv]-up1[iv]);
	      VAV=0.5*(u0[iv]+um1[iv]);
	      VMD=VAV-0.5*DM4JPH;
	      VLC=u0[iv]+0.5*(u0[iv]-up1[iv])+BB2*DM4JMH;
	      VMIN=my_max(my_min(u0[iv],my_min(um1[iv],VMD)),my_min(u0[iv],my_min(VUL,VLC)));
	      VMAX=my_min(my_max(u0[iv],my_max(um1[iv],VMD)),my_max(u0[iv],my_max(VUL,VLC)));
	      VL=VOR+DMM(VMIN-VOR,VMAX-VOR);	      	      
	    }
	  ul[iv]=VL;
	}      

    }

  return 0;
}

/**************************************************/
/* saves characteristic wavespeeds from aaa[] *****/
/* to arrays corresponding to ix,iy,iz cell *******/
/**************************************************/
int
save_wavespeeds(int ix,int iy,int iz, ldouble *aaa,ldouble* max_lws)
{
  ldouble aaaxhd,aaaxrad,aaayhd,aaayrad,aaazhd,aaazrad;

	      
  set_u_scalar(ahdxl,ix,iy,iz,aaa[0]);
  set_u_scalar(ahdxr,ix,iy,iz,aaa[1]);
  set_u_scalar(ahdyl,ix,iy,iz,aaa[2]);
  set_u_scalar(ahdyr,ix,iy,iz,aaa[3]);
  set_u_scalar(ahdzl,ix,iy,iz,aaa[4]);
  set_u_scalar(ahdzr,ix,iy,iz,aaa[5]);
  set_u_scalar(aradxl,ix,iy,iz,aaa[6]);
  set_u_scalar(aradxr,ix,iy,iz,aaa[7]);
  set_u_scalar(aradyl,ix,iy,iz,aaa[8]);
  set_u_scalar(aradyr,ix,iy,iz,aaa[9]);
  set_u_scalar(aradzl,ix,iy,iz,aaa[10]);
  set_u_scalar(aradzr,ix,iy,iz,aaa[11]);

  aaaxhd=my_max(fabs(aaa[0]),fabs(aaa[1]));
  aaayhd=my_max(fabs(aaa[2]),fabs(aaa[3]));
  aaazhd=my_max(fabs(aaa[4]),fabs(aaa[5]));

  aaaxrad=my_max(fabs(aaa[6]),fabs(aaa[7]));
  aaayrad=my_max(fabs(aaa[8]),fabs(aaa[9]));
  aaazrad=my_max(fabs(aaa[10]),fabs(aaa[11]));

  set_u_scalar(ahdx,ix,iy,iz,aaaxhd);
  set_u_scalar(ahdy,ix,iy,iz,aaayhd);
  set_u_scalar(ahdz,ix,iy,iz,aaazhd);
  set_u_scalar(aradx,ix,iy,iz,aaaxrad);
  set_u_scalar(arady,ix,iy,iz,aaayrad);
  set_u_scalar(aradz,ix,iy,iz,aaazrad);


#ifdef RADIATION
  if(my_max(aaaxhd,aaaxrad)>max_lws[0]) max_lws[0]=my_max(aaaxhd,aaaxrad);
  if(my_max(aaayhd,aaayrad)>max_lws[1]) max_lws[1]=my_max(aaayhd,aaayrad);
  if(my_max(aaazhd,aaazrad)>max_lws[2]) max_lws[2]=my_max(aaazhd,aaazrad);
#else
  if(aaaxhd>max_lws[0]) max_lws[0]=aaaxhd;
  if(aaayhd>max_lws[1]) max_lws[1]=aaayhd;
  if(aaazhd>max_lws[2]) max_lws[2]=aaazhd;
#endif

  ldouble dx=get_size_x(ix,0);
  ldouble dy=get_size_x(iy,1);
  ldouble dz=get_size_x(iz,2);

  ldouble wsx,wsy,wsz;

  #ifdef RADIATION
  wsx=my_max(aaaxhd,aaaxrad);
  wsy=my_max(aaayhd,aaayrad);
  wsz=my_max(aaazhd,aaazrad);
#else
  wsx=aaaxhd;
  wsy=aaayhd;
  wsz=aaazhd;
  #endif

  //determining the time step
  //only domain cells 

  if(if_indomain(ix,iy,iz)==1) 
    {      
      ldouble tstepden,ws_ph;
     
      ws_ph=wsx*sqrt(get_g(g,1,1,ix,iy,iz));
      //#pragma omp critical
      if(wsx>max_ws[0]) max_ws[0]=wsx;
      if(ws_ph>max_ws_ph) max_ws_ph=ws_ph;
      
      ws_ph=wsy*sqrt(get_g(g,2,2,ix,iy,iz));
      //#pragma omp critical
      if(wsy>max_ws[1]) max_ws[1]=wsy;
      if(ws_ph>max_ws_ph) max_ws_ph=ws_ph;

      ws_ph=wsz*sqrt(get_g(g,3,3,ix,iy,iz));
      //#pragma omp critical
      if(wsz>max_ws[2]) max_ws[2]=wsz;
      if(ws_ph>max_ws_ph) max_ws_ph=ws_ph;

      if(NZ>1 && NY>1)
	tstepden=max_ws[0]/dx + max_ws[1]/dy + max_ws[2]/dz;
      else if(NZ==1 && NY>1)
	tstepden=max_ws[0]/dx + max_ws[1]/dy;
      else if(NY==1 && NZ>1)
	tstepden=max_ws[0]/dx + max_ws[2]/dz;
      else
	tstepden=max_ws[0]/dx;   

      //#pragma omp critical
      if(tstepden>tstepdenmax) tstepdenmax=tstepden;  
    }
  
  return 0;
}

/***************************************************/
/* the advective plus metric operator  *************/
/* applied explicitly using Lax-Friedrcich fluxes **/
/* starts from *u, copies to *ubase, updates *u, ***/
/***************************************************/
int
op_explicit(ldouble t, ldouble dt,ldouble *ubase) 
{
  int ix,iy,iz,iv,ii;

  copy_u(1.,u,ubase);


  //global
  max_ws[0]=max_ws[1]=max_ws[2]=-1.;

  //local
  ldouble max_lws[3];
  max_lws[0]=max_lws[1]=max_lws[2]=-1.;

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  //determines treatment or not of ghost cells
  int xlim,ylim,zlim;
  if(NX>1) xlim=1; else xlim=0;  
  if(NY>1) ylim=1; else ylim=0;
  if(NZ>1) zlim=1; else zlim=0;
  
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
 
  //calculates the primitives
#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain only
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      calc_primitives(ix,iy,iz,0);
    }

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  //fixup here hd after inversions
  cell_fixup_hd();

  //correct the axis if needed
#ifdef CORRECT_POLARAXIS
  correct_polaraxis();
#endif
  
  //**********************************************************************
  //* MPI ****************************************************************
  //**********************************************************************
 
  mpi_exchangedata();
 
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  //projects primitives onto ghost cells at the boundaries of the total domain
  //or calulates conserved from exchanged primitives
  set_bc(t,0);

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
	
  int superverbose=0;

  ldouble pp[NV];
     
  //calculates and saves wavespeeds
#pragma omp parallel for private(ix,iy,iz,iv,max_lws,ii) schedule (static)
  for(ii=0;ii<Nloop_1;ii++) //domain plus some ghost cells
    {
      ix=loop_1[ii][0];
      iy=loop_1[ii][1];
      iz=loop_1[ii][2]; ldouble aaa[12];

      calc_wavespeeds_lr(ix,iy,iz,aaa);	
      save_wavespeeds(ix,iy,iz,aaa,max_lws);
    }
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************


  //interpolation and flux-calculation
#pragma omp parallel for private(iy,iz,iv,ix,ii)  schedule (static) 
  for(ii=0;ii<Nloop_1;ii++) //domain plus some ghost cells
    {
      ix=loop_1[ii][0];
      iy=loop_1[ii][1];
      iz=loop_1[ii][2]; ldouble aaa[12];

      //interpolating conserved quantities
      ldouble x0[3],x0l[3],x0r[3],xm1[3],xp1[3];
      ldouble fd_r0[NV],fd_rm1[NV],fd_rp1[NV];
      ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV];
      ldouble fd_p0[NV],fd_pp1[NV],fd_pp2[NV],fd_pm1[NV],fd_pm2[NV],fd_pm3[NV],fd_pp3[NV];
      ldouble fd_s0[NV],fd_sp1[NV],fd_sp2[NV],fd_sm1[NV],fd_sm2[NV];
      ldouble fd_uLr[NV],fd_uLl[NV],fd_uRl[NV],fd_uRr[NV];
      ldouble fd_pLr[NV],fd_pLl[NV],fd_pRl[NV],fd_pRr[NV];
      ldouble fd_pr[NV],fd_pl[NV];
      ldouble fd_sLr[NV],fd_sLl[NV],fd_sRl[NV],fd_sRr[NV];
      ldouble fd_fstarl[NV],fd_fstarr[NV],fd_dul[3*NV],fd_dur[3*NV],fd_pdiffl[NV],fd_pdiffr[NV];
      ldouble fd_der[NV];
      ldouble a0[2],am1[2],ap1[2],al,ar,amax,dx;  
      ldouble ffRl[NV],ffRr[NV],ffLl[NV],ffLr[NV];
      ldouble ffl[NV],ffr[NV];
      ldouble gloc[4][5],gl[4][5],gr[4][5],Gl[4][5],Gr[4][5];
      ldouble dx0, dxm2, dxm1, dxp1, dxp2;  
      struct geometry geom;
      int i;
	      
      pick_g(ix,iy,iz,gloc);

      //**********************************************************************
      //**********************************************************************
      //x 'sweep'

#ifdef MPI4CORNERS
      if(NX>1 && iy>=-1 && iy<NY+1 && iz>=-1 && iz<NZ+1)
#else
      if(NX>1 && iy>=0 && iy<NY && iz>=0 && iz<NZ)
#endif
	{
	  x0[0]=get_x(ix,0);

	  x0l[0]=get_xb(ix,0);
	  xm1[0]=get_x(ix-1,0);
	  x0l[1]=xm1[1]=get_x(iy,1); 
	  x0l[2]=xm1[2]=get_x(iz,2);

	  x0r[0]=get_xb(ix+1,0);
	  xp1[0]=get_x(ix+1,0);
	  x0r[1]=xp1[1]=get_x(iy,1);
	  x0r[2]=xp1[2]=get_x(iz,2);

	  dx0=get_size_x(ix,0);    
	  dxm1=get_size_x(ix-1,0);    
	  dxp1=get_size_x(ix+1,0);    
	  
	  if(INT_ORDER>1)
	    {
	      dxm2=get_size_x(ix-2,0);    
	      dxp2=get_size_x(ix+2,0);    
	    }
	  
	  for(i=0;i<NV;i++)
	    {
	      //resetting derivatives
	      fd_der[i]=0.;

	      //primitives - to be interpolated
	      fd_p0[i]=get_u(p,i,ix,iy,iz);
	      fd_pp1[i]=get_u(p,i,ix+1,iy,iz);
	      fd_pm1[i]=get_u(p,i,ix-1,iy,iz);
	      if(INT_ORDER>1)
		{
		  fd_pm2[i]=get_u(p,i,ix-2,iy,iz);
		  fd_pp2[i]=get_u(p,i,ix+2,iy,iz);
		}
	    }
	 		
	  avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2);   
	  
	  if(ix>=0) //no need to calculate at left face of first GC
	    {
	      fill_geometry_face(ix,iy,iz,0,&geom);
	      check_floors_mhd(fd_pl,VELPRIM,&geom);
	      f_flux_prime(fd_pl,0,ix,iy,iz,ffl);
	    }

	  if(ix<NX)
	    {
	      fill_geometry_face(ix+1,iy,iz,0,&geom);
	      check_floors_mhd(fd_pr,VELPRIM,&geom);
	      f_flux_prime(fd_pr,0,ix+1,iy,iz,ffr);   	  
	    }

	  //saving to memory
	  for(i=0;i<NV;i++)
	    {
	      set_ubx(pbRx,i,ix,iy,iz,fd_pl[i]);
	      set_ubx(pbLx,i,ix+1,iy,iz,fd_pr[i]);

	      set_ubx(flRx,i,ix,iy,iz,ffl[i]);
	      set_ubx(flLx,i,ix+1,iy,iz,ffr[i]);		  
	    }
	}

      //**********************************************************************
      //**********************************************************************
      //y 'sweep'
  
 
#ifdef MPI4CORNERS
      if(NY>1 && ix>=-1 && ix<NX+1 && iz>=-1 && iz<NZ+1)
#else
      if(NY>1 && ix>=0 && ix<NX && iz>=0 && iz<NZ)
#endif
	{    
	  x0l[1]=get_xb(iy,1);
	  xm1[1]=get_x(iy-1,1);
	  x0l[0]=xm1[0]=get_x(ix,0); 
	  x0l[2]=xm1[2]=get_x(iz,2);

	  x0r[1]=get_xb(iy+1,1);
	  xp1[1]=get_x(iy+1,1);
	  x0r[0]=xp1[0]=get_x(ix,0);
	  x0r[2]=xp1[2]=get_x(iz,2);

	  dx0=get_size_x(iy,1);    
	  dxm1=get_size_x(iy-1,1);    
  	  dxp1=get_size_x(iy+1,1);    
	
	  if(INT_ORDER>1)
	    {
	      dxm2=get_size_x(iy-2,1);  
	      dxp2=get_size_x(iy+2,1);
	    }    
		  
	  for(i=0;i<NV;i++)
	    {
	      fd_p0[i]=get_u(p,i,ix,iy,iz);
	      fd_pp1[i]=get_u(p,i,ix,iy+1,iz);
	      fd_pm1[i]=get_u(p,i,ix,iy-1,iz);
	      if(INT_ORDER>1)
		{
		  fd_pm2[i]=get_u(p,i,ix,iy-2,iz);
		  fd_pp2[i]=get_u(p,i,ix,iy+2,iz);
		}
	    }

	  avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2);   
	  
	  if(iy>=0)
	    {
	      fill_geometry_face(ix,iy,iz,1,&geom);
	      check_floors_mhd(fd_pl,VELPRIM,&geom);
	      f_flux_prime(fd_pl,1,ix,iy,iz,ffl);
	    }

	  if(iy<NY)
	    {
	      fill_geometry_face(ix,iy+1,iz,1,&geom);
	      check_floors_mhd(fd_pr,VELPRIM,&geom);
	      f_flux_prime(fd_pr,1,ix,iy+1,iz,ffr);   	          
	    }

	  //saving to memory
	  for(i=0;i<NV;i++)
	    {
	      set_uby(pbRy,i,ix,iy,iz,fd_pl[i]);
	      set_uby(pbLy,i,ix,iy+1,iz,fd_pr[i]);

	      set_uby(flRy,i,ix,iy,iz,ffl[i]);
	      set_uby(flLy,i,ix,iy+1,iz,ffr[i]);		  
	    }
	}

      //**********************************************************************
      //**********************************************************************
      //z 'sweep'	      
#ifdef MPI4CORNERS
      if(NZ>1 && ix>=-1 && ix<NX+1 && iy>=-1 && iy<NY+1)
#else
      if(NZ>1 && ix>=0 && ix<NX && iy>=0 && iy<NY)
#endif
	{
	  x0l[2]=get_xb(iz,2);
	  xm1[2]=get_x(iz-1,2);
	  x0l[0]=xm1[0]=get_x(ix,0); 
	  x0l[1]=xm1[1]=get_x(iy,1);
  
	  x0r[2]=get_xb(iz+1,2);
	  xp1[2]=get_x(iz+1,2);
	  x0r[0]=xp1[0]=get_x(ix,0);
	  x0r[1]=xp1[1]=get_x(iy,1);

	  dx0=get_size_x(iz,2);    
	  dxm1=get_size_x(iz-1,2);    
	  dxp1=get_size_x(iz+1,2);
    
	  if(INT_ORDER>1)
	    {
	      dxm2=get_size_x(iz-2,2);    
	      dxp2=get_size_x(iz+2,2);    
	    }
		 
	  for(i=0;i<NV;i++)
	    {
	      fd_p0[i]=get_u(p,i,ix,iy,iz);
	      fd_pp1[i]=get_u(p,i,ix,iy,iz+1);
	      fd_pm1[i]=get_u(p,i,ix,iy,iz-1);

	      if(INT_ORDER>1)
		{
		  fd_pm2[i]=get_u(p,i,ix,iy,iz-2);
		  fd_pp2[i]=get_u(p,i,ix,iy,iz+2);
		}
	    }

	  avg2point(fd_pm2,fd_pm1,fd_p0,fd_pp1,fd_pp2,fd_pl,fd_pr,dxm2,dxm1,dx0,dxp1,dxp2);   

	  if(iz>=0)
	    {
	      fill_geometry_face(ix,iy,iz,2,&geom);
	      check_floors_mhd(fd_pl,VELPRIM,&geom);
	      f_flux_prime(fd_pl,2,ix,iy,iz,ffl);
	    }

	  if(iz<NZ)
	    {
	      fill_geometry_face(ix,iy,iz+1,2,&geom);
	      check_floors_mhd(fd_pr,VELPRIM,&geom);
	      f_flux_prime(fd_pr,2,ix,iy,iz+1,ffr);   	          
	    }
	  

	  //saving to memory
	  for(i=0;i<NV;i++)
	    {
	      set_ubz(pbRz,i,ix,iy,iz,fd_pl[i]);
	      set_ubz(pbLz,i,ix,iy,iz+1,fd_pr[i]);

	      set_ubz(flRz,i,ix,iy,iz,ffl[i]);
	      set_ubz(flLz,i,ix,iy,iz+1,ffr[i]);		  
	    }
	  
	}


    }


  /*
  if(PROCID==0)
    {
      print_Nvector(&get_u(p,0,0,0,-1),NV);
    }
  */
  
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

#pragma omp parallel for private(iy,iz,ix,ii)  schedule (static) 
  for(ii=0;ii<Nloop_1;ii++) //domain plus some ghost cells
    {
      ix=loop_1[ii][0];
      iy=loop_1[ii][1];
      iz=loop_1[ii][2]; ldouble aaa[12];

      //combines right - left fluxes
      f_calc_fluxes_at_faces(ix,iy,iz);
    }

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

#ifdef MAGNFIELD
  flux_ct(); //constrained transport to preserve div.B=0
#endif

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  //calculating the derivatives
#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 
      
      ldouble fd_der[NV],t_der[NV],val,ms_der[NV],ss_der[NV],rho,uint,pp[NV],uu[NV],uuold[NV],duu[NV];
	      
      ldouble gg[4][5],GG[4][5];
      ldouble tup[4][4],tlo[4][4];
      pick_T(tmuup,ix,iy,iz,tup);
      pick_T(tmulo,ix,iy,iz,tlo);
      pick_g(ix,iy,iz,gg);
      pick_G(ix,iy,iz,GG);
      ldouble gdet=gg[3][4];
	      	      
      //updating u - fluxes
      for(iv=0;iv<NV;iv++)
	{
	  ldouble flxr,flyr,flzr,flxl,flyl,flzl;

	  ldouble dx=get_size_x(ix,0);
	  ldouble dy=get_size_x(iy,1);
	  ldouble dz=get_size_x(iz,2);

	  if(INT_ORDER==1 || INT_ORDER==2) //linear or parabolic
	    {
	      flxl=get_ub(flbx,iv,ix,iy,iz,0);
	      flxr=get_ub(flbx,iv,ix+1,iy,iz,0);
	      flyl=get_ub(flby,iv,ix,iy,iz,1);
	      flyr=get_ub(flby,iv,ix,iy+1,iz,1);
	      flzl=get_ub(flbz,iv,ix,iy,iz,2);
	      flzr=get_ub(flbz,iv,ix,iy,iz+1,2);
	    }
	  if(INT_ORDER==4) //MP5
	    {
	      if(1)
		{
		  //flux reconstruction based on 3-point stencil near the boundaries
		  flxr=13./12*get_ub(flbx,iv,ix+1,iy,iz,0) -1./24. * (get_ub(flbx,iv,ix,iy,iz,0) + get_ub(flbx,iv,ix+2,iy,iz,0));
		  flxl=13./12*get_ub(flbx,iv,ix,iy,iz,0) -1./24. * (get_ub(flbx,iv,ix-1,iy,iz,0) + get_ub(flbx,iv,ix+1,iy,iz,0));
		  flyr=13./12*get_ub(flby,iv,ix,iy+1,iz,1) -1./24. * (get_ub(flby,iv,ix,iy,iz,1) + get_ub(flby,iv,ix,iy+2,iz,1));
		  flyl=13./12*get_ub(flby,iv,ix,iy,iz,1) -1./24. * (get_ub(flby,iv,ix,iy-1,iz,1) + get_ub(flby,iv,ix,iy+1,iz,1));
		  flzr=13./12*get_ub(flbz,iv,ix,iy,iz+1,2) -1./24. * (get_ub(flbz,iv,ix,iy,iz,2) + get_ub(flbz,iv,ix,iy,iz+2,2));
		  flzl=13./12*get_ub(flbz,iv,ix,iy,iz,2) -1./24. * (get_ub(flbz,iv,ix,iy,iz-1,2) + get_ub(flbz,iv,ix,iy,iz+1,2));
		}
	      else
		{
		  //flux reconstruction based on 5-point stencil 			  
		  flxr=1067./960.*get_ub(flbx,iv,ix+1,iy,iz,0) -29./480. * (get_ub(flbx,iv,ix,iy,iz,0) + get_ub(flbx,iv,ix+2,iy,iz,0)) +3./640. *(get_ub(flbx,iv,ix-1,iy,iz,0) + get_ub(flbx,iv,ix+3,iy,iz,0));
		  flxl=1067./960.*get_ub(flbx,iv,ix,iy,iz,0) -29./480. * (get_ub(flbx,iv,ix-1,iy,iz,0) + get_ub(flbx,iv,ix+1,iy,iz,0))+3./640. *(get_ub(flbx,iv,ix-2,iy,iz,0) + get_ub(flbx,iv,ix+2,iy,iz,0));
		  flyr=1067./960.*get_ub(flby,iv,ix,iy+1,iz,1) -29./480. * (get_ub(flby,iv,ix,iy,iz,1) + get_ub(flby,iv,ix,iy+2,iz,1))+3./640. *(get_ub(flby,iv,ix,iy-1,iz,1) + get_ub(flby,iv,ix,iy+3,iz,1));
		  flyl=1067./960.*get_ub(flby,iv,ix,iy,iz,1) -29./480. * (get_ub(flby,iv,ix,iy-1,iz,1) + get_ub(flby,iv,ix,iy+1,iz,1))+3./640. *(get_ub(flby,iv,ix,iy-2,iz,1) + get_ub(flby,iv,ix,iy+2,iz,1));
		  flzr=1067./960.*get_ub(flbz,iv,ix,iy,iz+1,2) -29./480. * (get_ub(flbz,iv,ix,iy,iz,2) + get_ub(flbz,iv,ix,iy,iz+2,2))+3./640. *(get_ub(flbz,iv,ix,iy,iz-1,2) + get_ub(flbz,iv,ix,iy,iz+3,2));
		  flzl=1067./960.*get_ub(flbz,iv,ix,iy,iz,2) -29./480. * (get_ub(flbz,iv,ix,iy,iz-1,2) + get_ub(flbz,iv,ix,iy,iz+1,2))+3./640. *(get_ub(flbz,iv,ix,iy,iz-2,2) + get_ub(flbz,iv,ix,iy,iz+2,2));
		}
	    }
		  
	  //unsplit scheme
	  t_der[iv]=-(flxr-flxl)/dx - (flyr-flyl)/dy - (flzr-flzl)/dz;

	  val=get_u(u,iv,ix,iy,iz)+t_der[iv]*dt;

	  if(isnan(val) || isinf(val)) {printf("i: %4d %4d %4d %d der: %e %e %e %e %e %e %e %e %e %e %e %e\n",ix,iy,iz,iv,flxr,flxl,flyr,flyl,flzr,flzl,dx,dy,dz,
					       get_u(u,iv,ix,iy,iz),get_u(p,iv,ix,iy,iz),dt);getchar();}
	  
	  set_u(u,iv,ix,iy,iz,val);	

	}

    }     
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  
  //hd fixup after advection step
  //cell_fixup_hd();

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************


  //again over cells - source terms
#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
   for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      ldouble gg[4][5],GG[4][5];
      ldouble tup[4][4],tlo[4][4];
      pick_T(tmuup,ix,iy,iz,tup);
      pick_T(tmulo,ix,iy,iz,tlo);
      pick_g(ix,iy,iz,gg);
      pick_G(ix,iy,iz,GG);
      ldouble gdet=gg[3][4];
      ldouble pp[NV],uu[NV];
	      	      
      //**********************************************************************
      //**********************************************************************
      //**********************************************************************
      //updating u - geometrical source terms
ldouble ms[NV],val;
#ifndef COUPLEMETRICWITHRADIMPLICIT
      f_metric_source_term(ix,iy,iz,ms);
     		  
      for(iv=0;iv<NV;iv++)
	{
	  val=get_u(u,iv,ix,iy,iz)+ms[iv]*dt;
	  set_u(u,iv,ix,iy,iz,val);	
	  uu[iv]=val;
	} 
#endif


      /***************************/
      //testing if entropy increased after advection and calculating post-advection primitives
      /***************************/

      calc_primitives(ix,iy,iz,VERIFYENTROPYAFTERADVECTION);
    }

   //**********************************************************************
   //**********************************************************************
   //**********************************************************************
  
   //hd fixup after advection and source
   cell_fixup_hd();

   //**********************************************************************
   //**********************************************************************
   //**********************************************************************

#ifdef RADIATION

   /************************************************************************/
   /********* explicit *** RADIATION ***************************************/
   /************************************************************************/

#ifndef SKIPRADSOURCE
#ifdef EXPLICIT_LAB_RAD_SOURCE

#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 

      //no need for it - already updated
      //calc_primitives(ix,iy,iz,0);

      explicit_rad_source_term(ix,iy,iz,dt);
    } //source terms

  //no need of radiative fixup here after source term 
  //cell_fixup_rad();

#endif //EXPLICIT_LAB
#endif //SKIPRADSOURCE
#endif //RADIATION

   
   //**********************************************************************
   //* mimics alpha-dynamo in axisymmetric sims involvin MRI ***************
   //**********************************************************************

#ifdef MIMICDYNAMO
 
  //correlates ghost cells
  mpi_exchangedata();
  set_bc(t,0);

  //mimics dynamo
  mimic_dynamo(dt); 
  
  //must be the last one as it does not update the magn. field outside the inner domain

#endif

  return GSL_SUCCESS;
}


/***************************************************/
/* the radiative implcit source term operator ******/
/* starts from *u, copies to *ubase, updates *u, ***/
/***************************************************/
int
op_implicit(ldouble t, ldouble dt,ldouble *ubase) 
{
  int ix,iy,iz,iv,ii;
  copy_u(1., u, ubase);

  //to count the average number of iteration in the implicit solver
  for(ii=0;ii<GLOBALINTSLOT_ITERIMPLTE;ii++)
    global_int_slot[ii]=0.;

  /************************************************************************/
  /******** implicit **** RADIATION ***************************************/
  /************************************************************************/

#ifdef RADIATION
#ifndef SKIPRADSOURCE
#ifdef IMPLICIT_LAB_RAD_SOURCE
  
  //again over cells - source terms
#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_0;ii++) //domain 
    {
      ix=loop_0[ii][0];
      iy=loop_0[ii][1];
      iz=loop_0[ii][2]; 
      
      calc_primitives(ix,iy,iz,0);

      implicit_lab_rad_source_term(ix,iy,iz,dt);
    } //source terms

   //fixup here after source term 
  cell_fixup_rad();

#endif //IMPLICIT_LAB_RAD_SOURCE
#endif //SKIPRADSOURCE
#endif //RADIATION

  
  return 0;
}

//************************************************************************
//* calculates fluxes at faces using Lax-Friedrichs scheme ***************
//************************************************************************
ldouble f_calc_fluxes_at_faces(int ix,int iy,int iz)
{
  int i;
  struct geometry geom;
   
  ldouble a0[2],am1[2],ap1[2],am2[2],ap2[2],ag,al,ar,amax,cmin,cmax,csLl[2],csLr[2],csRl[2],csRr[2];
  ldouble am1l[2],am1r[2],ap1l[2],ap1r[2];
  ldouble fd_u0[NV],fd_up1[NV],fd_up2[NV],fd_um1[NV],fd_um2[NV],fd_r0[NV],fd_rm1[NV],fd_rp1[NV];
  ldouble fd_uLr[NV],fd_uLl[NV],fd_uRl[NV],fd_uRr[NV],fd_fstarl[NV],fd_fstarr[NV],fd_dul[3*NV],fd_dur[3*NV],fd_pdiffl[NV],fd_pdiffr[NV];
  ldouble gdet,gg[4][5],GG[4][5],eup[4][4],elo[4][4];


  //**********************************************************************
  //**********************************************************************
  //x 'sweep'
 
#ifdef MPI4CORNERS
    if(NX>1 && iy>=-1 && iy<NY+1 && iz>=-1 && iz<NZ+1)
#else
    if(NX>1 && iy>=0 && iy<NY && iz>=0 && iz<NZ)
#endif
      {
      //characteristic wave speeds for calculating the flux on both sides of a face
      ap1l[0]=get_u_scalar(ahdxl,ix,iy,iz);
      ap1r[0]=get_u_scalar(ahdxr,ix,iy,iz);
      ap1l[1]=get_u_scalar(aradxl,ix,iy,iz);
      ap1r[1]=get_u_scalar(aradxr,ix,iy,iz);
      ap1[0]=get_u_scalar(ahdx,ix,iy,iz);
      ap1[1]=get_u_scalar(aradx,ix,iy,iz);
      am1l[0]=get_u_scalar(ahdxl,ix-1,iy,iz);
      am1r[0]=get_u_scalar(ahdxr,ix-1,iy,iz);
      am1l[1]=get_u_scalar(aradxl,ix-1,iy,iz);
      am1r[1]=get_u_scalar(aradxr,ix-1,iy,iz);
      am1[0]=get_u_scalar(ahdx,ix-1,iy,iz);
      am1[1]=get_u_scalar(aradx,ix-1,iy,iz);

      /*
	ap2[0]=get_u_scalar(ahdx,ix+1,iy,iz);
	ap2[1]=get_u_scalar(aradx,ix+1,iy,iz);
	am2[0]=get_u_scalar(ahdx,ix-2,iy,iz);
	am2[1]=get_u_scalar(aradx,ix-2,iy,iz);
      */

      //primitives at faces
      for(i=0;i<NV;i++)
	{
	  fd_uLl[i]=get_ub(pbLx,i,ix,iy,iz,0);
	  fd_uRl[i]=get_ub(pbRx,i,ix,iy,iz,0);
	}

      //converting interpolated primitives to conserved
      fill_geometry_face(ix,iy,iz,0,&geom);

#ifdef WAVESPEEDSATFACES
      ldouble aaa[12];
      calc_wavespeeds_lr_pure(fd_uLl,&geom,aaa);
      am1l[0]=aaa[0];
      am1r[0]=aaa[1];
      am1l[1]=aaa[6];
      am1r[1]=aaa[7];
      am1[0]=my_max(fabs(aaa[0]),fabs(aaa[1]));
      am1[1]=my_max(fabs(aaa[6]),fabs(aaa[7]));
      calc_wavespeeds_lr_pure(fd_uRl,&geom,aaa);
      ap1l[0]=aaa[0];
      ap1r[0]=aaa[1];
      ap1l[1]=aaa[6];
      ap1r[1]=aaa[7];
      ap1[0]=my_max(fabs(aaa[0]),fabs(aaa[1]));
      ap1[1]=my_max(fabs(aaa[6]),fabs(aaa[7]));
#endif
   
      p2u(fd_uLl,fd_uLl,&geom);
      p2u(fd_uRl,fd_uRl,&geom);
  
      //variable loop
      for(i=0;i<NV;i++)
	{
	  //choosing the proper characteristic speed - radiation decoupled from hydro
#ifdef RADIATION
	  if(i<NVMHD)      
	    {
	      ag=my_max(ap1[0],am1[0]);
	      al=my_min(ap1l[0],am1l[0]);
	      ar=my_max(ap1r[0],am1r[0]);
	    }
	  else
	    {
	      ag=my_max(ap1[1],am1[1]); 
	      al=my_min(ap1l[1],am1l[1]);
	      ar=my_max(ap1r[1],am1r[1]);
	    }
#else
	  ag=my_max(ap1[0],am1[0]); 
	  al=my_min(ap1l[0],am1l[0]);
	  ar=my_max(ap1r[0],am1r[0]);
#endif

#ifdef FULLDISSIPATION
	  ag=max_ws[0];
#endif

	  if (FLUXMETHOD==LAXF_FLUX) //Lax-Fr
	    {
	      fd_fstarl[i] = .5*(get_ub(flRx,i,ix,iy,iz,0) + get_ub(flLx,i,ix,iy,iz,0) - ag * (fd_uRl[i] - fd_uLl[i]));
	    }
	  if (FLUXMETHOD==HLL_FLUX) //HLL
	    {
	      if(al>0.) 
		fd_fstarl[i] = get_ub(flLx,i,ix,iy,iz,0);
	      else if(ar<0.)
		fd_fstarl[i] = get_ub(flLx,i,ix,iy,iz,0);
	      else
		fd_fstarl[i] = (-al*get_ub(flRx,i,ix,iy,iz,0) + ar*get_ub(flLx,i,ix,iy,iz,0) + al*ar* (fd_uRl[i] - fd_uLl[i]))/(ar-al);
	    }

	  set_ubx(flbx,i,ix,iy,iz,fd_fstarl[i]);
	}
    }


  //**********************************************************************
  //**********************************************************************
  //y 'sweep'
#ifdef MPI4CORNERS
    if(NY>1 && ix>=-1 && ix<NX+1 && iz>=-1 && iz<NZ+1)
#else
    if(NY>1 && ix>=0 && ix<NX && iz>=0 && iz<NZ)
#endif
    {
      ap1l[0]=get_u_scalar(ahdyl,ix,iy,iz);
      ap1r[0]=get_u_scalar(ahdyr,ix,iy,iz);
      ap1l[1]=get_u_scalar(aradyl,ix,iy,iz);
      ap1r[1]=get_u_scalar(aradyr,ix,iy,iz);
      ap1[0]=get_u_scalar(ahdy,ix,iy,iz);
      ap1[1]=get_u_scalar(arady,ix,iy,iz);

      am1l[0]=get_u_scalar(ahdyl,ix,iy-1,iz);
      am1r[0]=get_u_scalar(ahdyr,ix,iy-1,iz);
      am1l[1]=get_u_scalar(aradyl,ix,iy-1,iz);
      am1r[1]=get_u_scalar(aradyr,ix,iy-1,iz);
      am1[0]=get_u_scalar(ahdy,ix,iy-1,iz);
      am1[1]=get_u_scalar(arady,ix,iy-1,iz);

      for(i=0;i<NV;i++)
	{
	  fd_uLl[i]=get_ub(pbLy,i,ix,iy,iz,1);
	  fd_uRl[i]=get_ub(pbRy,i,ix,iy,iz,1);
	}

      fill_geometry_face(ix,iy,iz,1,&geom);
 
#ifdef WAVESPEEDSATFACES
      ldouble aaa[12];
      calc_wavespeeds_lr_pure(fd_uLl,&geom,aaa);
      am1l[0]=aaa[2];
      am1r[0]=aaa[3];
      am1l[1]=aaa[8];
      am1r[1]=aaa[9];
      am1[0]=my_max(fabs(aaa[2]),fabs(aaa[3]));
      am1[1]=my_max(fabs(aaa[8]),fabs(aaa[9]));
      calc_wavespeeds_lr_pure(fd_uRl,&geom,aaa);
      ap1l[0]=aaa[2];
      ap1r[0]=aaa[3];
      ap1l[1]=aaa[8];
      ap1r[1]=aaa[9];
      ap1[0]=my_max(fabs(aaa[2]),fabs(aaa[3]));
      ap1[1]=my_max(fabs(aaa[8]),fabs(aaa[9]));
#endif
 	    
      p2u(fd_uLl,fd_uLl,&geom);
      p2u(fd_uRl,fd_uRl,&geom);


      for(i=0;i<NV;i++)
	{
#ifdef RADIATION
	  if(i<NVMHD)      
	    {
	      ag=my_max(ap1[0],am1[0]);
	      al=my_min(ap1l[0],am1l[0]);
	      ar=my_max(ap1r[0],am1r[0]);
	    }
	  else
	    {
	      ag=my_max(ap1[1],am1[1]); 
	      al=my_min(ap1l[1],am1l[1]);
	      ar=my_max(ap1r[1],am1r[1]);
	    }
#else
	  ag=my_max(ap1[0],am1[0]); 
	  al=my_min(ap1l[0],am1l[0]);
	  ar=my_max(ap1r[0],am1r[0]);
#endif


#ifdef FULLDISSIPATION
	  ag=max_ws[1];
#endif
	  
	  if (FLUXMETHOD==LAXF_FLUX) //Lax-Fr
	    fd_fstarl[i] = .5*(get_ub(flRy,i,ix,iy,iz,1) + get_ub(flLy,i,ix,iy,iz,1) - ag * (fd_uRl[i] - fd_uLl[i]));
	  if (FLUXMETHOD==HLL_FLUX) //HLL
	    {
	      if(al>0.) 
		fd_fstarl[i] = get_ub(flLy,i,ix,iy,iz,1);
	      else if(ar<0.)
		fd_fstarl[i] = get_ub(flLy,i,ix,iy,iz,1);
	      else
		fd_fstarl[i] = (-al*get_ub(flRy,i,ix,iy,iz,1) + ar*get_ub(flLy,i,ix,iy,iz,1) + al*ar* (fd_uRl[i] - fd_uLl[i]))/(ar-al);
	    }
      
	  set_uby(flby,i,ix,iy,iz,fd_fstarl[i]);
 	}
    }
  else for(i=0;i<NV;i++) set_uby(flby,i,ix,iy,iz,0.);

  //**********************************************************************
  //**********************************************************************
  //z 'sweep'
#ifdef MPI4CORNERS
    if(NZ>1 && ix>=-1 && ix<NX+1 && iy>=-1 && iy<NY+1)
#else
    if(NZ>1 && ix>=0 && ix<NX && iy>=0 && iy<NY)
#endif
    {
      ap1l[0]=get_u_scalar(ahdzl,ix,iy,iz);
      ap1r[0]=get_u_scalar(ahdzr,ix,iy,iz);
      ap1l[1]=get_u_scalar(aradzl,ix,iy,iz);
      ap1r[1]=get_u_scalar(aradzr,ix,iy,iz);
      ap1[0]=get_u_scalar(ahdz,ix,iy,iz);
      ap1[1]=get_u_scalar(aradz,ix,iy,iz);

      am1l[0]=get_u_scalar(ahdzl,ix,iy,iz-1);
      am1r[0]=get_u_scalar(ahdzr,ix,iy,iz-1);
      am1l[1]=get_u_scalar(aradzl,ix,iy,iz-1);
      am1r[1]=get_u_scalar(aradzr,ix,iy,iz-1);
      am1[0]=get_u_scalar(ahdz,ix,iy,iz-1);
      am1[1]=get_u_scalar(aradz,ix,iy,iz-1);
 
      for(i=0;i<NV;i++)
	{
	  fd_uLl[i]=get_ub(pbLz,i,ix,iy,iz,2);
	  fd_uRl[i]=get_ub(pbRz,i,ix,iy,iz,2);
	}

      fill_geometry_face(ix,iy,iz,2,&geom);

#ifdef WAVESPEEDSATFACES
      ldouble aaa[12];
      calc_wavespeeds_lr_pure(fd_uLl,&geom,aaa);
      am1l[0]=aaa[4];
      am1r[0]=aaa[5];
      am1l[1]=aaa[10];
      am1r[1]=aaa[11];
      am1[0]=my_max(fabs(aaa[4]),fabs(aaa[5]));
      am1[1]=my_max(fabs(aaa[10]),fabs(aaa[11]));
      calc_wavespeeds_lr_pure(fd_uRl,&geom,aaa);
      ap1l[0]=aaa[4];
      ap1r[0]=aaa[5];
      ap1l[1]=aaa[10];
      ap1r[1]=aaa[11];
      ap1[0]=my_max(fabs(aaa[4]),fabs(aaa[5]));
      ap1[1]=my_max(fabs(aaa[10]),fabs(aaa[11]));
#endif

      p2u(fd_uLl,fd_uLl,&geom);
      p2u(fd_uRl,fd_uRl,&geom);

      for(i=0;i<NV;i++)
	{
#ifdef RADIATION
	  if(i<NVMHD)      
	    {
	      ag=my_max(ap1[0],am1[0]);
	      al=my_min(ap1l[0],am1l[0]);
	      ar=my_max(ap1r[0],am1r[0]);
	    }
	  else
	    {
	      ag=my_max(ap1[1],am1[1]);
	      al=my_min(ap1l[1],am1l[1]);
	      ar=my_max(ap1r[1],am1r[1]);
	    }
#else
	  ag=my_max(ap1[0],am1[0]); 
	  al=my_min(ap1l[0],am1l[0]);
	  ar=my_max(ap1r[0],am1r[0]);
#endif

#ifdef FULLDISSIPATION
	  ag=max_ws[2];
#endif

	  if (FLUXMETHOD==LAXF_FLUX) //Lax-Fr
	    fd_fstarl[i] = .5*(get_ub(flRz,i,ix,iy,iz,2) + get_ub(flLz,i,ix,iy,iz,2) - ag * (fd_uRl[i] - fd_uLl[i]));
	  if (FLUXMETHOD==HLL_FLUX) //HLL
	    {
	      if(al>0.) 
		fd_fstarl[i] = get_ub(flLz,i,ix,iy,iz,2);
	      else if(ar<0.)
		fd_fstarl[i] = get_ub(flLz,i,ix,iy,iz,2);
	      else
		fd_fstarl[i] = (-al*get_ub(flRz,i,ix,iy,iz,2) + ar*get_ub(flLz,i,ix,iy,iz,2) + al*ar* (fd_uRl[i] - fd_uLl[i]))/(ar-al);
	    }
  	        
	  set_ubz(flbz,i,ix,iy,iz,fd_fstarl[i]);
 	}

    }
  else for(i=0;i<NV;i++) set_ubz(flbz,i,ix,iy,iz,0.);
	
  
  return 0;
}

//***********************************************************************
//* set up grid basing on boundaries from calc_xb ***********************
//***********************************************************************
int
set_grid(ldouble *mindx,ldouble *mindy, ldouble *mindz, ldouble *maxdtfac)
{
  int i1,i2,ix,iy,iz;
  ldouble mdx,mdy,mdz,dx,dy,dz,gloc[4][5],xx[4];
  mdx=mdy=mdz=-1;
  ldouble maxdt=-1;
      

  //x
  for(i1=-NG;i1<=NX+NG;i1++)
    {
      set_xb(i1,0,calc_xb(i1,0));  
      if(i1>-NG) set_x(i1-1,0,.5*(get_xb(i1,0)+get_xb(i1-1,0)));
     }
  //y
  for(i1=-NG;i1<=NY+NG;i1++)
    {
      set_xb(i1,1,calc_xb(i1,1));  
      if(i1>-NG) set_x(i1-1,1,.5*(get_xb(i1,1)+get_xb(i1-1,1)));
   }
  //z
  for(i1=-NG;i1<=NZ+NG;i1++)
    {
      set_xb(i1,2,calc_xb(i1,2));  
      if(i1>-NG) set_x(i1-1,2,.5*(get_xb(i1,2)+get_xb(i1-1,2)));
    }

  //consistency check
#ifdef SCHWARZSCHILD
  if(get_x(-1,0)<=r_horizon_BL(0.))
    {
      printf("ix %d > %f\n",-1,get_x(-1,0));
      my_err("-1 cell inside horizon\n");
    }
#endif
#ifdef KERR
  if(get_x(-1,0)<=r_horizon_BL(BHSPIN))
    {
      printf("ix %d > %f\n",-1,get_x(-1,0));
      my_err("-1 cell inside horizon\n");

    }
#endif

  //minimal time step
  for(ix=1;ix<=NX;ix++)
    {
      for(iy=1;iy<=NY;iy++)
	{
	  for(iz=1;iz<=NZ;iz++)
	    {
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);

	      dx=fabs(get_xb(ix,0)-get_xb(ix-1,0));
	      dy=fabs(get_xb(iy,1)-get_xb(iy-1,1));
	      dz=fabs(get_xb(iz,2)-get_xb(iz-1,2));

	      if((dx<mdx || mdx<0.)) mdx=dx;
	      if((dy<mdx || mdy<0.)) mdy=dy;
	      if((dz<mdx || mdz<0.)) mdz=dz;
	      
	      //unused
	      /*
	      if((dx/sqrt(-gloc[0][0])<maxdt || maxdt<0.)) maxdt=dx/sqrt(-gloc[0][0]);
	      if((dy/sqrt(-gloc[0][0])<maxdt || maxdt<0.)) maxdt=dy/sqrt(-gloc[0][0]);
	      if((dz/sqrt(-gloc[0][0])<maxdt || maxdt<0.)) maxdt=dz/sqrt(-gloc[0][0]);
	      */
	    }
	}
    }  

  *mindx=mdx;
  *mindy=mdy;
  *mindz=mdz;
  *maxdtfac=maxdt;
  
  //**********************************************************************
  //**********************************************************************
  //auxiliary arrays to speed up parallel for loops
  //**********************************************************************
  //**********************************************************************

  //inside domain only
  Nloop_0=0;
  loop_0=(int **)malloc(sizeof(int*));
  loop_0[0]=(int *)malloc(3*sizeof(int));

  for(ix=0;ix<NX;ix++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {	
	      loop_0[Nloop_0][0]=ix;
	      loop_0[Nloop_0][1]=iy;
	      loop_0[Nloop_0][2]=iz;

	      Nloop_0++;
	      
	      loop_0=(int **)realloc(loop_0,(Nloop_0+1)*sizeof(int*));
	      loop_0[Nloop_0]=(int *)malloc(3*sizeof(int));	      
	    }
	}
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_0,Nloop_0);
#endif

  //**********************************************************************
  //**********************************************************************
  //inside + ghost cells - number depending on the order of reconstruction
  //used to indicate where calculate fluxes
  int xlim,ylim,zlim;
  int lim;

  if(INT_ORDER==1) lim=1;
  if(INT_ORDER==2) lim=1;
  if(INT_ORDER==4) lim=2;

  if(NX>1) xlim=lim; else xlim=0;  
  if(NY>1) ylim=lim; else ylim=0;
  if(NZ>1) zlim=lim; else zlim=0;

  Nloop_1=0;
  loop_1=(int **)malloc(sizeof(int*));
  loop_1[0]=(int *)malloc(3*sizeof(int));

  for(ix=-xlim;ix<NX+xlim;ix++)
    {
      for(iy=-ylim;iy<NY+ylim;iy++)
	{
	  for(iz=-zlim;iz<NZ+zlim;iz++)
	    {	
	      //if(if_outsidegc(ix,iy,iz)==1) continue; //avoid corners

	      loop_1[Nloop_1][0]=ix;
	      loop_1[Nloop_1][1]=iy;
	      loop_1[Nloop_1][2]=iz;

	      Nloop_1++;
	      
	      loop_1=(int **)realloc(loop_1,(Nloop_1+1)*sizeof(int*));
	      loop_1[Nloop_1]=(int *)malloc(3*sizeof(int));	      
	    }
	}
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_1,Nloop_1);
#endif


  //**********************************************************************
  //**********************************************************************
  //only ghost cells (with no corners)

  //reduction of size basing on the dimension
  //for constrained transport fluxes copied onto missing dimensions

  if(NX>1) xlim=NG; else xlim=0;  
  if(NY>1) ylim=NG; else ylim=0;
  if(NZ>1) zlim=NG; else zlim=0;

  Nloop_2=0;
  loop_2=(int **)malloc(sizeof(int*));
  loop_2[0]=(int *)malloc(3*sizeof(int));

  for(ix=-xlim;ix<NX+xlim;ix++)
    {
      for(iy=-ylim;iy<NY+ylim;iy++)
	{
	  for(iz=-zlim;iz<NZ+zlim;iz++)
	    {	 
	      //within domain:
	      if(if_indomain(ix,iy,iz)==1) continue;

	      //at the corners
	      //I commented it out to precalculate metric in the corners but may affect something, who knows? :)
	      //if(if_outsidegc(ix,iy,iz)==1) continue;

	      loop_2[Nloop_2][0]=ix;
	      loop_2[Nloop_2][1]=iy;
	      loop_2[Nloop_2][2]=iz;

	      Nloop_2++;
	      
	      loop_2=(int **)realloc(loop_2,(Nloop_2+1)*sizeof(int*));
	      loop_2[Nloop_2]=(int *)malloc(3*sizeof(int));
	    }
	}
    }	
  
  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_2,Nloop_2);
#endif

  //**********************************************************************
  //**********************************************************************
  //domain and all ghost cells (no corners)
  Nloop_02=Nloop_0+Nloop_2;
  loop_02=(int **)malloc(Nloop_02*sizeof(int*));
  for(ix=0;ix<Nloop_0;ix++)
    {
      loop_02[ix]=(int *)malloc(3*sizeof(int));
      loop_02[ix][0]=loop_0[ix][0];
      loop_02[ix][1]=loop_0[ix][1];
      loop_02[ix][2]=loop_0[ix][2];
    }
  for(ix=0;ix<Nloop_2;ix++)
    {
      loop_02[ix+Nloop_0]=(int *)malloc(3*sizeof(int));
      loop_02[ix+Nloop_0][0]=loop_2[ix][0];
      loop_02[ix+Nloop_0][1]=loop_2[ix][1];
      loop_02[ix+Nloop_0][2]=loop_2[ix][2];
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_02,Nloop_02);
#endif

  //**********************************************************************
  //**********************************************************************
  //1-deep surfaces on corners only
  if(NX>1) xlim=NG; else xlim=0;  
  if(NY>1) ylim=NG; else ylim=0;
  if(NZ>1) zlim=NG; else zlim=0;

  Nloop_3=0;
  loop_3=(int **)malloc(sizeof(int*));
  loop_3[0]=(int *)malloc(3*sizeof(int));

  for(ix=-xlim;ix<NX+xlim;ix++)
    {
      for(iy=-ylim;iy<NY+ylim;iy++)
	{
	  for(iz=-zlim;iz<NZ+zlim;iz++)
	    {	  
	      //if outside the corners skip
	      if(if_outsidegc(ix,iy,iz)==0) continue;

	      //now check if in the surface layer of a corner
	      int dix,diy,diz;
	      if(ix<0) dix=-ix;
	      else dix=ix-NX+1;
	      if(iy<0) diy=-iy;
	      else diy=iy-NY+1;
	      if(iz<0) diz=-iz;
	      else diz=iz-NZ+1;	      
	      if(dix!=1 && diy!=1 && diz!=1) continue;

	      //printf("%4d %4d %4d %d %d %d\n",ix,iy,iz,dix,diy,diz);
	      
	      loop_3[Nloop_3][0]=ix;
	      loop_3[Nloop_3][1]=iy;
	      loop_3[Nloop_3][2]=iz;

	      Nloop_3++;
	      
	      loop_3=(int **)realloc(loop_3,(Nloop_3+1)*sizeof(int*));
	      loop_3[Nloop_3]=(int *)malloc(3*sizeof(int));	      
	    }
	}
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_3,Nloop_3);
#endif

  //**********************************************************************
  //**********************************************************************
  //all corners of the domain 
  
  Nloop_4=0;
  loop_4=(int **)malloc(sizeof(int*));
  loop_4[0]=(int *)malloc(3*sizeof(int));

  if(NY>1) ylim=NY; else ylim=0;
  if(NZ>1) zlim=NZ; else zlim=0;
  for(ix=0;ix<=NX;ix++)
    {
      for(iy=0;iy<=ylim;iy++)
	{
	  for(iz=0;iz<=zlim;iz++)
	    {	
	      loop_4[Nloop_4][0]=ix;
	      loop_4[Nloop_4][1]=iy;
	      loop_4[Nloop_4][2]=iz;

	      Nloop_4++;
	      
	      loop_4=(int **)realloc(loop_4,(Nloop_4+1)*sizeof(int*));
	      loop_4[Nloop_4]=(int *)malloc(3*sizeof(int));	      
	    }
	}
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_4,Nloop_4);
#endif
      
  //**********************************************************************
  //**********************************************************************
  //all corners
  
  Nloop_5=0;
  loop_5=(int **)malloc(sizeof(int*));
  loop_5[0]=(int *)malloc(3*sizeof(int));

   for(ix=-NGCX;ix<=NX+NGCX;ix++)
    {
      for(iy=-NGCY;iy<=NY+NGCY;iy++)
	{
	  for(iz=-NGCZ;iz<=NZ+NGCZ;iz++)
	    {	
	      loop_5[Nloop_5][0]=ix;
	      loop_5[Nloop_5][1]=iy;
	      loop_5[Nloop_5][2]=iz;

	      Nloop_5++;
	      
	      loop_5=(int **)realloc(loop_5,(Nloop_5+1)*sizeof(int*));
	      loop_5[Nloop_5]=(int *)malloc(3*sizeof(int));	      
	    }
	}
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_5,Nloop_5);
#endif
      
  //**********************************************************************
  //**********************************************************************
  //inner domain plus 1-cell layer including corners
  
  Nloop_6=0;
  loop_6=(int **)malloc(sizeof(int*));
  loop_6[0]=(int *)malloc(3*sizeof(int));

  if(NY>1) ylim=1 ; else ylim=0;
  if(NZ>1) zlim=1 ; else zlim=0;

  for(ix=-1;ix<NX+1;ix++)
    {
      for(iy=-ylim;iy<NY+ylim;iy++)
	{
	  for(iz=-zlim;iz<=NZ+zlim;iz++)
	    {
	      loop_6[Nloop_6][0]=ix;
	      loop_6[Nloop_6][1]=iy;
	      loop_6[Nloop_6][2]=iz;

	      Nloop_6++;
	      
	      loop_6=(int **)realloc(loop_6,(Nloop_6+1)*sizeof(int*));
	      loop_6[Nloop_6]=(int *)malloc(3*sizeof(int));	      
	    }
	}
    }

  //shuffling:
#if (SHUFFLELOOPS)
  shuffle_loop(loop_6,Nloop_6);
#endif
      

  return 0;
}

//**********************************************************************
//* prints grid **********************************************************
//**********************************************************************

//prints grid
int
print_grid(ldouble min_dx, ldouble min_dy, ldouble min_dz)
{
  int i1;
  for(i1=-NG;i1<=NX+NG-1;i1++)
    printf("x: %6d %8.3f|%8.3f|%8.3f (%8.3f)\n",i1,get_xb(i1,0),get_x(i1,0),get_xb(i1+1,0),get_size_x(i1,0));
  printf("\n");
  for(i1=-NG;i1<=NY+NG-1;i1++)
    printf("y: %6d %8.3f|%8.3f|%8.3f (%8.3f)\n",i1,get_xb(i1,1),get_x(i1,1),get_xb(i1+1,1),get_size_x(i1,1));
  printf("\n");
  for(i1=-NG;i1<=NZ+NG-1;i1++)
    printf("z: %6d %8.3f|%8.3f|%8.3f (%8.3f)\n",i1,get_xb(i1,2),get_x(i1,2),get_xb(i1+1,2),get_size_x(i1,2));

  printf("\n min_dx = %8.3f\n min_dy = %8.3f\n min_dz = %8.3f\n",min_dx,min_dy,min_dz);
      
  //getchar(); 

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************

//routines to handle arrays - many of them replaced by defs in ko.h

//**********************************************************************
//**********************************************************************
//**********************************************************************


//returns four-vector of coordinates
//so far neglects time
int 
get_xx(int ix,int iy,int iz,ldouble *xx)
{
  xx[0]=0.;
  xx[1]=get_x(ix,0);
  xx[2]=get_x(iy,1);
  xx[3]=get_x(iz,2);
  return 0;
}

//sets cell center location
int set_x(int ic, int idim,ldouble val)
{  
  if(idim==0)
    x[ic+NG]=val;
  if(idim==1)
    x[ic+NG + NX+2*NG]=val;
  if(idim==2)
    x[ic+NG + NX+2*NG + NY+2*NG]=val;
  return 0;
}
/*
//returns location of cell (i) boundaries (i,i+1)
ldouble get_xb(int ic, int idim)
{
  if(ic<-NG || (idim==0 && ic>NX+NG) || (idim==1 && ic>NY+NG) || (idim==2 && ic>NZ+NG)) my_err("blont w get_xb - index ouf of range");
  if(idim==0)
    return xb[ic+NG];
  if(idim==1)
    return xb[ic+NG + NX+2*NG + 1];
  if(idim==2)
    return xb[ic+NG + NX+2*NG + 1 + NY+2*NG + 1];

  return -1;
}
*/

//returns size of a cell
ldouble get_size_x(int ic, int idim)
{
  if(ic<-NG || (idim==0 && ic>NX-1+NG) || (idim==1 && ic>NY-1+NG) || (idim==2 && ic>NZ-1+NG)) my_err("blont w get_size_x - index ouf of range");
  return get_xb(ic+1,idim)-get_xb(ic,idim);
}

//sets cell boundaries
int set_xb(int ic, int idim,ldouble val)
{  
  if(ic<-NG || (idim==0 && ic>NX+NG) || (idim==1 && ic>NY+NG) || (idim==2 && ic>NZ+NG)) my_err("blont w set_xb - index ouf of range");
  if(idim==0)
    xb[ic+NG]=val;
  if(idim==1)
    xb[ic+NG + NX+2*NG +1]=val;
  if(idim==2)
    xb[ic+NG + NX+2*NG + 1 + NY+2*NG + 1]=val;
  return 0;
}

//deals with arrays [NX+NG x NY+NG x NZ+NG x NV] - cell centers 
/*
ldouble get_u(ldouble* uarr, int iv, int ix, int iy, int iz)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_u - index ouf of range");
  
  return uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV];
}
*/
//deals with arrays [NX+NG x NY+NG x NZ+NG x NV] - cell centers 
 /*
int set_u(ldouble* uarr,int iv,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_u - index ouf of range");
  
  uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV] = value;
  return 0;
}
 */
  /*
//deals with arrays [NX+NG x NY+NG x NZ+NG x gSIZE] - cell centers metric
ldouble get_g(ldouble* uarr, int i, int j, int ix, int iy, int iz)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_g - index ouf of range");
  
  return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE];
}
*/

//deals with arrays [NX+NG x NY+NG x NZ+NG x gSIZE] - cell centers metric
int set_g(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_g - index ouf of range");
  
  uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZ(iz)+NGCZ)*(SY)*(SX)*gSIZE] = value;
  return 0;
}

//deals with arrays [NX+NG x NY+NG x NZ+NG x 16] - 4x4 tensors
int set_T(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_T - index ouf of range");
  
  uarr[i*4+j + (ix+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY)*(SX)*16] = value;
  return 0;
}


//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x NV] - cell boundaries in idim
int set_ub(ldouble* uarr,int iv,int ix,int iy,int iz,ldouble value,int idim)
{
  if(idim==0)
    {
      uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX+1)*NV + (iZ(iz)+NGCZ)*(SY)*(SX+1)*NV] = value;
    }
  if(idim==1)
    {
      uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY+1)*(SX)*NV] = value;
    }
  if(idim==2)
    {
      uarr[iv + (iX(ix)+NGCX)*NV + (iY(iy)+NGCY)*(SX)*NV + (iZ(iz)+NGCZ)*(SY)*(SX)*NV] = value;
    }
  return 0;
}


/*
//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x NV] - cell boundaries in idim
ldouble get_ub(ldouble* uarr,int iv,int ix,int iy,int iz,int idim)
{
  if(idim==0)
    {
      if(ix<-NG || ix>NX+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_ub x - index ouf of range");  
      return uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG+1)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*NV];
    }
  if(idim==1)
    {
      if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_ub y - index ouf of range");  
      return uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*NV];
    }
  if(idim==2)
    {
      if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ+NG) my_err("blont w get_ub z - index ouf of range");  
      return uarr[iv + (ix+NG)*NV + (iy+NG)*(NX+2*NG)*NV + (iz+NG)*(NY+2*NG)*(NX+2*NG)*NV];
    }
  return 0;
}
*/
//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x gSIZE] - cell boundaries in idim metric
int set_gKrb(int i,int j,int k,int ix,int iy,int iz,ldouble val,int idim)
{
  if(idim==0)
    {
      gKrbx[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX+1)*64 + (iZ(iz)+NGCZ)*(SY)*(SX+1)*64]=val;
    }
  if(idim==1)
    {
      gKrby[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX)*64 + (iZ(iz)+NGCZ)*(SY+1)*(SX)*64]=val;
    }
  if(idim==2)
    {
      gKrbz[i*4*4+j*4+k + (iX(ix)+NGCX)*64 + (iY(iy)+NGCY)*(SX)*64 + (iZ(iz)+NGCZ)*(SY)*(SX)*64]=val;
    }
  return 0;
}


//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x gSIZE] - cell boundaries in idim metric
int set_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim)
{
  if(idim==0)
    {
      uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX+1)*gSIZE + (iZ(iz)+NGCZ)*(SY)*(SX+1)*gSIZE] = value;
    }
  if(idim==1)
    {
      uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZ(iz)+NGCZ)*(SY+1)*(SX)*gSIZE] = value;
    }
  if(idim==2)
    {
      uarr[i*5+j + (iX(ix)+NGCX)*gSIZE + (iY(iy)+NGCY)*(SX)*gSIZE + (iZ(iz)+NGCZ)*(SY)*(SX)*gSIZE] = value;
    }
  return 0;
}

//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x 16] - tensors at cell boundaries in idim metric
int set_Tb(ldouble* uarr,int i,int j,int ix,int iy,int iz,ldouble value,int idim)
{
  if(idim==0)
    {
      uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX+1)*16 + (iZ(iz)+NGCZ)*(SY)*(SX+1)*16] = value;
    }
  if(idim==1)
    {
      uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY+1)*(SX)*16] = value;
    }
  if(idim==2)
    {
      uarr[i*4+j + (iX(ix)+NGCX)*16 + (iY(iy)+NGCY)*(SX)*16 + (iZ(iz)+NGCZ)*(SY)*(SX)*16] = value;
    }
  return 0;
}

/*
//deals with arrays ~[NX+NG+1 x NY+NG x NZ+NG x gSIZE] - cell boundaries in idim metric
ldouble get_gb(ldouble* uarr,int i,int j,int ix,int iy,int iz,int idim)
{
  if(idim==0)
    {
      if(ix<-NG || ix>NX+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_gb x - index ouf of range");  
      return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG+1)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG+1)*gSIZE];
    }
  if(idim==1)
    {
      if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_gb y - index ouf of range");  
      return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG+1)*(NX+2*NG)*gSIZE];
    }
  if(idim==2)
    {
      if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ+NG) my_err("blont w get_gb z - index ouf of range");  
      return uarr[i*5+j + (ix+NG)*gSIZE + (iy+NG)*(NX+2*NG)*gSIZE + (iz+NG)*(NY+2*NG)*(NX+2*NG)*gSIZE];
    }
  return 0;
}
*/
/*
//deals with arrays [NX+NG x NY+NG x NZ+NG] - cell centers 
int set_u_scalar(ldouble* uarr,int ix,int iy,int iz,ldouble value)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w set_u_scalar - index ouf of range");

  uarr[ix+NG + (iy+NG)*(NX+2*NG) + (iz+NG)*(NY+2*NG)*(NX+2*NG)] = value;
  return 0;
}
*/

//deals with arrays [NX+NG x NY+NG x NZ+NG x NV] - cell centers 
 /*
ldouble get_u_scalar(ldouble* uarr,int ix,int iy,int iz)
{
  if(ix<-NG || ix>NX-1+NG || iy<-NG || iy>NY-1+NG || iz<-NG || iz>NZ-1+NG) my_err("blont w get_u_scalar - index ouf of range");
  
  //TODO something better, so far skipping calculating wave speeds at ghosts
  //as it is easier due to extrapolation of primitives quantities 

  //printf("%4d %4d %4d\n",ix,iy,iz); getchar();
  
  return uarr[ix+NG + (iy+NG)*(NX+2*NG) + (iz+NG)*(NY+2*NG)*(NX+2*NG)];
}
*/
//**********************************************************************
//**********************************************************************
//**********************************************************************

//array multiplication
//uu2=factor*uu1
int 
copy_u_core(ldouble factor,ldouble *uu1,ldouble* uu2, int N )
{
  int i;
#pragma omp parallel for private (i) 
  for (i=0;i<N;i++)
    uu2[i]=uu1[i]*factor;
  return 0;
}

int 
copy_u(ldouble factor,ldouble *uu1,ldouble* uu2 )
{
  copy_u_core(factor,uu1,uu2,SX*SY*SZ*NV);
  return 0;
}

//array multiplication plus addition
//uu3=f1*uu1+f2*uu2
int 
add_u_core(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3,int N)
{
  int i;
#pragma omp parallel for private (i) 
  for (i=0;i<N;i++)
    uu3[i]=uu1[i]*f1+uu2[i]*f2;
  return 0;
}

int 
add_u(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble *uu3)
{
  add_u_core(f1,uu1,f2,uu2,uu3,SX*SY*SZ*NV);
  return 0;
}

//array multiplication plus addition on 3 matrices
//uu3=f1*uu1+f2*uu2
int 
add_u_core_3(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble f3, ldouble *uu3, ldouble *uu4,int N)
{
  int i;
#pragma omp parallel for private (i) 
  for (i=0;i<N;i++)
    uu4[i]=uu1[i]*f1+uu2[i]*f2+uu3[i]*f3;
  return 0;
}

int 
add_u_3(ldouble f1, ldouble* uu1, ldouble f2, ldouble *uu2, ldouble f3, ldouble *uu3, ldouble *uu4)
{
  add_u_core_3(f1,uu1,f2,uu2,f3,uu3,uu4,SX*SY*SZ*NV);
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************

//checks if cell is inside main domain
int 
if_indomain(int ix,int iy,int iz)
{
  if(ix>=0 && ix<NX && iy>=0 && iy<NY && iz>=0 && iz<NZ) return 1;
  else return 0;
}

//checks if cell outside both domain and ghostcells, i.e. if cells in corners
int 
if_outsidegc(int ix,int iy,int iz)
{  
  if(((ix<0 || ix>=NX) && (iy>=0 && iy<NY) && (iz>=0 && iz<NZ)) || 
     ((ix>=0 && ix<NX) && (iy<0 || iy>=NY) && (iz>=0 && iz<NZ)) || 
     ((ix>=0 && ix<NX) && (iy>=0 && iy<NY) && (iz<0 || iz>=NZ)) ||
     (ix>=0 && ix<NX && iy>=0 && iy<NY && iz>=0 && iz<NZ))
    return 0;
  else
    return 1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************

//boundary conditions - sets conserved in the ghost cells
int set_bc_core(int ix,int iy,int iz,double t,ldouble *uval,ldouble *pval,int ifinit,int BCtype)
{
  int iix,iiy,iiz,iv;
  iix=ix;
  iiy=iy;
  iiz=iz;
          
#ifdef SPECIFIC_BC  //BC specific for given problem
  calc_bc(ix,iy,iz,t,uval,pval,ifinit,BCtype);

#else  
  //standard BC  
  if(BCtype==XBCLO || BCtype==XBCHI)
    {       
#ifdef PERIODIC_XBC
      iix=ix;
      if(ix<0) iix=ix+NX;
      if(ix>NX-1) iix=ix-NX;
#endif
#ifdef COPY_XBC
      iix=ix;
      if(ix<0) iix=0;
      if(ix>NX-1) iix=NX-1;
#endif
    }

  if(BCtype==YBCLO || BCtype==YBCHI)
    {       
#ifdef PERIODIC_YBC
      iiy=iy;
      if(iy<0) iiy=iy+NY;
      if(iy>NY-1) iiy=iy-NY;
      if(NY<NG) iiy=0;
#endif
#ifdef COPY_YBC
      iiy=iy;
      if(iy<0) iiy=0;
      if(iy>NY-1) iiy=NY-1;
#endif
    }

  if(BCtype==ZBCLO || BCtype==ZBCHI)
    {       
#ifdef PERIODIC_ZBC
      iiz=iz;
      if(iz<0) iiz=iz+NZ;
      if(iz>NZ-1) iiy=iz-NZ;
      if(NZ<NG) iiz=0;
#endif
#ifdef COPY_ZBC
      iiz=iz;
      if(iz<0) iiz=0;
      if(iz>NZ-1) iiz=NZ-1;
#endif
    }


  for(iv=0;iv<NV;iv++)
    pval[iv]=get_u(p,iv,iix,iiy,iiz);

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  p2u(pval,uval,&geom);
     
#endif //SPECIFIC_BC   

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************

//boundary conditions - sets conserved in the ghost cells
int set_bc(ldouble t,int ifinit)
{
  int ix,iy,iz,ii,iv;


  //first fill the GC with no corners
#pragma omp parallel for private(ix,iy,iz,iv,ii) schedule (static)
  for(ii=0;ii<Nloop_2;ii++) //ghost cells only, no corners
    {
      ix=loop_2[ii][0];
      iy=loop_2[ii][1];
      iz=loop_2[ii][2];

      //type of BC
      int BCtype;
      if(ix<0) BCtype=XBCLO;
      if(ix>=NX) BCtype=XBCHI;
      if(iy<0) BCtype=YBCLO;
      if(iy>=NY) BCtype=YBCHI;
      if(iz<0) BCtype=ZBCLO;
      if(iz>=NZ) BCtype=ZBCHI;

      if(mpi_isitBC(BCtype)==0) //this border exchanged through MPI
	{
	  struct geometry geom;
	  fill_geometry(ix,iy,iz,&geom);
	  p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	}
      else //need for real BC
	{
	  ldouble uval[NV],pval[NV];
	  set_bc_core(ix,iy,iz,t,uval,pval,ifinit,BCtype);
	  
	  for(iv=0;iv<NV;iv++)
	    {
	      set_u(u,iv,ix,iy,iz,uval[iv]);
	      set_u(p,iv,ix,iy,iz,pval[iv]);	      
	    }
	}
    }

#ifdef MPI4CORNERS
  //now calculate conserved in corners - 
  struct geometry geom;
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCLO)==0) 
    {
      for(ix=-NG;ix<0;ix++)
	for(iy=-NG;iy<0;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {	    
	      fill_geometry(ix,iy,iz,&geom);
	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
    }
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCLO)==0) 
    {
      for(ix=NX;ix<NX+NG;ix++)
	for(iy=-NG;iy<0;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {	    
	      fill_geometry(ix,iy,iz,&geom);
	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
    }
  if(mpi_isitBC(XBCLO)==0 && mpi_isitBC(YBCHI)==0) 
    {
      for(ix=-NG;ix<0;ix++)
	for(iy=NY;iy<NY+NG;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {	    
	      fill_geometry(ix,iy,iz,&geom);
	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
    }
  if(mpi_isitBC(XBCHI)==0 && mpi_isitBC(YBCHI)==0) 
    {
      for(ix=NX;ix<NX+NG;ix++)
	for(iy=NY;iy<NY+NG;iy++)
	  for(iz=0;iz<NZ;iz++)
	    {	    
	      fill_geometry(ix,iy,iz,&geom);
	      p2u(&get_u(p,0,ix,iy,iz),&get_u(u,0,ix,iy,iz),&geom);
	    }
    }
  
  //corners of the whole domain are never real BC so need to fill them with something
  int xlim,ylim,zlim;
  int lim,i,j;

  if(NZ>1)
    {
      printf("filling corners not implemented for NZ>1\n");exit(-1);
    }

  iz=0;

  //total corners, filling one cell deep surfaces
  if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCLO)==1)
    {

      for(i=0;i<NG-1;i++)
	{
	  PLOOP(iv)
	  {
	    set_u(p,iv,-NG+i,-1,iz,get_u(p,iv,-NG+1,0,iz));
	    set_u(p,iv,-1,-NG+i,iz,get_u(p,iv,0,-NG+1,iz));
	  }
	  fill_geometry(-NG+i,-1,iz,&geom);
	  p2u(&get_u(p,0,-NG+i,-1,iz),&get_u(u,0,-NG+i,-1,iz),&geom);
	  fill_geometry(-1,-NG+i,iz,&geom);
	  p2u(&get_u(p,0,-1,-NG+i,iz),&get_u(u,0,-1,-NG+i,iz),&geom);
	}
      
      //copying (0,0) -> (-1,-1) 
      /*
      PLOOP(iv)
	set_u(p,iv,-1,-1,iz,get_u(p,iv,0,0,iz));
      fill_geometry(-1,-1,iz,&geom);
      p2u(&get_u(p,0,-1,-1,iz),&get_u(u,0,-1,-1,iz),&geom);
      */

      //averaging <(-1,0),(0,-1)> -> (-1,-1)
      PLOOP(iv)
	set_u(p,iv,-1,-1,iz,.5*(get_u(p,iv,-1,0,iz)+get_u(p,iv,0,-1,iz)));
      fill_geometry(-1,-1,iz,&geom);
      p2u(&get_u(p,0,-1,-1,iz),&get_u(u,0,-1,-1,iz),&geom);

    }

  if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCHI)==1)
    {
 
      for(i=0;i<NG-1;i++)
	{
	  PLOOP(iv)
	  {
	    set_u(p,iv,-NG+i,NY,iz,get_u(p,iv,-NG+1,NY-1,iz));
	    set_u(p,iv,-1,NY+i+1,iz,get_u(p,iv,0,NY+i+1,iz));
	  }
	  fill_geometry(-NG+i,NY,iz,&geom);
	  p2u(&get_u(p,0,-NG+i,NY,iz),&get_u(u,0,-NG+i,NY,iz),&geom);
	  fill_geometry(-1,NY+i+1,iz,&geom);
	  p2u(&get_u(p,0,-1,NY+i+1,iz),&get_u(u,0,-1,NY+i+1,iz),&geom);
	}

      PLOOP(iv)
	set_u(p,iv,-1,NY,iz,.5*(get_u(p,iv,-1,NY-1,iz)+get_u(p,iv,0,NY,iz)));
      fill_geometry(-1,NY,iz,&geom);
      p2u(&get_u(p,0,-1,NY,iz),&get_u(u,0,-1,NY,iz),&geom);
   }

  if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCLO)==1)
    {
     for(i=0;i<NG-1;i++)
	{
	  PLOOP(iv)
	  {
	    set_u(p,iv,NX+i+1,-1,iz,get_u(p,iv,NX+i+1,0,iz));
	    set_u(p,iv,NX,-NG+i,iz,get_u(p,iv,NX-1,-NG+1,iz));
	  }
	  fill_geometry(NX+i+1,-1,iz,&geom);
	  p2u(&get_u(p,0,NX+i+1,-1,iz),&get_u(u,0,NX+i+1,-1,iz),&geom);
	  fill_geometry(NX,-NG+i,iz,&geom);
	  p2u(&get_u(p,0,NX,-NG+i,iz),&get_u(u,0,NX,-NG+i,iz),&geom);
	}

      PLOOP(iv)
	set_u(p,iv,NX,-1,iz,.5*(get_u(p,iv,NX-1,-1,iz)+get_u(p,iv,NX,0,iz)));
      fill_geometry(NX,-1,iz,&geom);
      p2u(&get_u(p,0,NX,-1,iz),&get_u(u,0,NX,-1,iz),&geom);

    }

 if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCHI)==1)
    {
      for(i=0;i<NG-1;i++)
	{
	  PLOOP(iv)
	  {
	    set_u(p,iv,NX+i+1,NY,iz,get_u(p,iv,NX+i+1,NY-1,iz));
	    set_u(p,iv,NX,NY+i+1,iz,get_u(p,iv,NX-1,NY+i+1,iz));
	  }
	  fill_geometry(NX+i+1,NY,iz,&geom);
	  p2u(&get_u(p,0,NX+i+1,NY,iz),&get_u(u,0,NX+i+1,NY,iz),&geom);
	  fill_geometry(NX,NY+i+1,iz,&geom);
	  p2u(&get_u(p,0,NX,NY+i+1,iz),&get_u(u,0,NX,NY+i+1,iz),&geom);
	}

      PLOOP(iv)
	set_u(p,iv,NX,NY,iz,.5*(get_u(p,iv,NX-1,NY,iz)+get_u(p,iv,NX,NY-1,iz)));
      fill_geometry(NX,NY,iz,&geom);
      p2u(&get_u(p,0,NX,NY,iz),&get_u(u,0,NX,NY,iz),&geom);
    }

 ldouble uval[NV],pval[NV];
 iz=0;

 //corners in the midda - apply boundary condition on what is already in ghost cells
 if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=-NG;j<0;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCLO);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(XBCLO)==1 && mpi_isitBC(YBCHI)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=NY;j<NY+NG;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCLO);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCLO)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=-NG;j<0;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCHI);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(XBCHI)==1 && mpi_isitBC(YBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=NY;j<NY+NG;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,XBCHI);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(YBCLO)==1 && mpi_isitBC(XBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=-NG;j<0;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCLO);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(YBCLO)==1 && mpi_isitBC(XBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=-NG;j<0;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCLO);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(YBCHI)==1 && mpi_isitBC(XBCLO)==0)
    {
      for(i=-NG;i<0;i++)
	for(j=NY;j<NY+NG;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCHI);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
 if(mpi_isitBC(YBCHI)==1 && mpi_isitBC(XBCHI)==0)
    {
      for(i=NX;i<NX+NG;i++)
	for(j=NY;j<NY+NG;j++)
	  {
	    set_bc_core(i,j,iz,t,uval,pval,ifinit,YBCHI);
	    PLOOP(iv)
	    {
	      set_u(u,iv,i,j,iz,uval[iv]);
	      set_u(p,iv,i,j,iz,pval[iv]);	      
	    }
	  }
    }
  
#endif

  return 0;
}

//fixing after failed MHD main inversion
int
cell_fixup_hd()
{
  if(DOFIXUPS==0)
    return 0;

  int ix,iy,iz,iv;
  int in,ii;
  int verbose=2;

  copy_u(1.,u,u_bak);
  copy_u(1.,p,p_bak);

  //gets the neiboring the primitives
#pragma omp parallel for private(ix,iy,iz,iv,ii,in) schedule (static)
  for(ix=0;ix<NX;ix++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {	      
	      if(get_cflag(HDFIXUPFLAG,ix,iy,iz)==1)
		{
		  //total fixups  
		  struct geometry geom;
		  fill_geometry(ix,iy,iz,&geom);

		  ldouble ppn[6][NV],pp[NV],uu[NV];

		  //int gix,giy,giz;
		  //mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);

		  //should care about global but at the stage where it is called knowns not about the boundaries

		  in=0; //number of successfull neighbors
		  
		  if(ix-1>=0 &&  get_cflag(HDFIXUPFLAG,ix-1,iy,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix-1,iy,iz);
		    }

		  if(ix+1<NX && get_cflag(HDFIXUPFLAG,ix+1,iy,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix+1,iy,iz);
		    }

		  if(iy-1>=0 && get_cflag(HDFIXUPFLAG,ix,iy-1,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy-1,iz);
		    }

		  if(iy+1<NY && get_cflag(HDFIXUPFLAG,ix,iy+1,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy+1,iz);
		    }

		  if(iz-1>=0 && get_cflag(HDFIXUPFLAG,ix,iy,iz-1)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy,iz-1);
		    }

		  if(iz+1<NZ && get_cflag(HDFIXUPFLAG,ix,iy,iz+1)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy,iz+1);
		    }

		  if((NZ==1 && NY==1 && in>=1) ||
		     (NZ==1 && in>=2) ||
		     (NY==1 && in>=1) ||
		     in>3) //sufficient number of neighbors
		    {
		      for(iv=0;iv<NV;iv++)
			{
			  pp[iv]=0;
			  for(ii=0;ii<in;ii++)
			    pp[iv]+=ppn[ii][iv];
			  pp[iv]/=(ldouble)in;  
			}
		      p2u(pp,uu,&geom);

		      if(verbose>1) 
			{
			  printf("%4d > %4d %4d %4d > MHDFIX > fixing up mhd with %d neighbors\n",PROCID,ix,iy,iz,in);
			}

		      //save to updated arrays memory
		      for(iv=0;iv<NVMHD;iv++)
			{
			  set_u(u_bak,iv,ix,iy,iz,uu[iv]);
			  set_u(p_bak,iv,ix,iy,iz,pp[iv]);
			}
		    }
		  else
		    {
		      #ifndef MPI
		      fprintf(fout_fail,"%4d > %4d %4d %4d > MHDFIXFAIL > didn't manage to hd fixup\n",PROCID,ix,iy,iz);
		      #endif
		      printf("%4d > %4d %4d %4d > MHDFIXFAIL > didn't manage to hd fixup \n",PROCID,ix,iy,iz);
		    }
		}
	    }
	}
    }

  //restoring to memory
  copy_u(1.,u_bak,u);
  copy_u(1.,p_bak,p);

  return 0;
}


//fixing after radiative implicit solver
int
cell_fixup_rad()
{
  if(DOFIXUPS==0)
    return 0;

  int ix,iy,iz,iv;
  int in,ii;
  int verbose=0;

  copy_u(1.,u,u_bak);
  copy_u(1.,p,p_bak);

  //gets the neighboring the primitives
#pragma omp parallel for private(ix,iy,iz,iv,ii,in) schedule (static)
  for(ix=0;ix<NX;ix++)
    {
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {	      
	      if(get_cflag(RADFIXUPFLAG,ix,iy,iz)<0)
		{
		  ldouble ppn[26][NV],pp[NV],uu[NV];

		  //total fixups  
		  struct geometry geom;
		  fill_geometry(ix,iy,iz,&geom); 
		  
		  //int gix,giy,giz;
		  //mpi_local2globalidx(ix,iy,iz,&gix,&giy,&giz);
		  //as above

		  in=0; //number of successfull neighbors
		  
		  if(ix-1>=0 &&  get_cflag(RADFIXUPFLAG,ix-1,iy,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix-1,iy,iz);
		    }
		  
		  if(ix+1<NX && get_cflag(RADFIXUPFLAG,ix+1,iy,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix+1,iy,iz);
		    }

		  if(iy-1>=0 && get_cflag(RADFIXUPFLAG,ix,iy-1,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy-1,iz);
		    }

		  if(iy+1<NY && get_cflag(RADFIXUPFLAG,ix,iy+1,iz)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy+1,iz);
		    }

		  if(iz-1>=0 && get_cflag(RADFIXUPFLAG,ix,iy,iz-1)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy,iz-1);
		    }

		  if(iz+1<NZ && get_cflag(RADFIXUPFLAG,ix,iy,iz+1)==0)
		    {
		      in++;
		      for(iv=0;iv<NV;iv++)
			ppn[in-1][iv]=get_u(p,iv,ix,iy,iz+1);
		    }

		  //try corners as well
		  //unnecessary + typos
		  //TODO: so far only for NZ==1
		  /*
		  if(NZ==1 && NY>1)
		    {
		      if(gix-1>=0 && giy-1>=0 && get_cflag(RADFIXUPFLAG,ix-1,iy-1,iz)==0)
			{
			  in++;
			  for(iv=0;iv<NV;iv++)
			    ppn[in-1][iv]=get_u(p,iv,ix-1,iy-1,iz);
			  idx[in-1][0]=ix-1;
			  idx[in-1][1]=iy-1;
			  idx[in-1][2]=iz;
			}

		      if(gix+1>=0 && giy-1>=0 && get_cflag(RADFIXUPFLAG,ix+1,iy-1,iz)==0)
			{
			  in++;
			  for(iv=0;iv<NV;iv++)
			    ppn[in-1][iv]=get_u(p,iv,ix+1,iy-1,iz);
			  idx[in-1][0]=ix+1;
			  idx[in-1][1]=iy-1;
			  idx[in-1][2]=iz;
			}

		      if(gix+1>=0 && giy+1>=0 && get_cflag(RADFIXUPFLAG,ix+1,iy+1,iz)==0)
			{
			  in++;
			  for(iv=0;iv<NV;iv++)
			    ppn[in-1][iv]=get_u(p,iv,ix+1,iy+1,iz);
			  idx[in-1][0]=ix+1;
			  idx[in-1][1]=iy+1;
			  idx[in-1][2]=iz;
			}

		      if(gix-1>=0 && giy+1>=0 && get_cflag(RADFIXUPFLAG,ix-1,iy+1,iz)==0)
			{
			  in++;
			  for(iv=0;iv<NV;iv++)
			    ppn[in-1][iv]=get_u(p,iv,ix-1,iy+1,iz);
			  idx[in-1][0]=ix-1;
			  idx[in-1][1]=iy+1;
			  idx[in-1][2]=iz;
			}			 
		    }
		  */
		  
		  if((NZ==1 && NY==1 && in>=1) ||
		     (NZ==1 && in>=1) ||
		     (NY==1 && in>=1) ||
		     in>3) //sufficient number of neighbors
		    {
		      for(iv=0;iv<NV;iv++)
			{
			  pp[iv]=0;
			  for(ii=0;ii<in;ii++)
			    pp[iv]+=ppn[ii][iv];
			  pp[iv]/=(ldouble)in;  
			}
		      p2u(pp,uu,&geom);
		      
		      if(verbose>1 || 1) 
			  printf("%4d > %4d %4d %4d > RADFIX > fixing up rad with %d neighbors\n",PROCID,ix,iy,iz,in);

		      //save to updated arrays memory
		      //all the primitives!
		      for(iv=0;iv<NV;iv++)
			{
			  //if(iv!=UU && iv!=EE0) continue; //why?
			  set_u(u_bak,iv,ix,iy,iz,uu[iv]);
			  set_u(p_bak,iv,ix,iy,iz,pp[iv]);
			}
		    }
		  else
		    {
                      #ifndef MPI
		      fprintf(fout_fail,"%4d > %4d %4d %4d > RADFIXFAIL > didn't manage to rad fixup\n",PROCID,ix,iy,iz);
                      #endif
		      printf("%4d > %4d %4d %4d > RADFIXFAIL > didn't manage to rad fixup \n",PROCID,ix,iy,iz);
		    }
		  
		}
	    }
	}
    }

  //restoring to memory
  copy_u(1.,u_bak,u);
  copy_u(1.,p_bak,p);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//outdated
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************


//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame, returns ultimate deltas ***********************
//******* the fiducial approach *****************************************
//**********************************************************************

struct f_implicit_metric_params
{
  void* ggg;
  ldouble *uu0;
  ldouble dt;
  int n;
};
     
int f_implicit_metric(const gsl_vector * x, void *paramsp,
	      gsl_vector * f)
//(ldouble *uu0,ldouble *uu,ldouble *pp,ldouble dt,void* ggg,ldouble *f)
{

 struct f_implicit_metric_params *params 
   = (struct f_implicit_metric_params *) paramsp;

  struct geometry *geom
    = (struct geometry *) params->ggg;

  ldouble dt=params->dt;
  
  ldouble *uu0=params->uu0;

  int n=params->n;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;

  ldouble pp[NV],uu[NV];
  int iv;
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=gsl_vector_get (x, iv);
      if(isnan(pp[iv])) return -1;
    }
  //  print_Nvector(pp,NV);
  
 //calculating conserved
  int corr,fixup[2];
  p2u(pp,uu,geom);

   ldouble ms[NV],fval[NV];
  //metric source terms
   //  print_Nvector(uu,NV);
   //  print_Nvector(uu0,NV);
  
  f_metric_source_term_arb(pp,geom,ms);

  for(iv=0;iv<n;iv++)
    {      
      fval[iv]=uu[iv] - uu0[iv] - ms[iv]*dt;
      gsl_vector_set (f, iv, fval[iv]);
    }

  //  print_Nvector(ms,NV);
  //  print_Nvector(fval,NV);
  //  getchar();

  return GSL_SUCCESS;
} 

int
print_state_metric (int iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .5e % .5e % .5e % .5e % .5e  % .5e "
	  "f(x) = % .3e % .3e % .3e % .3e % .3e % .3e\n",
	  iter,
	  gsl_vector_get (s->x, 0),
	  gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2),
	  gsl_vector_get (s->x, 3),
	  gsl_vector_get (s->x, 4),
	  gsl_vector_get (s->x, 5),

	  gsl_vector_get (s->f, 0),
	  gsl_vector_get (s->f, 1),
	  gsl_vector_get (s->f, 2),
	  gsl_vector_get (s->f, 3),
	  gsl_vector_get (s->f, 4),
	  gsl_vector_get (s->f, 5));

  return 0;
}


int
solve_implicit_metric(int ix,int iy,int iz,ldouble dt,ldouble *ubase)
{
  int i1,i2,i3,iv,i,j;

  int verbose=0;

  ldouble pp[NV],uu[NV],uu0[NV],uu00[NV],uup[NV]; 

  ldouble (*gg)[5],(*GG)[5];

  //  verbose=1;

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
      uu[iv]=get_u(ubase,iv,ix,iy,iz);  
      uu00[iv]=uu[iv];
      uu0[iv]=uu[iv];
   }

  const gsl_multiroot_fsolver_type *T;
  gsl_multiroot_fsolver *s;
     
  int status;
  int iter = 0;
     
  const size_t n = NV;
  struct f_implicit_metric_params par = {&geom,uu0,dt,n};
  gsl_multiroot_function f = {&f_implicit_metric, n, &par};

  gsl_vector *x = gsl_vector_alloc (n);
     

  double x_init[NV];
  for (iv=0;iv<n;iv++)
    {
      x_init[iv]=1.*pp[iv];     
      gsl_vector_set (x, iv, x_init[iv]);
    }
     
  T = gsl_multiroot_fsolver_hybrids;
  s = gsl_multiroot_fsolver_alloc (T, n);
  gsl_multiroot_fsolver_set (s, &f, x);

  //  print_state_metric (iter, s);  
  
  do
    {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);

      if (status)   /* check if solver is stuck */
	break;
      //print_state_metric (iter, s);
    
      status =gsl_multiroot_test_residual (s->f, 1e-5);

      status = gsl_multiroot_test_delta (s->dx, s->x,
					1e-4, 1e-4);
    }
  while (status == GSL_CONTINUE && iter < 1000);
     
  //printf ("status = %s\n", gsl_strerror (status)); getchar();
  if(status!=GSL_SUCCESS) {if(verbose) printf("complain : %s\n",gsl_strerror (status));return -1;}

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=gsl_vector_get (s->x, iv);
    }

  p2u(pp,uu,&geom);

  for(iv=0;iv<NV;iv++)
    {
      set_u(p,iv,ix,iy,iz,pp[iv]);      
      set_u(u,iv,ix,iy,iz,uu[iv]);      
    }
     
  gsl_multiroot_fsolver_free (s);
  gsl_vector_free (x);
  return 0;

}

//treats the most polar cells is a special way, correcting them, not evolving them
int
correct_polaraxis()
{
  int nc=NCCORRECTPOLAR; //correct velocity in nc most polar cells;

  int ix,iy,iz,iv,ic,iysrc,ixsrc;

  //spherical like coords
  if (MYCOORDS==SCHWCOORDS || MYCOORDS==KSCOORDS || MYCOORDS==KERRCOORDS || MYCOORDS==SPHCOORDS || MYCOORDS==MKS1COORDS || MYCOORDS==MKS2COORDS)
    {
#pragma omp parallel for private(ic,ix,iy,iz,iv,iysrc) schedule (static)
      for(ix=0;ix<NX;ix++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {
	      ldouble th,thsrc,thaxis;
	      ldouble pp[NV],uu[NV];
	      struct geometry geom;

	      //upper
	      if(TJ==0) //tile number
		{
		  thaxis=get_xb(0,1);
		  for(ic=0;ic<nc;ic++)
		    {
		      iy=ic;iysrc=nc;
		      th=get_x(iy,1);
		      thsrc=get_x(iysrc,1);	      
	      	  
		      fill_geometry(ix,iy,iz,&geom);
	  
		      PLOOP(iv)
			pp[iv]=get_u(p,iv,ix,iy,iz);
		  
		      //gas densities
		      pp[RHO]=get_u(p,RHO,ix,iysrc,iz);
		      pp[UU]=get_u(p,UU,ix,iysrc,iz);
		      pp[ENTR]=get_u(p,ENTR,ix,iysrc,iz);		  
		  		  		  
		      //gas velocities
		      pp[VX]=get_u(p,VX,ix,iysrc,iz);
		      pp[VZ]=get_u(p,VZ,ix,iysrc,iz);
		      pp[VY]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,VY,ix,iysrc,iz);

#ifdef MAGNFIELD
		      //gas velocities
		      pp[B1]=get_u(p,B1,ix,iysrc,iz);
		      pp[B3]=get_u(p,B3,ix,iysrc,iz);
		      pp[B2]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,B2,ix,iysrc,iz);
#endif

#ifdef RADIATION
		      //rad density
		      pp[EE0]=get_u(p,EE0,ix,iysrc,iz);

		      //rad velocities
		      pp[FX0]=get_u(p,FX0,ix,iysrc,iz);
		      pp[FZ0]=get_u(p,FZ0,ix,iysrc,iz);
		      pp[FY0]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,FY0,ix,iysrc,iz);
#endif 

		      p2u(pp,uu,&geom);

		      PLOOP(iv)
		      {
			set_u(p,iv,ix,iy,iz,pp[iv]);  
			set_u(u,iv,ix,iy,iz,uu[iv]);
		      }
		    }
		}
  
	      //lower
	      if(TJ==NTY-1)
		{
		  thaxis=get_xb(NY,1);
		  for(ic=0;ic<nc;ic++)
		    {
		      iy=NY-1-ic;iysrc=NY-1-nc;
		      th=get_x(iy,1);
		      thsrc=get_x(iysrc,1);	      
	      	  
		      fill_geometry(ix,iy,iz,&geom);
	  
		      PLOOP(iv)
			pp[iv]=get_u(p,iv,ix,iy,iz);
  
		      //gas densities
		      pp[RHO]=get_u(p,RHO,ix,iysrc,iz);
		      pp[UU]=get_u(p,UU,ix,iysrc,iz);
		      pp[ENTR]=get_u(p,ENTR,ix,iysrc,iz);		  

		      //gas velocities
		      pp[VX]=get_u(p,VX,ix,iysrc,iz);
		      pp[VZ]=get_u(p,VZ,ix,iysrc,iz);
		      pp[VY]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,VY,ix,iysrc,iz);

#ifdef MAGNFIELD
		      pp[B1]=get_u(p,B1,ix,iysrc,iz);
		      pp[B3]=get_u(p,B3,ix,iysrc,iz);
		      pp[B2]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,B2,ix,iysrc,iz);
#endif

#ifdef RADIATION
		      //rad density
		      pp[EE0]=get_u(p,EE0,ix,iysrc,iz);

		      //rad velocities
		      pp[FX0]=get_u(p,FX0,ix,iysrc,iz);
		      pp[FZ0]=get_u(p,FZ0,ix,iysrc,iz);
		      pp[FY0]=fabs((th-thaxis)/(thsrc-thaxis))*get_u(p,FY0,ix,iysrc,iz);
#endif 
		  
		      p2u(pp,uu,&geom);

		      PLOOP(iv)
		      {
			set_u(p,iv,ix,iy,iz,pp[iv]);  
			set_u(u,iv,ix,iy,iz,uu[iv]);
		      }
		    }
		}
	    }
	}
    }

  //cylindrical like coords
  if (MYCOORDS==CYLCOORDS || MYCOORDS==MCYL1COORDS)
    {
#pragma omp parallel for private(ic,ix,iy,iz,iv,ixsrc) schedule (static)
      for(iy=0;iy<NY;iy++)
	{
	  for(iz=0;iz<NZ;iz++)
	    {
	      ldouble R,Rsrc,Raxis;
	      ldouble pp[NV],uu[NV];
	      struct geometry geom;

	      //upper
	      Raxis=get_xb(0,0);
	      for(ic=0;ic<nc;ic++)
		{
		  ix=ic;ixsrc=nc;
		  R=get_x(ix,0);
		  Rsrc=get_x(ixsrc,0);	      
	      	  
		  fill_geometry(ix,iy,iz,&geom);
	  
		  PLOOP(iv)
		    pp[iv]=get_u(p,iv,ix,iy,iz);

		  //gas densities
		  pp[RHO]=get_u(p,RHO,ixsrc,iy,iz);
		  pp[UU]=get_u(p,UU,ixsrc,iy,iz);
		  pp[ENTR]=get_u(p,ENTR,ixsrc,iy,iz);		  

		  //gas velocities
		  pp[VY]=get_u(p,VY,ixsrc,iy,iz);
		  pp[VZ]=get_u(p,VZ,ixsrc,iy,iz);
		  pp[VX]=fabs((R-Raxis)/(Rsrc-Raxis))*get_u(p,VX,ixsrc,iy,iz);

		  #ifdef MAGNFIELD
		  pp[B2]=get_u(p,B2,ixsrc,iy,iz);
		  pp[B3]=get_u(p,B3,ixsrc,iy,iz);
		  pp[B1]=fabs((R-Raxis)/(Rsrc-Raxis))*get_u(p,B1,ixsrc,iy,iz);
		  #endif

		  #ifdef RADIATION
		  //rad density
		  pp[EE0]=get_u(p,EE0,ixsrc,iy,iz);

		  //rad velocities
		  pp[FY0]=get_u(p,FY0,ixsrc,iy,iz);
		  pp[FZ0]=get_u(p,FZ0,ixsrc,iy,iz);
		  pp[FX0]=fabs((R-Raxis)/(Rsrc-Raxis))*get_u(p,FX0,ixsrc,iy,iz);
		  #endif 

		  p2u(pp,uu,&geom);

		  PLOOP(iv)
		  {
		    set_u(p,iv,ix,iy,iz,pp[iv]);  
		    set_u(u,iv,ix,iy,iz,uu[iv]);
		  }
		}
	    }
	}
    }


  return 0; 
}

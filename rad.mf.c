//KORAL - rad.c
//radiation-related routines adjusted for multi rad fluids

#include "ko.h"

//***********************************************************************************
//******* redistributes radiation fluids wrapper ***************************************
//***********************************************************************************
int
redistribute_radfluids_at_cell(int ix,int iy,int iz)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  int iv;
  ldouble pp[NV],uu[NV];

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  if(MFREDISTRIBUTEMETHOD==1)
    redistribute_radfluids_m1(pp,uu,&geom);
  if(MFREDISTRIBUTEMETHOD==2)
    redistribute_radfluids_m2(pp,uu,&geom);
  if(MFREDISTRIBUTEMETHOD==3)
    redistribute_radfluids_m3(pp,uu,&geom);

  u2p_rad(uu,pp,&geom,&iv);

  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
      set_u(p,iv,ix,iy,iz,pp[iv]);
    }
	
  return 0;
}

//***********************************************************************************
//******* redistributes radiation fluids wrapper ***************************************
//***********************************************************************************
int
mf_correct_in_azimuth_at_cell(int ix,int iy,int iz,ldouble dt)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  int iv;
  ldouble pp[NV],uu[NV];

  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz);
      pp[iv]=get_u(p,iv,ix,iy,iz);
    }

  mf_correct_in_azimuth(pp,uu,&geom,dt);

  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
      set_u(p,iv,ix,iy,iz,pp[iv]);
    }
	
  return 0;
}

//***********************************************************************************
//******* corrects the distribution in azimuth by splitting highly azimuthal *******
//******* beams into almost-radial and more azimuthal one to avoid the *************
//******* the inner funnel *********************************************************
//***********************************************************************************
int
mf_correct_in_azimuth(ldouble *pp, ldouble *uu, void* ggg, ldouble dt)
{
  int NDIM=2;
  int verbose=0,ii,jj,irf;
  
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble radius,xxvec[4],xxvecCYL[4];
  get_xx(geom->ix,geom->iy,geom->iz,xxvec);
  coco_N(xxvec,xxvecCYL,MYCOORDS,CYLCOORDS);
  radius=xxvecCYL[1];

  //  if(pp[FX(0)]!=0.) verbose=1;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  //to ortonormal basis
  ldouble ppon[NV],ppon2[NV],pptemp[NV];
  prad_lab2on(pp,ppon,ggg);

  //total flux and energy
  ldouble Eon,Fon[3]={0.,0.,0.};
  Eon=0.;
  for(ii=0;ii<NRF;ii++)
    {
      Eon+=ppon[EE(ii)];
      Fon[0]+=ppon[FX(ii)];
      Fon[1]+=ppon[FY(ii)];
      Fon[2]+=ppon[FZ(ii)];
    }

  for(ii=0;ii<NVHD;ii++)
    {
      ppon2[ii]=ppon[ii];
      pptemp[ii]=ppon[ii];
    }

  for(ii=NVHD;ii<NV;ii++)
    ppon2[ii]=0.;

  if(verbose)
    {
      printf("\noooooooooo %d %d %d oooooooooo\n\n",geom->ix,geom->iy,geom->iz);
      printf("=== ppon ===\n");
      print_Nvector(ppon,NV);
    }

  for(irf=0;irf<NRF;irf++)
    {
      //within each fluid we check if the flux falls into the forbidden region
      //and we split onto fixed angle r-component along ksi

      ldouble ksi; 

      ksi=radius/10.*(M_PI/2./9.);//10 deg at r=10
      //ksi=(M_PI/2./9.);
      //      ksi=0.;

      ldouble EF[4]={ppon[EE(irf)],ppon[FX(irf)],ppon[FY(irf)],ppon[FZ(irf)]};
      ldouble f0=sqrt(EF[1]*EF[1]+EF[2]*EF[2]+EF[3]*EF[3])/EF[0];

      //TODO: generalize to 3d
      ldouble phi = atan2(EF[3],EF[1]); //angle on the r-phi plane
      
      if(phi<0.) phi+=2.*M_PI;

      if(verbose)
	{
	  printf("===\nfluid:  %d\n",irf);
	  printf("phi: %f ksi: %f\n",phi*180./M_PI,ksi*180./M_PI);
	}

      if(phi<=M_PI/2. || phi>=3.*M_PI/2. || f0<1.e-8 || fabs(phi-M_PI)<ksi)
	//flux pointing outwards, small or close to r-axis
	{
	  ppon2[EE(irf)]+=EF[0];
	  ppon2[FX(irf)]+=EF[1];
	  ppon2[FY(irf)]+=EF[2];
	  ppon2[FZ(irf)]+=EF[3];
	}
      else
	{
	  //inclined tetrad
	  ldouble erp[3]={-cos(ksi),0,sin(ksi)};
	  ldouble efp[3]={sin(ksi),0,cos(ksi)};

	  //components along the inclined tetrad
	  ldouble Frp = EF[1]*erp[0] + EF[3]*erp[2];
	  ldouble Ffp = EF[1]*efp[0] + EF[3]*efp[2];

	  //length of the original
	  ldouble Fzero = sqrt(EF[1]*EF[1]+EF[3]*EF[3]);
	  
	  //fraction applied:
	  ldouble frac;
	  //TODO: better estimate the velocity?
	  if(dt<0.)
	    frac=1.;
	  else
	    {
	      frac = dt / (radius / (1./3.)) * 10.;
	      if(frac>1.) frac=1.;
	    }


	  

	  if(verbose) printf("frac applied: %e\n",frac);

	  ppon2[EE(irf)]+=(1.-frac)*EF[0];
	  ppon2[FX(irf)]+=(1.-frac)*EF[1];
	  ppon2[FY(irf)]+=(1.-frac)*EF[2];
	  ppon2[FZ(irf)]+=(1.-frac)*EF[3];

	  //new component
	  ldouble En,Fn[3];
	  En=Frp/(Frp+Ffp)*EF[0];
	  Fn[0]=Frp*erp[0];
	  Fn[1]=Frp/(Frp+Ffp)*EF[2];
	  Fn[2]=Frp*erp[2];

	  if(verbose)
	    {
	      printf("component 1:\n %e %e %e %e\n",En,Fn[0],Fn[1],Fn[2]);
	    }

	  //distributing it over wedges
	  ldouble Avec[NRF];

	  //correcting for zero flux limit
	  ldouble f;
	  if(En<EEFLOOR) f=0.;
	  else
	  f=sqrt(Fn[0]*Fn[0]+Fn[1]*Fn[1]+Fn[2]*Fn[2])/En;
	  if(f<1.e-8)
	    {
	      for(ii=0;ii<NRF;ii++)
		Avec[ii]=1./(ldouble)NRF;
	    }
	  else
	    {      
	      assign_wedge_discrete_m1(Fn,Avec);	  
	      if(verbose) print_Nvector(Avec,NRF);
	      //transition from opticaly thin to thick
	      ldouble ftrans=0.0001;
	      ldouble Atrans=step_function(f-ftrans,ftrans/10.);
	      //to get rid of zeros
	      ldouble MINMIX=1.e-6;
	      if(Atrans>(1.-MINMIX)) Atrans=(1.-MINMIX);
	      for(ii=0;ii<NRF;ii++)
		Avec[ii]=Avec[ii]*Atrans + 1./(ldouble)NRF*(1.-Atrans);
	      //if(verbose) print_Nvector(Avec,NRF);
	    }	

	  for(ii=0;ii<NRF;ii++)
	    {
	      ppon2[EE(ii)]+=frac*Avec[ii]*En;
	      ppon2[FX(ii)]+=frac*Avec[ii]*Fn[0];
	      ppon2[FY(ii)]+=frac*Avec[ii]*Fn[1];
	      ppon2[FZ(ii)]+=frac*Avec[ii]*Fn[2];
	    }

	  //the other component
	  En=Ffp/(Frp+Ffp)*EF[0];
	  Fn[0]=Ffp*efp[0];
	  Fn[1]=Ffp/(Frp+Ffp)*EF[2];
	  Fn[2]=Ffp*efp[2];

	  if(verbose)
	    {
	      printf("component 2:\n %e %e %e %e\n",En,Fn[0],Fn[1],Fn[2]);
	    }
	 
	  //correcting for zero flux limit
	  if(En<EEFLOOR) f=0.; 
	  else
	  f=sqrt(Fn[0]*Fn[0]+Fn[1]*Fn[1]+Fn[2]*Fn[2])/En;
	  if(f<1.e-8)
	    {
	      for(ii=0;ii<NRF;ii++)
		Avec[ii]=1./(ldouble)NRF;
	    }
	  else
	    {      
	      assign_wedge_discrete_m1(Fn,Avec);
	      if(verbose) print_Nvector(Avec,NRF);
	      //transition from opticaly thin to thick
	      ldouble ftrans=0.0001;
	      ldouble Atrans=step_function(f-ftrans,ftrans/10.);
	      //to get rid of zeros
	      ldouble MINMIX=1.e-6;
	      if(Atrans>(1.-MINMIX)) Atrans=(1.-MINMIX);
	      for(ii=0;ii<NRF;ii++)
		Avec[ii]=Avec[ii]*Atrans + 1./(ldouble)NRF*(1.-Atrans);
	      //	      if(verbose) print_Nvector(Avec,NRF);
	    }
	  
	  for(ii=0;ii<NRF;ii++)
	    {
	      ppon2[EE(ii)]+=frac*Avec[ii]*En;
	      ppon2[FX(ii)]+=frac*Avec[ii]*Fn[0];
	      ppon2[FY(ii)]+=frac*Avec[ii]*Fn[1];
	      ppon2[FZ(ii)]+=frac*Avec[ii]*Fn[2];
	    }
	}
    }

  //in ppon2[] target ortonormal fluids
  if(verbose)
    {
      printf("=== ppon2 ===\n");
      print_Nvector(ppon2,NV);
      getchar();
    }

  //end total flux and energy
  ldouble Eon2,Fon2[3]={0.,0.,0.};
  Eon2=0.;
  for(ii=0;ii<NRF;ii++)
    {
      Eon2+=ppon2[EE(ii)];
      Fon2[0]+=ppon2[FX(ii)];
      Fon2[1]+=ppon2[FY(ii)];
      Fon2[2]+=ppon2[FZ(ii)];
    }

  if((fabs(1.-Eon2/Eon)>1.e-5) ||
     (Fon[0]!=0. && (fabs(1.-Fon2[0]/Fon[0])>1.e-5)) || 
     (Fon[1]!=0. && (fabs(1.-Fon2[1]/Fon[1])>1.e-5)) || 
     (Fon[2]!=0. && (fabs(1.-Fon2[2]/Fon[2])>1.e-5)))
    {
      printf("EF ratios: %f %f %f %f\n",Eon2/Eon,Fon2[0]/Fon[0],Fon2[1]/Fon[1],Fon2[2]/Fon[2]);
      getchar();
    }



  //to code basis
  prad_on2lab(ppon2,pp,ggg);

  //to conserved
  p2u_rad(pp,uu,gg,GG);

  return 0;
}



//***********************************************************************************
//******* assignes no. of wedge for discrete mixing basing on flux direction *******
//******* flux distributed always between adjacent wedges *******
//***********************************************************************************
int
assign_wedge_discrete_m1(ldouble F[], ldouble Avec[NRF])
{
  int NDIM=2; 
  int irf;

  if(NDIM==2)
    {
      ldouble factor=50.;
      ldouble nrf=(ldouble)NRF;
      ldouble phi,n0,n1,n,r,p;
      
      if(NZ==1)
	phi=atan2(F[1],F[0]);
      else if(NY==1)
	phi=atan2(F[2],F[0]);
      else my_err("assign_wedge_discrete_m1() not working for this setup\n");

      if(phi<0.) phi+=2.*M_PI;
      n0=(phi-2.*M_PI/2./nrf)/(2.*M_PI/nrf);
      n1=floor(n0);
      r=n0-n1;
      
      /*
      //pure power
      ldouble power=10.;
      if(r<0.5) 
	p=-pow(((1.-r)-.5)*2.,power)/2.;
      else
	p=pow((r-.5)*2.,power)/2.;
      */
      
      //exponent based smoothed step function
      r-=0.5;
      p=my_sign(r)*0.5*(exp(factor*fabs(r))-1.)/(exp(0.5*factor)-1.);

      n=n1+p+1.;

      //mixing coefficients
      for(irf=0;irf<NRF;irf++)
	Avec[irf]=0.;

      int i1,i2;
      i1=(int)floor(n);
      i2=(int)floor(n)+1;

      if(i1>NRF-1) i1-=NRF;
      if(i2>NRF-1) i2-=NRF;

      Avec[i1]=1.-(n-floor(n));
      Avec[i2]=(n-floor(n));
      
      /*
      print_Nvector(F,3);
      printf("%f %f %f %f %f\n",phi,n0,n1,n,nrf);      
      print_Nvector(Avec,NRF);
      */
    }

  return 0;
}

//***********************************************************************************
//******* assignes no. of wedge for discrete mixing basing on flux direction *******
//******* flux distributed basing on Lorentz boost of isotropic 1/3 along urf
//***********************************************************************************
int
assign_wedge_discrete_m2(ldouble flux[], ldouble Avec[NRF])
{
  int NDIM=2; 
  int irf;

  if(NDIM==2)
    {
      if(NZ!=1) my_err("assign_wedge_discrete() not working for 3D\n");

      ldouble nrf=(ldouble)NRF;
      ldouble phi[3],n0,n1,n,r,p;
      
      phi[0]=atan2(flux[1],flux[0]);
      if(phi[0]<0.) phi[0]+=2.*M_PI;
      
      //estimation of the radiative velocity
      ldouble urf=(flux[0]*flux[0]+flux[1]*flux[1]+flux[2]*flux[2]);
      urf=pow(urf,.05);
      
      ldouble vpar=(1./3.+urf)/(1.+1./3.*urf);
      ldouble vperp=sqrt(1./9.+urf*urf-1./9.*urf*urf);

      ldouble fan=atan(vperp/vpar);
      //fan=2.*M_PI/nrf/2.;
      
      phi[1]=phi[0]-fan;
      phi[2]=phi[0]+fan;

      /*
      if(phi[0]<0.) phi[0]+=2.*M_PI;
      if(phi[1]<0.) phi[1]+=2.*M_PI;
      if(phi[2]<0.) phi[2]+=2.*M_PI;

      if(phi[1]>phi[2]) {fan=phi[1];phi[1]=phi[2];phi[2]=fan;}
      */
      
      //wedges edges
      ldouble wedge[NRF+1];
      for(irf=1;irf<NRF+1;irf++)
	{
	  wedge[irf]=((ldouble)irf - 0.5) * 2.*M_PI/nrf;
	}
      wedge[0]=wedge[NRF]-2.*M_PI;
      // wedge[NRF]=wedge[0]+2.*M_PI;

      //print_Nvector(wedge,NRF+1);

      //mixing coefficients
      for(irf=0;irf<NRF;irf++)
	Avec[irf]=0.;

      //indices of wedges for phi[1-2]
      int phiedge[3]={-1,-1,-1};
      for(irf=-NRF;irf<2*NRF;irf++)
	{
	  //printf("%d %f\n",irf,((ldouble)(irf+1) - 0.5) * 2.*M_PI/nrf);
	  if(phi[1]>((ldouble)(irf+1) - 0.5) * 2.*M_PI/nrf)
	    phiedge[1]=irf+1;
	  if(phi[2]>((ldouble)(irf+1) - 0.5) * 2.*M_PI/nrf)
	    phiedge[2]=irf+1;
	}

      //printf("phi in wedges : %d %d\n",phiedge[1],phiedge[2]);

      int wno(int n)
      {
	while(n<0) n+=NRF;
	while(n>=NRF) n-=NRF;
	return n;
      }

      if(phiedge[2]-phiedge[1]>=2)
	{
	  for(irf=phiedge[1]+1;irf<=phiedge[2]-1;irf++)
	    Avec[wno(irf)]=1.;
	  Avec[wno(phiedge[1])]=((((ldouble)(phiedge[1]+1) - 0.5) * 2.*M_PI/nrf)-phi[1])/
	    (2.*M_PI/nrf);
	  Avec[wno(phiedge[2])]=(phi[2]-(((ldouble)(phiedge[2]+1-1)- 0.5) * 2.*M_PI/nrf))/
	    (2.*M_PI/nrf);
	  
	}
      if(phiedge[2]-phiedge[1]==1)
	{
	  Avec[wno(phiedge[1])]=((((ldouble)(phiedge[1]+1) - 0.5) * 2.*M_PI/nrf)-phi[1])/
	    (2.*M_PI/nrf);
	  Avec[wno(phiedge[2])]=(phi[2]-(((ldouble)(phiedge[2]+1-1) - 0.5) * 2.*M_PI/nrf))/
	    (2.*M_PI/nrf);
	}
      if(phiedge[2]-phiedge[1]==0)
	{
	  Avec[wno(phiedge[1])]=(phi[2]-phi[1])/
	    (2.*M_PI/nrf);	 
	}
    
      //print_Nvector(Avec,NRF);

      ldouble sum=0;
      for(irf=0;irf<NRF;irf++)
	sum+=Avec[irf];
       
      for(irf=0;irf<NRF;irf++)
	Avec[irf]=Avec[irf]/sum;

      /*
      print_Nvector(flux,3);
      printf("phi %f %f %f\n",phi[0],phi[1],phi[2]);      
      print_Nvector(Avec,NRF);
      */
    }

  return 0;
}

int
redistribute_radfluids(ldouble *pp, ldouble *uu, void* ggg)
{
  if(MFREDISTRIBUTEMETHOD==1)
    redistribute_radfluids_m1(pp,uu,ggg);
  if(MFREDISTRIBUTEMETHOD==2)
    redistribute_radfluids_m2(pp,uu,ggg);
  if(MFREDISTRIBUTEMETHOD==3)
    redistribute_radfluids_m3(pp,uu,ggg);
  return 0; 
}

//***********************************************************************************
//******* redistributes radiation fluids ***** discrete ********************************
//***********************************************************************************
int
redistribute_radfluids_m3(ldouble *pp, ldouble *uu0, void* ggg)
{
  int NDIM=2;
  int verbose=0,ii,jj,irf;
  ldouble A[NRF][NRF],uu1[NV];

  struct geometry *geom
   = (struct geometry *) ggg;

  //(geom->ix==IXDOT1-1 && geom->iy==IYDOT1-1) verbose=1;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  //to ortonormal basis
  ldouble ppon[NV],pp1[NV];
  prad_lab2on(pp,ppon,ggg);

  if(verbose)
    {
      printf("\noooooooooo %d %d %d oooooooooo\n\n",geom->ix,geom->iy,geom->iz);
      printf("=== uu0 ===\n");
      print_Nvector(uu0,NV);
      printf("=== ppon ===\n");
      print_Nvector(ppon,NV);
    }


  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      A[ii][jj]=0.;

  ldouble flux[3],f;
  for(irf=0;irf<NRF;irf++)
    {
      flux[0]=ppon[FX(irf)]/ppon[EE(irf)];
      flux[1]=ppon[FY(irf)]/ppon[EE(irf)];
      flux[2]=ppon[FZ(irf)]/ppon[EE(irf)];
      //f=|F|/cE
      f=sqrt(flux[0]*flux[0]+flux[1]*flux[1]+flux[2]*flux[2]);

      if(f<1.e-8)
	{
	  for(ii=0;ii<NRF;ii++)
	    A[irf][ii]=1./(ldouble)NRF;
	}
      else
	{      
	  //coefficients telling where irf fluid should be moved
	  ldouble Avec[NRF];
	  assign_wedge_discrete_m1(flux,Avec);

	  /*
	  ldouble phi;	  
	  for(phi=0.;phi<2.*M_PI;phi+=M_PI/30.)
	    {
	    flux[0]=1.*cos(phi);
	    flux[1]=1.*sin(phi);
	    flux[2]=0.;
	    assign_wedge_discrete_m2(flux,Avec);getchar();
	    }
	  */
	    

	  //transition from opticaly thin to thick
	  ldouble ftrans=0.1;
	  ldouble Atrans=step_function(f-ftrans,ftrans/10.);

	  //to get rid of zeros
	  ldouble MINMIX=1.e-8;
	  if(Atrans>(1.-MINMIX)) Atrans=(1.-MINMIX);

	  if(verbose) 
	    {
	      print_Nvector(flux,3);
	      printf("=== coeff %d ===\n\n",irf);
	      print_Nvector(Avec,NRF);
	      printf("Atrans: %e\n",Atrans);
	    }
			       
	  for(ii=0;ii<NRF;ii++)
	    A[irf][ii]=Avec[ii]*Atrans + 1./(ldouble)NRF*(1.-Atrans);

	}
    }

  
  for(ii=0;ii<NVHD;ii++)
    uu1[ii]=uu0[ii];
  for(ii=NVHD;ii<NV;ii++)
    uu1[ii]=0.;

  if(verbose) 
    {
      printf("=== coefficients ===\n");
    }
  
  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      {
	if(verbose)
	  printf(" %d -> %d : %e\n",jj,ii,A[jj][ii]);

	uu1[EE(ii)]+=uu0[EE(jj)]*A[jj][ii];
	uu1[FX(ii)]+=uu0[FX(jj)]*A[jj][ii];
	uu1[FY(ii)]+=uu0[FY(jj)]*A[jj][ii];
	uu1[FZ(ii)]+=uu0[FZ(jj)]*A[jj][ii];
      }

  if(verbose) 
    {
      printf("=== uu1 ===\n");
      print_Nvector(uu1,NV);
      getchar();
    }

  for(ii=NVHD;ii<NV;ii++)
    uu0[ii]=uu1[ii];

  return 0;
}
			     

//***********************************************************************************
//******* redistributes radiation fluids ***********************************************
//***********************************************************************************
int
redistribute_radfluids_m1(ldouble *pp, ldouble *uu0, void* ggg)
{
  int NDIM=2;
  int method=4;
  ldouble power=10.;
  int verbose=0;
  
  struct geometry *geom
   = (struct geometry *) ggg;

  //if(geom->ix==IXDOT1+1 && geom->iy==IYDOT1+2 ) verbose=1;

  if(verbose)
    {
      printf("\noooooooooo %d %d %d oooooooooo\n\n",geom->ix,geom->iy,geom->iz);
      printf("=== uu0 ===\n");
      print_Nvector(uu0,NV);
      printf("=== pp0 ===\n");
      print_Nvector(pp,NV);
    }

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  //calculates wavespeed for each of the fluids
  ldouble aval[NRF][6];

  int irf,ii,jj;
  ldouble uu1[NV],A[NRF][NRF];

  calc_rad_wavespeeds_pure_mf_each(pp,geom,aval);

  for(ii=0;ii<NRF;ii++)
    {
      aval[ii][0]*=sqrt(gg[1][1]);
      aval[ii][1]*=sqrt(gg[1][1]);
      aval[ii][2]*=sqrt(gg[2][2]);
      aval[ii][3]*=sqrt(gg[2][2]);
      aval[ii][4]*=sqrt(gg[3][3]);
      aval[ii][5]*=sqrt(gg[3][3]);
    }

  if(verbose)
    {
      printf("=== wavespeeds ===\n");
      for(ii=0;ii<NRF;ii++)
	printf("%d : [%e %e] [%e %e] [%e %e]\n",ii,aval[ii][0],aval[ii][1],aval[ii][2],aval[ii][3],aval[ii][4],aval[ii][5]);
    }

  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      A[ii][jj]=0.;

  ldouble MINVEL=1.e-3;

  for(irf=0;irf<NRF;irf++)
    {
      //in aval[NRF][6] one has rad.char.wavespeeds for [irf] fluid
      //aval[irf][0] - min left going in x
      //aval[irf][1] - max right going in x
      //aval[irf][2] - min left going in y 
      //etc...

      if(NDIM==1)
	{
	  ldouble vxl,vxr;
	  vxl=my_min(aval[irf][0],-MINVEL);
	  vxr=my_max(aval[irf][1],MINVEL);
	  
	  if(method==0)
	    {
	      //mixing linear in characteristic velocities
	      vxl=fabs(vxl);
	      A[irf][0]=vxr/(vxl+vxr);
	      A[irf][1]=vxl/(vxl+vxr);
	    }

	  if(method==1)
	    {
	      //mixing to arb. power in characteristic velocities	  
	      vxl=fabs(vxl);
	      A[irf][0]=pow(vxr,power)/(pow(vxl,power)+pow(vxr,power));
	      A[irf][1]=pow(vxl,power)/(pow(vxl,power)+pow(vxr,power));
	    }
	  
	  if(method==2)
	    {
	      //discrete mixing	      
	      if(fabs((vxr-vxl)/vxr)<1.e-5)
		{
		  A[irf][0]=.5;
		  A[irf][1]=.5;
		}
	      else if(vxr>vxl)
		{
		  A[irf][0]=1.-MINVEL;
		  A[irf][1]=0.+MINVEL;
		}
	      else
		{
		  A[irf][1]=1.-MINVEL;
		  A[irf][0]=0.+MINVEL;
		}
	    }
	  
	}

      if(NDIM==2)
	{
	  ldouble vxl,vxr,vyl,vyr;
	  if(method==0 || method==1)
	    //direct vchar
	    {
	      if(NZ==1)
		{
		  vxl=my_min(aval[irf][0],-MINVEL);
		  vxr=my_max(aval[irf][1],MINVEL);
		  vyl=my_min(aval[irf][2],-MINVEL);
		  vyr=my_max(aval[irf][3],MINVEL);
		}
	      else if(NY==1)
		{
		  vxl=my_min(aval[irf][0],-MINVEL);
		  vxr=my_max(aval[irf][1],MINVEL);
		  vyl=my_min(aval[irf][4],-MINVEL);
		  vyr=my_max(aval[irf][5],MINVEL);
		}
	      vxl=fabs(vxl);
	      vyl=fabs(vyl);	      
	    }

	  if(method==3 || method==4)
	    //diagonal velocities
	    {
	      ldouble vxl0,vxr0,vyl0,vyr0;
	      if(NZ==1)
		{
		  vxl0=my_min(aval[irf][0],-MINVEL);
		  vxr0=my_max(aval[irf][1],MINVEL);
		  vyl0=my_min(aval[irf][2],-MINVEL);
		  vyr0=my_max(aval[irf][3],MINVEL);
		}
	      else if(NY==1)
		{
		  vxl0=my_min(aval[irf][0],-MINVEL);
		  vxr0=my_max(aval[irf][1],MINVEL);
		  vyl0=my_min(aval[irf][4],-MINVEL);
		  vyr0=my_max(aval[irf][5],MINVEL);
		}
	      
	      vxr=sqrt(vxr0*vxr0+vyl0*vyl0);
	      vyr=sqrt(vxr0*vxr0+vyr0*vyr0);
	      vxl=sqrt(vxl0*vxl0+vyr0*vyr0);
	      vyl=sqrt(vxl0*vxl0+vyl0*vyl0);
	    }
	    
	  if(method==0 || method==3)
	    {
	      //mixing linear in characteristic velocities
	      
	      A[irf][0]=vxr/(vxl+vxr)*vyr/(vyl+vyr);
	      A[irf][1]=vxl/(vxl+vxr)*vyr/(vyl+vyr);
	      A[irf][2]=vxl/(vxl+vxr)*vyl/(vyl+vyr);
	      A[irf][3]=vxr/(vxl+vxr)*vyl/(vyl+vyr);
	    }

	  if(method==1 || method==4)
	    {
	      //arbitrary power
	      A[irf][0]=pow(vxr,power)/(pow(vxl,power)+pow(vxr,power))*pow(vyr,power)/(pow(vyl,power)+pow(vyr,power));
	      A[irf][1]=pow(vxl,power)/(pow(vxl,power)+pow(vxr,power))*pow(vyr,power)/(pow(vyl,power)+pow(vyr,power));
	      A[irf][2]=pow(vxl,power)/(pow(vxl,power)+pow(vxr,power))*pow(vyl,power)/(pow(vyl,power)+pow(vyr,power));
	      A[irf][3]=pow(vxr,power)/(pow(vxl,power)+pow(vxr,power))*pow(vyl,power)/(pow(vyl,power)+pow(vyr,power));
	    }

	  if(method==2)
	    {
	      ldouble DUMPEDGE=0.001;
	      ldouble MINMIXING=1.e-5;
		  
	      if(NZ==1)
		{
		  //to ortonormal basis
		  ldouble pp0[NV];
		  prad_lab2on(pp,pp0,ggg);
		  ldouble fdump=(pp0[FX(irf)]*pp0[FX(irf)]+pp0[FY(irf)]*pp0[FY(irf)]+pp0[FZ(irf)]*pp0[FZ(irf)])/pp0[EE(irf)]/pp0[EE(irf)];

		  ldouble dumping=exp(-fdump/DUMPEDGE);

		  int wedgeno=-1;
		  if(pp0[FX(irf)]>fabs(pp0[FY(irf)]))
		    wedgeno=0;
		  if(pp0[FX(irf)]<-fabs(pp0[FY(irf)]))
		    wedgeno=2;
		  if(pp0[FY(irf)]>fabs(pp0[FX(irf)]))
		    wedgeno=1;
		  if(pp0[FY(irf)]<-fabs(pp0[FX(irf)]))
		    wedgeno=3;

		  if(wedgeno>=0) //not aligned with axes
		    {
		      A[irf][wedgeno]=1.-dumping*3./4.;
		      for(ii=0;ii<NRF;ii++)
			{
			  if(ii==wedgeno) continue;
			  A[irf][ii]=dumping*1./4.;
			  if(A[irf][ii]<MINMIXING) {
			    A[irf][ii]=MINMIXING;
			    A[irf][wedgeno]-=MINMIXING;
			  }
			}
		    }
		  else //special handling of aligned fluxes - not sure if necessary
		    {
		      for(ii=0;ii<NRF;ii++)
			A[irf][ii]=MINMIXING;
		      A[irf][irf]=1.-3.*MINMIXING;
		    }	   

		  if(verbose) printf("=== dumping for irf=%d -> %f\n",irf,dumping);
		}
	      if(NY==1)
		{
		  //to ortonormal basis
		  ldouble pp0[NV];
		  prad_lab2on(pp,pp0,ggg);
		  ldouble fdump=(pp0[FX(irf)]*pp0[FX(irf)]+pp0[FY(irf)]*pp0[FY(irf)]+pp0[FZ(irf)]*pp0[FZ(irf)])/pp0[EE(irf)]/pp0[EE(irf)];

		  ldouble dumping=exp(-fdump/DUMPEDGE);

		  int wedgeno=-1;
		  if(pp0[FX(irf)]>fabs(pp0[FZ(irf)]))
		    wedgeno=0;
		  if(pp0[FX(irf)]<-fabs(pp0[FZ(irf)]))
		    wedgeno=2;
		  if(pp0[FZ(irf)]>fabs(pp0[FX(irf)]))
		    wedgeno=1;
		  if(pp0[FZ(irf)]<-fabs(pp0[FX(irf)]))
		    wedgeno=3;

		  if(wedgeno>=0) //not aligned with axes
		    {
		      A[irf][wedgeno]=1.-dumping*3./4.;
		      for(ii=0;ii<NRF;ii++)
			{
			  if(ii==wedgeno) continue;
			  A[irf][ii]=dumping*1./4.;
			  if(A[irf][ii]<MINMIXING) {
			    A[irf][ii]=MINMIXING;
			    A[irf][wedgeno]-=MINMIXING;
			  }
			}
		    }
		  else //special handling of aligned fluxes - not sure if necessary
		    {
		      for(ii=0;ii<NRF;ii++)
			A[irf][ii]=MINMIXING;
		      A[irf][irf]=1.-3.*MINMIXING;
		    }	   

		  if(verbose) printf("=== dumping for irf=%d -> %f\n",irf,dumping);
		}
	    }
	}

      if(NDIM==3)
	{
	  my_err("NDIM==3 not implemented in redistribute()\n");
	}
    }  

  for(ii=0;ii<NVHD;ii++)
    uu1[ii]=uu0[ii];
  for(ii=NVHD;ii<NV;ii++)
    uu1[ii]=0.;

  if(verbose) 
    {
      printf("=== coefficients ===\n");
    }
  
  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      {
	if(verbose)
	  printf(" %d -> %d : %e\n",jj,ii,A[jj][ii]);

	uu1[EE(ii)]+=uu0[EE(jj)]*A[jj][ii];
	uu1[FX(ii)]+=uu0[FX(jj)]*A[jj][ii];
	uu1[FY(ii)]+=uu0[FY(jj)]*A[jj][ii];
	uu1[FZ(ii)]+=uu0[FZ(jj)]*A[jj][ii];
      }

  if(verbose) 
    {
      printf("=== uu1 ===\n");
      print_Nvector(uu1,NV);
      getchar();
    }

  for(ii=NVHD;ii<NV;ii++)
    uu0[ii]=uu1[ii];
  
  return 0;
}


//***********************************************************************************
//******* redistributes radiation fluids ***********************************************
//******* basing on velocities along axes **********************************************
//***********************************************************************************
int
redistribute_radfluids_m2(ldouble *pp, ldouble *uu0, void* ggg)
{
#ifdef MULTIRADFLUID

  ldouble skew=10.;
  int verbose=0;
  
  struct geometry *geom
   = (struct geometry *) ggg;

  //if(geom->ix==IXDOT1+1 && geom->iy==IYDOT1   ) verbose=1;

  if(verbose)
    {
      printf("\noooooooooo %d %d %d oooooooooo\n\n",geom->ix,geom->iy,geom->iz);
      printf("=== uu0 ===\n");
      print_Nvector(uu0,NV);
      printf("=== pp0 ===\n");
      print_Nvector(pp,NV);
    }

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  //calculates wavespeed for each of the fluids
  ldouble aval[NRF][6];

  int irf,ii,jj;
  ldouble uu1[NV],A[NRF][NRF];

  calc_rad_wavespeeds_pure_mf_each(pp,geom,aval);

  //to ortonormal - too simple?
  for(ii=0;ii<NRF;ii++)
    {
      aval[ii][0]*=sqrt(gg[1][1]);
      aval[ii][1]*=sqrt(gg[1][1]);
      aval[ii][2]*=sqrt(gg[2][2]);
      aval[ii][3]*=sqrt(gg[2][2]);
      aval[ii][4]*=sqrt(gg[3][3]);
      aval[ii][5]*=sqrt(gg[3][3]);
    }

  if(verbose)
    {
      printf("=== wavespeeds ===\n");
      for(ii=0;ii<NRF;ii++)
	printf("%d : [%e %e] [%e %e] [%e %e]\n",ii,aval[ii][0],aval[ii][1],aval[ii][2],aval[ii][3],aval[ii][4],aval[ii][5]);
    }

  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      A[ii][jj]=0.;

  ldouble MINVEL=1.e-4;

  for(irf=0;irf<NRF;irf++)
    {
      //in aval[NRF][6] one has rad.char.wavespeeds for [irf] fluid
      //aval[irf][0] - min left going in x
      //aval[irf][1] - max right going in x
      //aval[irf][2] - min left going in y 
      //etc...

      ldouble vxl,vxr,vyl,vyr,vzr,vzl,sumvel;
      vxl=fabs(my_min(aval[irf][0],-MINVEL));
      vxr=fabs(my_max(aval[irf][1],MINVEL));
      vyl=fabs(my_min(aval[irf][2],-MINVEL));
      vyr=fabs(my_max(aval[irf][3],MINVEL));
      vzl=fabs(my_min(aval[irf][4],-MINVEL));
      vzr=fabs(my_max(aval[irf][5],MINVEL));
    
      //wedges (x,y,z)
      //0 (x+)
      //1 (x-)
      //2 (y+)
      //3 (y-)
      //4 (z+)
      //5 (z-)

      /*
      vxl=pow(fabs(vxl),power);
      vxr=pow(fabs(vxr),power);
      vyl=pow(fabs(vyl),power);
      vyr=pow(fabs(vyr),power);
      vzl=pow(fabs(vzl),power);
      vzr=pow(fabs(vzr),power);
      */
      
      ldouble expskew(ldouble vel, ldouble skew) 
      {
	return (exp(skew*vel)-1.0)/(exp(skew)-1.0);
      }

      vxl=expskew(vxl,skew);
      vxr=expskew(vxr,skew);
      vyl=expskew(vyl,skew);
      vyr=expskew(vyr,skew);
      vzl=expskew(vzl,skew);
      vzr=expskew(vzr,skew);

      sumvel=vxl+vxr+vyl+vyr+vzl+vzr;

      if(NY==1)
	sumvel=vxl+vxr+vzl+vzr;
      if(NZ==1)
	sumvel=vxl+vxr+vyl+vyr;
      if(NY==1 && NZ==1)
	sumvel=vxl+vxr+2.*expskew(MINVEL,skew);

      //arbitrary power
      A[irf][0]=vxr/sumvel;
      A[irf][1]=vxl/sumvel;
      A[irf][2]=vyr/sumvel;
      A[irf][3]=vyl/sumvel;

      if(NY==1)
	{
	  A[irf][2]=vzr/sumvel;
	  A[irf][3]=vzl/sumvel;
	}

      if(NY==1 && NZ==1)
	{
	  A[irf][2]=expskew(MINVEL,skew)/sumvel;
	  A[irf][3]=expskew(MINVEL,skew)/sumvel;
	}

      if(NRF==6)
	{
	  A[irf][4]=vzr/sumvel;
	  A[irf][5]=vzl/sumvel;
	}

    }  

  for(ii=0;ii<NVHD;ii++)
    uu1[ii]=uu0[ii];
  for(ii=NVHD;ii<NV;ii++)
    uu1[ii]=0.;

  if(verbose) 
    {
      printf("=== coefficients ===\n");
    }
  
  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      {
	if(verbose)
	  printf(" %d -> %d : %e\n",jj,ii,A[jj][ii]);

	uu1[EE(ii)]+=uu0[EE(jj)]*A[jj][ii];
	uu1[FX(ii)]+=uu0[FX(jj)]*A[jj][ii];
	uu1[FY(ii)]+=uu0[FY(jj)]*A[jj][ii];
	uu1[FZ(ii)]+=uu0[FZ(jj)]*A[jj][ii];
      }

  if(verbose) 
    {
      printf("=== uu1 ===\n");
      print_Nvector(uu1,NV);
      getchar();
    }

  for(ii=NVHD;ii<NV;ii++)
    uu0[ii]=uu1[ii];
  
#endif

  return 0;
}

//***********************************************************************************
//******* redistributes radiation fluids projecting on the axes *********************
//******* should not be used throughout the simulation only for initial cond ********
//***********************************************************************************
int
redistribute_radfluids_along_axes(ldouble *pp, ldouble *uu, void* ggg)
{
  int verbose=0;
  
  struct geometry *geom
   = (struct geometry *) ggg;

  //  if(geom->ix==IXDOT1+1 && geom->iy==IYDOT1+1) verbose=1;

  if(verbose)
    {
      printf("\noooooooooo %d %d %d oooooooooo\n\n",geom->ix,geom->iy,geom->iz);
      printf("=== uu ===\n");
      print_Nvector(uu,NV);
      printf("=== pp ===\n");
      print_Nvector(pp,NV);
    }

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  //to ortonormal basis
  ldouble pp0[NV],pp1[NV];
  prad_lab2on(pp,pp0,ggg);

  if(verbose)
    {
      printf("=== pp on ===\n");
      print_Nvector(pp0,NV);
    }

  int irf,ii,jj,i2;
  ldouble Fmag2,Etot,phi,FXtot,FYtot,FZtot;
  ldouble FRAC=1.;
  for(ii=0;ii<NV;ii++)
    pp1[ii]=pp0[ii];

  Etot=pp0[EE(0)]+pp0[EE(1)]+pp0[EE(2)]+pp0[EE(3)];
  FXtot=pp0[FX(0)]+pp0[FX(1)]+pp0[FX(2)]+pp0[FX(3)];
  FYtot=pp0[FY(0)]+pp0[FY(1)]+pp0[FY(2)]+pp0[FY(3)];
  FZtot=pp0[FZ(0)]+pp0[FZ(1)]+pp0[FZ(2)]+pp0[FZ(3)];

  if(verbose) printf("Tots: %e %e %e %e\n",Etot,FXtot,FYtot,FZtot);
  
  for(irf=0;irf<NRF;irf++)
    {
      if(NZ==1)
	{
	  //order corresponding to method=2 from above, i.e. not to quadrants but wedges
	  //irf = 0 - x+ axis
	  //irf = 1 - y+ axis
	  //irf = 2 - x- axis
	  //irf = 3 - y- axis

	  Fmag2=(pp0[FX(irf)]*pp0[FX(irf)]+pp0[FY(irf)]*pp0[FY(irf)]);

	  //FX
	  i2=-1;
	  if((irf==1 || irf==3 || irf==2) && pp0[FX(irf)]>0.) i2=0;
	  if((irf==1 || irf==3 || irf==0) && pp0[FX(irf)]<0.) i2=2;
	  
	  if(i2>=0)
	    {
	      pp1[FX(i2)]+=FRAC*pp0[FX(irf)];
	      pp1[EE(i2)]+=FRAC*pp0[EE(irf)]*pp0[FX(irf)]*pp0[FX(irf)]/Fmag2;
	      pp1[FX(irf)]-=FRAC*pp0[FX(irf)];
	      pp1[EE(irf)]-=FRAC*pp0[EE(irf)]*pp0[FX(irf)]*pp0[FX(irf)]/Fmag2;
	    }

	  //FY
	  i2=-1;
	  if((irf==0 || irf==2 || irf==3) && pp0[FY(irf)]>0.) i2=1;
	  if((irf==0 || irf==2 || irf==1) && pp0[FY(irf)]<0.) i2=3;

	  if(i2>=0)
	    {
	      pp1[FY(i2)]+=FRAC*pp0[FY(irf)];
	      pp1[EE(i2)]+=FRAC*pp0[EE(irf)]*pp0[FY(irf)]*pp0[FY(irf)]/Fmag2;
	      pp1[FY(irf)]-=FRAC*pp0[FY(irf)];
	      pp1[EE(irf)]-=FRAC*pp0[EE(irf)]*pp0[FY(irf)]*pp0[FY(irf)]/Fmag2;
	    }	
	}

      if(NY==1)
	{
	  //order corresponding to method=2 from above, i.e. not to quadrants but wedges
	  //irf = 0 - x+ axis
	  //irf = 1 - z+ axis
	  //irf = 2 - x- axis
	  //irf = 3 - z- axis

	  Fmag2=(pp0[FX(irf)]*pp0[FX(irf)]+pp0[FZ(irf)]*pp0[FZ(irf)]);

	  //FX
	  i2=-1;
	  if((irf==1 || irf==3 || irf==2) && pp0[FX(irf)]>0.) i2=0;
	  if((irf==1 || irf==3 || irf==0) && pp0[FX(irf)]<0.) i2=2;
	  
	  if(i2>=0)
	    {
	      pp1[FX(i2)]+=FRAC*pp0[FX(irf)];
	      pp1[EE(i2)]+=FRAC*pp0[EE(irf)]*pp0[FX(irf)]*pp0[FX(irf)]/Fmag2;
	      pp1[FX(irf)]-=FRAC*pp0[FX(irf)];
	      pp1[EE(irf)]-=FRAC*pp0[EE(irf)]*pp0[FX(irf)]*pp0[FX(irf)]/Fmag2;
	    }

	  //FZ
	  i2=-1;
	  if((irf==0 || irf==2 || irf==3) && pp0[FZ(irf)]>0.) i2=1;
	  if((irf==0 || irf==2 || irf==1) && pp0[FZ(irf)]<0.) i2=3;

	  if(i2>=0)
	    {
	      pp1[FZ(i2)]+=FRAC*pp0[FZ(irf)];
	      pp1[EE(i2)]+=FRAC*pp0[EE(irf)]*pp0[FZ(irf)]*pp0[FZ(irf)]/Fmag2;
	      pp1[FZ(irf)]-=FRAC*pp0[FZ(irf)];
	      pp1[EE(irf)]-=FRAC*pp0[EE(irf)]*pp0[FZ(irf)]*pp0[FZ(irf)]/Fmag2;
	    }	
	}
    }

  //redistribution to make it more uniform
  ldouble REDISTR=.2;
  ldouble MINCONTRAST=1.e-7;
  ldouble invsum[4]={0.,0.,0.,0.};
  ldouble dEE,dFX,dFY,dFZ;
  dEE=dFX=dFY=dFZ=0.;
 if(verbose) 
    {
      printf("=== pp 0 ===\n");
      print_Nvector(pp1,NV);      
    }
  for(ii=NVHD;ii<NV;ii++)
    pp0[ii]=pp1[ii];

  for(irf=0;irf<NRF;irf++)
    {
      dEE+=REDISTR*pp1[EE(irf)];
      dFX+=REDISTR*pp1[FX(irf)];
      dFY+=REDISTR*pp1[FY(irf)];
      dFZ+=REDISTR*pp1[FZ(irf)];
      pp1[EE(irf)]-=REDISTR*pp1[EE(irf)];
      pp1[FX(irf)]-=REDISTR*pp1[FX(irf)];
      pp1[FY(irf)]-=REDISTR*pp1[FY(irf)];
      pp1[FZ(irf)]-=REDISTR*pp1[FZ(irf)];
      invsum[0]+=1./(pp1[EE(irf)]+MINCONTRAST*Etot);
      invsum[1]+=1./(pp1[FX(irf)]+MINCONTRAST*Etot);
      invsum[2]+=1./(pp1[FY(irf)]+MINCONTRAST*Etot);
      invsum[3]+=1./(pp1[FZ(irf)]+MINCONTRAST*Etot);
    }
 if(verbose) 
    {
      printf("=== pp 1 ===\n");
      print_Nvector(pp1,NV);      
    }

  for(irf=0;irf<NRF;irf++)
    {
      /*
      ldouble frdEE=(dEE/invsum/(pp1[EE(irf)]+MINCONTRAST*Etot)-pp0[EE(irf)])/dEE;
      pp1[FX(irf)]+=dFX * frdEE;
      pp1[FY(irf)]+=dFY * frdEE;
      pp1[FZ(irf)]+=dFZ * frdEE;
      */
      pp1[EE(irf)]+=dEE/invsum[0]/(pp1[EE(irf)]+MINCONTRAST*Etot);
      pp1[FX(irf)]+=dFX/invsum[1]/(pp1[FX(irf)]+MINCONTRAST*Etot);
      pp1[FY(irf)]+=dFY/invsum[2]/(pp1[FY(irf)]+MINCONTRAST*Etot);
      pp1[FZ(irf)]+=dFZ/invsum[3]/(pp1[FZ(irf)]+MINCONTRAST*Etot);
    }

 if(verbose) 
    {
      printf("=== pp 2 ===\n");
      print_Nvector(pp1,NV);     
    }

  for(ii=NVHD;ii<NV;ii++)
    pp0[ii]=pp1[ii];

  /*
  int i,j;
  for(i=0;i<NRF;i++)
    {
      pp0[FX(irf)]=pp0[FY(irf)]=pp0[FZ(irf)]=0.;
      for(j=0;j<NRF;j++)
	{	  
	  pp0[FX(i)]+=pp1[EE(i)]/Etot*pp1[FX(j)];
	  pp0[FY(i)]+=pp1[EE(i)]/Etot*pp1[FY(j)];
	  pp0[FZ(i)]+=pp1[EE(i)]/Etot*pp1[FZ(j)];
	}
    } 
  */

  if(verbose) 
    {
      printf("=== pp0 on ===\n");
      print_Nvector(pp0,NV);      
    }

  if(verbose)
    {
      Etot=pp0[EE(0)]+pp0[EE(1)]+pp0[EE(2)]+pp0[EE(3)];
      FXtot=pp0[FX(0)]+pp0[FX(1)]+pp0[FX(2)]+pp0[FX(3)];
      FYtot=pp0[FY(0)]+pp0[FY(1)]+pp0[FY(2)]+pp0[FY(3)];
      FZtot=pp0[FZ(0)]+pp0[FZ(1)]+pp0[FZ(2)]+pp0[FZ(3)];

      if(verbose) printf("Tots: %e %e %e %e\n",Etot,FXtot,FYtot,FZtot);
 
    }



  //back to code coordinates
  prad_on2lab(pp0,pp,ggg);
  p2u(pp,uu,gg,GG);

  if(verbose) 
    {
      printf("=== pp ===\n");
      print_Nvector(pp,NV);
      getchar();
    }
  
  return 0;
}

//***********************************************************************************
//******* takes primitives and closes Rij in arbitrary frame for all fluids *********
//***********************************************************************************
int
calc_Rij_mf(ldouble *pp0, ldouble gg[][5], ldouble GG[][5], ldouble Rij[][4][4])
{
#ifdef MULTIRADFLUID
  ldouble pp[NV],Erf;
  int verbose=0;
  int i,j,irf;
 //relative velocity
  ldouble urfcon[4];
  //covariant formulation

#ifdef LABRADFLUXES
  my_err("LABRADFLUXES not implemented in calc_Rij_mf\n");
#else
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];
#endif


  for(irf=0;irf<NRF;irf++)
    {
      //radiative energy density in the radiation rest frame
      Erf=pp[EE(irf)];

      urfcon[0]=0.;
      urfcon[1]=pp[FX(irf)];
      urfcon[2]=pp[FY(irf)];
      urfcon[3]=pp[FZ(irf)];
      //converting to lab four-velocity
      conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
      //lab frame stress energy tensor:
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  Rij[irf][i][j]=4./3.*Erf*urfcon[i]*urfcon[j]+1./3.*Erf*GG[i][j];

    }

  

  return 0;
#endif
}

/************************************************************************/
/******* calculates wavespeeds in the lab frame takin 1/@3 in ************/
/******* radiative rest frame and boosting it to lab frame ***************/
/******* using the HARM algorithm for all fluids ***********************/
/************************************************************************/
/************************************************************************/
/******* returns one wavespeed for all fluids ****************************/
/************************************************************************/
int
calc_rad_wavespeeds_mf_total(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble tautot[3],ldouble *aval)
{
#ifdef MULTIRADFLUID
  int i,j,irf;
  
  //metric
  ldouble g00=gg[0][0];
  ldouble g03=gg[0][3];
  ldouble g30=g03;
  ldouble g11=gg[1][1];
  ldouble g22=gg[2][2];
  ldouble g33=gg[3][3];

  //inversed metric
  ldouble G00=GG[0][0];
  ldouble G03=GG[0][3];
  ldouble G11=GG[1][1];
  ldouble G22=GG[2][2];
  ldouble G33=GG[3][3];
  ldouble G30=G03;

#ifdef LABRADFLUXES
  my_err("LABRADFLUXES not implemented for MULTIRADFLUID\n");
#endif

  for(i=0;i<6;i++)
    aval[i]=0.;
  
  for(irf=0;irf<NRF;irf++)
    {
      //radiative energy density in the radiation rest frame
      ldouble Erf=pp[EE(irf)];
      //relative four-velocity
      ldouble urfcon[4];
      urfcon[0]=0.;
      urfcon[1]=pp[FX(irf)];
      urfcon[2]=pp[FY(irf)];
      urfcon[3]=pp[FZ(irf)];

      //converting to lab four-velocity
      conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

      //square of radiative wavespeed in radiative rest frame
      ldouble rv2rad = 1./3.;
      ldouble rv2,rv2tau;

      //**********************************************************************
      //algorithm from HARM to transform the fluid frame wavespeed into lab frame
      //**********************************************************************

      ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,wspeed2;
      ldouble axl,axr,ayl,ayr,azl,azr;
      axl=axr=ayl=ayr=azl=azr=1.;
   
      //**********************************************************************
      //**********************************************************************
      int dim;
      for(dim=0;dim<3;dim++)
	{
	  //characterisitic limiter based on the optical depth
	  //TODO: validate against opt.thick tests
	  if(tautot[dim]>0.) 
	    {
	      rv2tau=4./3./tautot[dim]*4./3./tautot[dim];
	      rv2=my_min(rv2rad,rv2tau);		     
	    }
	  else
	    rv2=rv2rad;
      
	  Acov[0]=0.;
	  Acov[1]=0.;
	  Acov[2]=0.;
	  Acov[3]=0.;
	  Acov[dim+1]=1.;
	  indices_12(Acov,Acon,GG);
  
	  Bcov[0]=1.;
	  Bcov[1]=0.;
	  Bcov[2]=0.;
	  Bcov[3]=0.;
	  indices_12(Bcov,Bcon,GG);

	  Asq = dot(Acon,Acov);
	  Bsq = dot(Bcon,Bcov);
	  Au = dot(Acov, urfcon);
	  Bu = dot(Bcov, urfcon);
	  AB = dot(Acon, Bcov);
	  Au2 = Au * Au;
	  Bu2 = Bu * Bu;
	  AuBu = Au * Bu;

	  wspeed2=rv2;
	  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
	  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
	  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
	  if(discr<0.) {printf("x1discr in ravespeeds lt 0\n"); discr=0.;}
	  discr = sqrt(discr);
	  ldouble cst1 = -(-B + discr) / (2. * A);
	  ldouble cst2 = -(-B - discr) / (2. * A);  

	  axl = my_min(cst1,cst2);
	  axr = my_max(cst1,cst2);

	  aval[dim*2+0]=my_min(axl,aval[dim*2+0]);
	  aval[dim*2+1]=my_max(axr,aval[dim*2+1]);
	}
    }
 

  return 0;
#endif
}

/************************************************************************/
/******* calculates wavespeeds in the lab frame takin 1/@3 in ************/
/******* radiative rest frame and boosting it to lab frame ***************/
/******* using the HARM algorithm for all fluids ***********************/
/******* no tau limiting **************************************************/
/************************************************************************/
/******* returns one wavespeed for each fluids ****************************/
/************************************************************************/
int
calc_rad_wavespeeds_pure_mf_each(ldouble *pp,void *ggg,ldouble aval[][6])
{
#ifdef MULTIRADFLUID

  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  int i,j,irf;
  
  //metric
  ldouble g00=gg[0][0];
  ldouble g03=gg[0][3];
  ldouble g30=g03;
  ldouble g11=gg[1][1];
  ldouble g22=gg[2][2];
  ldouble g33=gg[3][3];

  //inversed metric
  ldouble G00=GG[0][0];
  ldouble G03=GG[0][3];
  ldouble G11=GG[1][1];
  ldouble G22=GG[2][2];
  ldouble G33=GG[3][3];
  ldouble G30=G03;

#ifdef LABRADFLUXES
  my_err("LABRADFLUXES not implemented for MULTIRADFLUID\n");
#endif
  
  for(irf=0;irf<NRF;irf++)
    {      
      for(i=0;i<6;i++)
	aval[irf][i]=0.;

      //radiative energy density in the radiation rest frame
      ldouble Erf=pp[EE(irf)];
      //relative four-velocity
      ldouble urfcon[4];
      urfcon[0]=0.;
      urfcon[1]=pp[FX(irf)];
      urfcon[2]=pp[FY(irf)];
      urfcon[3]=pp[FZ(irf)];

      //converting to lab four-velocity
      conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

      //square of radiative wavespeed in radiative rest frame
      ldouble rv2rad = 1./3.;
      ldouble rv2,rv2tau;

      //**********************************************************************
      //algorithm from HARM to transform the fluid frame wavespeed into lab frame
      //**********************************************************************

      ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,wspeed2;
      ldouble axl,axr,ayl,ayr,azl,azr,cst1,cst2;
      axl=axr=ayl=ayr=azl=azr=1.;
   
      //**********************************************************************
      //**********************************************************************
      int dim;
      for(dim=0;dim<3;dim++)
	{
	  rv2=rv2rad;
      
	  Acov[0]=0.;
	  Acov[1]=0.;
	  Acov[2]=0.;
	  Acov[3]=0.;
	  Acov[dim+1]=1.;
	  indices_12(Acov,Acon,GG);
  
	  Bcov[0]=1.;
	  Bcov[1]=0.;
	  Bcov[2]=0.;
	  Bcov[3]=0.;
	  indices_12(Bcov,Bcon,GG);

	  Asq = dot(Acon,Acov);
	  Bsq = dot(Bcon,Bcov);
	  Au = dot(Acov, urfcon);
	  Bu = dot(Bcov, urfcon);
	  AB = dot(Acon, Bcov);
	  Au2 = Au * Au;
	  Bu2 = Bu * Bu;
	  AuBu = Au * Bu;

	  wspeed2=rv2;
	  B = 2. * (AuBu * (1.0-wspeed2)  - AB*wspeed2);
	  A = Bu2 * (1.0 - wspeed2) - Bsq * wspeed2;
	  discr = 4.0 * wspeed2 * ((AB * AB - Asq * Bsq) * wspeed2 + (2.0 * AB * Au * Bu - Asq * Bu2 - Bsq * Au2) * (wspeed2 - 1.0));
	  if(discr<0.) {printf("x1discr in ravespeeds lt 0\n"); discr=0.;}
	  discr = sqrt(discr);
	  cst1 = -(-B + discr) / (2. * A);
	  cst2 = -(-B - discr) / (2. * A);  

	  axl = my_min(cst1,cst2);
	  axr = my_max(cst1,cst2);

	  aval[irf][dim*2+0]=my_min(axl,aval[irf][dim*2+0]);
	  aval[irf][dim*2+1]=my_max(axr,aval[irf][dim*2+1]);
	}

      if(geom->ix==100 && (pp[7]!=0. || pp[FX(1)]!=0.) && 0)
	{
	  printf("--- %d\n",irf);
	  print_4vector(&pp[EE(irf)]);
	  print_4vector(urfcon);
	  printf("%e %e\n",cst1,cst2);
	  print_Nvector(aval[irf],6);
	  getchar();
	}
    }

 
 

  return 0;
#endif
}

//**********************************************************************
//******* takes E and F^i from primitives (artificial) **********************
//******* and calculates radiation stress ******************************
//******* tensor R^ij in fluid frame using M1 closure scheme ***********
//**********************************************************************
int
calc_Rij_ff_mf(ldouble *pp, ldouble Rij[][4][4])
{
#ifdef MULTIRADFLUID
  int irf;

  for(irf=0;irf<NRF;irf++)
    {

      ldouble E=pp[EE(irf)];
      ldouble F[3]={pp[FX(irf)],pp[FY(irf)],pp[FZ(irf)]};

      ldouble nx,ny,nz,nlen,f;

      nx=F[0]/E;
      ny=F[1]/E;
      nz=F[2]/E;

      nlen=sqrt(nx*nx+ny*ny+nz*nz);
 
#ifdef EDDINGTON_APR
      f=1./3.;
#else  
      if(nlen>=1.)
	{
	  f=1.;
	}
      else //M1
	f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
#endif
  
      if(nlen>0) 
	{
	  nx/=nlen;
	  ny/=nlen;
	  nz/=nlen;
	}
      else
	{
	  ;
	}
 
      Rij[irf][0][0]=E;
      Rij[irf][0][1]=Rij[irf][1][0]=F[0];
      Rij[irf][0][2]=Rij[irf][2][0]=F[1];
      Rij[irf][0][3]=Rij[irf][3][0]=F[2];

      Rij[irf][1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
      Rij[irf][1][2]=E*(.5*(3.*f - 1.)*nx*ny);
      Rij[irf][1][3]=E*(.5*(3.*f - 1.)*nx*nz);

      Rij[irf][2][1]=E*(.5*(3.*f - 1.)*ny*nx);
      Rij[irf][2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
      Rij[irf][2][3]=E*(.5*(3.*f - 1.)*ny*nz);

      Rij[irf][3][1]=E*(.5*(3.*f - 1.)*nz*nx);
      Rij[irf][3][2]=E*(.5*(3.*f - 1.)*nz*ny);
      Rij[irf][3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

    }

  return 0;
#endif
}


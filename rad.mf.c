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

  //if(geom->iy==0 &&  fabs(Fon[0]/Eon)>1.e-2) verbose=1;
  //if(geom->iy==2 && radius>4. && radius<4.5 && fabs(Fon[0]/Eon)>1.e-2) verbose=1;

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
      ldouble EF[4]={ppon[EE(irf)],ppon[FX(irf)],ppon[FY(irf)],ppon[FZ(irf)]};
      ldouble f0=sqrt(EF[1]*EF[1]+EF[2]*EF[2]+EF[3]*EF[3])/EF[0];

      if(verbose)
	{
	  printf("===\nfluid:  %d\n",irf);
	}

      if(f0<1.e-8 || EF[0]<1.e-8*Eon)
	//flux small or empty wedge
	{
	  ppon2[EE(irf)]+=EF[0];
	  ppon2[FX(irf)]+=EF[1];
	  ppon2[FY(irf)]+=EF[2];
	  ppon2[FZ(irf)]+=EF[3];
	}
      else
	{
	  //components along the tetrad
	  ldouble Frp = EF[1];
	  ldouble Ftp = EF[2];
	  ldouble Ffp = EF[3];
	  ldouble sumF = fabs(Frp)+fabs(Ftp)+fabs(Ffp);

	  //fraction applied:
	  ldouble frac;
	  ldouble fraccorr[3]={1.,1.,1.};
	  //TODO: better estimate the velocity?
	  if(dt<0.)
	    frac=1.;
	  else
	    {
	      frac = dt / (radius / (1./3.)) * 1./5.;
	      if(frac>1.) frac=1.;
	    }
	  
	  if(verbose) printf("frac applied: %e\n",frac);

	  ppon2[EE(irf)]+=(1.-frac)*EF[0];
	  ppon2[FX(irf)]+=(1.-frac)*EF[1];
	  ppon2[FY(irf)]+=(1.-frac)*EF[2];
	  ppon2[FZ(irf)]+=(1.-frac)*EF[3];

	  //new components
	  for(jj=0;jj<3;jj++)
	    {
	      ldouble En,Fn[3],avals[6],fracloc;

	      if(jj==0) //along r
		{
		  En=fabs(Frp)/sumF*EF[0];
		  Fn[0]=Frp;
		  Fn[1]=0.;
		  Fn[2]=0.;
		  if(Frp==0.) continue;

		  if(Frp>0.) ii=1;
		  if(Frp<0.) ii=0;
		}

	      if(jj==1) //along theta/z
		{
		  En=fabs(Ftp)/sumF*EF[0];
		  Fn[0]=0.;
		  Fn[1]=Ftp;
		  Fn[2]=0.;
		  if(Ftp==0.) continue;
		  
		  if(Ftp>0.) ii=3;
		  if(Ftp<0.) ii=2;
		}

	      if(jj==2) //along phi
		{
		  En=fabs(Ffp)/sumF*EF[0];
		  Fn[0]=0.;
		  Fn[1]=0.;
		  Fn[2]=Ffp;
		  if(Ffp==0.) continue;

		  if(Ffp>0.) ii=5;
		  if(Ffp<0.) ii=4;
		}


	      //to check if new fluid goes out of bounds
	      ldouble Fx0=ppon2[FX(ii)];
	      ldouble Fy0=ppon2[FY(ii)];
	      ldouble Fz0=ppon2[FZ(ii)];
	      ldouble E0=ppon2[EE(ii)];
	      ldouble Fx1=Fn[0];
	      ldouble Fy1=Fn[1];
	      ldouble Fz1=Fn[2];
	      ldouble E1=En;

	      fracloc=(-2*E0*E1 + 2*Fx0*Fx1 + 2*Fy0*Fy1 + 2*Fz0*Fz1 + 
		       Sqrt(Power(2*E0*E1 - 2*Fx0*Fx1 - 2*Fy0*Fy1 - 2*Fz0*Fz1,2) - 
			    4*(Power(E0,2) - Power(Fx0,2) - Power(Fy0,2) - Power(Fz0,2))*
			    (Power(E1,2) - Power(Fx1,2) - Power(Fy1,2) - Power(Fz1,2))))/
		(2.*(Power(E1,2) - Power(Fx1,2) - Power(Fy1,2) - Power(Fz1,2)));

	      if(fracloc<0. || isnan(fracloc) || fracloc<frac)
		fracloc=frac;

	      ppon2[EE(ii)]+=fracloc*En;
	      ppon2[FX(ii)]+=fracloc*Fn[0];
	      ppon2[FY(ii)]+=fracloc*Fn[1];
	      ppon2[FZ(ii)]+=fracloc*Fn[2];

	      /*
	      ppon2[EE(irf)]+=(1.-(fracloc-frac))*EF[0];
	      ppon2[FX(irf)]+=(1.-(fracloc-frac))*EF[1];
	      ppon2[FY(irf)]+=(1.-(fracloc-frac))*EF[2];
	      ppon2[FZ(irf)]+=(1.-(fracloc-frac))*EF[3];
	      */
	      /*
	      if(sqrt(ppon2[FX(ii)]*ppon2[FX(ii)]+ppon2[FX(ii)]*ppon2[FX(ii)]+ppon2[FX(ii)]*ppon2[FX(ii)])
		 /ppon2[EE(ii)]>1.0001) 
		printf("F/cE exceeded!\n");
	      */
	   
	      if(verbose)
		{
		  printf("\n component %d:\n %e %e %e %e\n",jj,En,Fn[0]/En,Fn[1]/En,Fn[2]/En);
		  printf("fracs:\n %e %e \n",frac,fracloc);
		  printf("new %d fluid:\n %e %e %e %e \n\n",ii,ppon2[EE(ii)],ppon2[FX(ii)]/ppon2[EE(ii)],ppon2[FY(ii)]/ppon2[EE(ii)],ppon2[FZ(ii)]/ppon2[EE(ii)]);
		}

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
//******* corrects the distribution in azimuth by splitting highly azimuthal *******
//******* beams into almost-radial and more azimuthal one to avoid the *************
//******* the inner funnel *********************************************************
//***********************************************************************************
int
mf_correct_in_azimuth_old(ldouble *pp, ldouble *uu, void* ggg, ldouble dt)
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

  //if(geom->iy==0 &&  fabs(Fon[0]/Eon)>1.e-2) verbose=1;
  //if(geom->iy==2 && radius>4. && radius<4.5 && fabs(Fon[0]/Eon)>1.e-2) verbose=1;

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
      //ksi=(M_PI/2./9.)/5.;
      //ksi=0.;

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

      if(phi<=1.01*M_PI/2. || phi>=.99*3.*M_PI/2. || f0<1.e-8 || fabs(phi-M_PI)<ksi)
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
	  //and z-component
	  ldouble Ftp = EF[2];

	  //length of the original
	  ldouble Fzero = sqrt(EF[1]*EF[1]+EF[3]*EF[3]);
	  
	  //fraction applied:
	  ldouble frac;
	  //TODO: better estimate the velocity?
	  if(dt<0.)
	    frac=1.;
	  else
	    {
	      frac = dt / (radius / (1./3.)) * 3.;
	      if(frac>1.) frac=1.;
	    }

	  
	  if(verbose) printf("frac applied: %e\n",frac);

	  ppon2[EE(irf)]+=(1.-frac)*EF[0];
	  ppon2[FX(irf)]+=(1.-frac)*EF[1];
	  ppon2[FY(irf)]+=(1.-frac)*EF[2];
	  ppon2[FZ(irf)]+=(1.-frac)*EF[3];

	  //new components
	  for(jj=0;jj<3;jj++)
	    {
	      ldouble En,Fn[3],avals[6];

	      if(jj==0) //roughly along r
		{
		  En=   Frp/(Frp+Ffp+Ftp)*EF[0];
		  Fn[0]=Frp*erp[0];
		  Fn[1]=0.;
		  Fn[2]=Frp*erp[2];
		}

	      if(jj==1) //roughly along phi
		{
		  //the other component
		  En=   Ffp/(Frp+Ffp+Ftp)*EF[0];
		  Fn[0]=Ffp*efp[0];
		  Fn[1]=0.;
		  Fn[2]=Ffp*efp[2];
		}

	      if(jj==2) //z component
		{
		  //the other component
		  En=   Ftp/(Frp+Ffp+Ftp)*EF[0];
		  Fn[0]=0.;
		  Fn[1]=Ftp;
		  Fn[2]=0.;
		}

	      //distributing it over wedges
	      ldouble Avec[NRF];
	      ldouble SKEW,MINVEL;
	      if(jj<2) {SKEW=50.;MINVEL=1.e-3;}
	      else {SKEW=10.;MINVEL=1.e-2;}
	      calc_rad_wavespeeds_on(Fn[0]/En,Fn[1]/En,Fn[2]/En,avals);
	      redistribute_with_velocities(avals,Avec,SKEW,MINVEL);

	      if(verbose)
		{
		  printf("component %d:\n %e %e %e %e\n",jj,En,Fn[0]/En,Fn[1]/En,Fn[2]/En);
		  printf("wavespeeds %d:\n %e %e | %e %e | %e %e\n",jj,avals[0],avals[1],avals[2],avals[3],avals[4],avals[5]);
		  printf("coefficients %d:\n %e %e %e %e %e %e\n",jj,Avec[0],Avec[1],Avec[2],Avec[3],Avec[4],Avec[5]);
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
//******* basing on skewed characteristic velocitites **********************************
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
redistribute_with_velocities(ldouble avals[6],ldouble A[NRF],ldouble skew,ldouble MINVEL)
{
  ldouble expskew(ldouble vel, ldouble skew) 
  {
    return (exp(skew*vel)-1.0)/(exp(skew)-1.0);
  }

  ldouble vxl,vxr,vyl,vyr,vzl,vzr,sumvel;

  vxl=fabs(my_min(avals[0],-MINVEL));
  vxr=fabs(my_max(avals[1],MINVEL));
  vyl=fabs(my_min(avals[2],-MINVEL));
  vyr=fabs(my_max(avals[3],MINVEL));
  if(NRF==6)
    {
      vzl=fabs(my_min(avals[4],-MINVEL));
      vzr=fabs(my_max(avals[5],MINVEL));
    }
 
  vxl=expskew(vxl,skew);
  vxr=expskew(vxr,skew);
  vyl=expskew(vyl,skew);
  vyr=expskew(vyr,skew);
  vzl=expskew(vzl,skew);
  vzr=expskew(vzr,skew);

  //wedges (x,y,z)
  //0 (x-)
  //1 (x+)
  //2 (y-)
  //3 (y+)
  //4 (z-)
  //5 (z+)

  if(NRF==6)
    {
      sumvel=vxl+vxr+vyl+vyr+vzl+vzr;
      A[0]=vxl/sumvel;
      A[1]=vxr/sumvel;
      A[2]=vyl/sumvel;
      A[3]=vyr/sumvel;
      A[4]=vzl/sumvel;
      A[5]=vzr/sumvel;
    }

  if(NRF==4)
    {
      if(NY==1)
	sumvel=vxl+vxr+vzl+vzr;
      if(NZ==1)
	sumvel=vxl+vxr+vyl+vyr;
      if(NY==1 && NZ==1)
	sumvel=vxl+vxr+2.*expskew(MINVEL,skew);

      //arbitrary power
      A[0]=vxl/sumvel;
      A[1]=vxr/sumvel;

      if(NZ==1)
	{
	  A[2]=vyl/sumvel;
	  A[3]=vyr/sumvel;
	}

      if(NY==1)
	{
	  A[2]=vzl/sumvel;
	  A[3]=vzr/sumvel;
	}

      if(NY==1 && NZ==1)
	{
	  A[2]=expskew(MINVEL,skew)/sumvel;
	  A[3]=expskew(MINVEL,skew)/sumvel;
	}
    }
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

  int verbose=0;
  
  struct geometry *geom
   = (struct geometry *) ggg;

  //  if(geom->ix==49   ) verbose=1;

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

 
  for(ii=0;ii<NRF;ii++)
    for(jj=0;jj<NRF;jj++)
      A[ii][jj]=0.;

  for(irf=0;irf<NRF;irf++)
    {
      //in aval[NRF][6] one has rad.char.wavespeeds for [irf] fluid
      //aval[irf][0] - min left going in x
      //aval[irf][1] - max right going in x
      //aval[irf][2] - min left going in y 
      //etc...

      ldouble Avec[NRF];
      ldouble SKEW=50.;
      ldouble MINVEL=1.e-4;
      redistribute_with_velocities(&aval[irf][0],Avec,SKEW,MINVEL);

      for(ii=0;ii<NRF;ii++)
	A[irf][ii]=Avec[ii];
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


//returns wavespeeds basing on ortonormal coefficients nx=Fx/E
//returns avals={vxr,vxl,vyr,vyl,vzr,vzl}
int calc_rad_wavespeeds_on(ldouble nx,ldouble ny,ldouble nz, ldouble *avals)
{
  int ii;
  ldouble j21,j22,j23,j24,j31,j32,j33,j34,j41,j42,j43,j44;
  ldouble a,b,c,d,e;
  ldouble cccbrtb,cccrt,x1,x2,x3,x4;
  gsl_complex z1,z2,z3,z4;
  ldouble e1,e2,e3,e4;

  ldouble nlen=sqrt(nx*nx+ny*ny+nz*nz);
  if(nlen>1.)
    {
      nx/=nlen;
      ny/=nlen;
      nz/=nlen;
    }
  
  if(isnan(nx)) nx=0.;
  if(isnan(ny)) ny=0.;
  if(isnan(nz)) nz=0.;

  if(fabs(nx)<1.e-20 && fabs(ny)<1.e-20 && fabs(nz)<1.e-20)
    {
      avals[0]=-sqrt(1./3.);
      avals[1]=sqrt(1./3.);
      avals[2]=-sqrt(1./3.);
      avals[3]=sqrt(1./3.);
      avals[4]=-sqrt(1./3.);
      avals[5]=sqrt(1./3.);
      return 0;
    }
  
  int idim;
  for(idim=0;idim<3;idim++)
    {
      if(idim==0)
	{
	  j21=(0. + 24.*Power(nx,2) - 68.*Power(nx,4) + 28.*Power(ny,2) - 64.*Power(nx,2)*Power(ny,2) + 4.*Power(ny,4) + 28.*Power(nz,2) - 64.*Power(nx,2)*Power(nz,2) + 8.*Power(ny,2)*Power(nz,2) + 4.*Power(nz,4) + 15.*Power(nx,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 20.*Power(nx,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(ny,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(ny,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(ny,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 8.*Power(ny,2)*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nz,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j22=nx*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nx,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2) + (1.*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/(Power(nx,2) + Power(ny,2) + Power(nz,2)));
	  j23=ny*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nx,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
	  j24=nz*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nx,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
	  j31=(0. - 4.*nx*ny - 72.*Power(nx,3)*ny - 72.*nx*Power(ny,3) - 72.*nx*ny*Power(nz,2) + 2.*nx*ny*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*ny*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,3)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*ny*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j32=(ny*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j33=(nx*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j34=(nx*ny*((0.5*nz*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nz*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j41=(0. - 4.*nx*nz - 72.*Power(nx,3)*nz - 72.*nx*Power(ny,2)*nz - 72.*nx*Power(nz,3) + 2.*nx*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,2)*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(nz,3)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j42=(nz*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j43=(nx*nz*((0.5*ny*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*ny*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j44=(nx*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);

	  a=1;
	  b=-j22 - j33 - j44;
	  c=-j21 - j23*j32 + j22*j33 - j24*j42 - j34*j43 + j22*j44 + j33*j44;
	  d=-(j24*j41) + j24*j33*j42 - j24*j32*j43 + j22*j34*j43 - j22*j33*j44 + j21*(j33 + j44) - j23*(j31 + j34*j42 - j32*j44);
	  e=j24*j33*j41 - j23*j34*j41 - j24*j31*j43 + j21*j34*j43 + j23*j31*j44 - j21*j33*j44;
	}   

      if(idim==1)
	{
	  j21=(0. - 4.*nx*ny - 72.*Power(nx,3)*ny - 72.*nx*Power(ny,3) - 72.*nx*ny*Power(nz,2) + 2.*nx*ny*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*ny*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,3)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*ny*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j22=(ny*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j23=(nx*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j24=(nx*ny*((0.5*nz*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nz*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j31=(0. + 28.*Power(nx,2) + 4.*Power(nx,4) + 24.*Power(ny,2) - 64.*Power(nx,2)*Power(ny,2) - 68.*Power(ny,4) + 28.*Power(nz,2) + 8.*Power(nx,2)*Power(nz,2) - 64.*Power(ny,2)*Power(nz,2) + 4.*Power(nz,4) + 13.*Power(nx,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nx,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 15.*Power(ny,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(ny,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 20.*Power(ny,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 8.*Power(nx,2)*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(ny,2)*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nz,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j32=nx*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(ny,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
	  j33=ny*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(ny,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2) + (1.*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/(Power(nx,2) + Power(ny,2) + Power(nz,2)));
	  j34=nz*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(ny,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
	  j41=(0. - 4.*ny*nz - 72.*Power(nx,2)*ny*nz - 72.*Power(ny,3)*nz - 72.*ny*Power(nz,3) + 2.*ny*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,2)*ny*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(ny,3)*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*ny*Power(nz,3)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j42=(ny*nz*((0.5*nx*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nx*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j43=(nz*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j44=(ny*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);

	  a=1;
	  b=-j22 - j33 - j44;
	  c=-j31 - j23*j32 + j22*j33 - j24*j42 - j34*j43 + j22*j44 + j33*j44;
	  d=-(j21*j32) - j34*j41 + j24*j33*j42 - j23*j34*j42 - j24*j32*j43 + j31*j44 + j23*j32*j44 + j22*(j31 + j34*j43 - j33*j44);
	  e=-(j24*j32*j41) + j22*j34*j41 + j24*j31*j42 - j21*j34*j42 - j22*j31*j44 + j21*j32*j44;

	}

      if(idim==2)
	{
	  j21=(0. - 4.*nx*nz - 72.*Power(nx,3)*nz - 72.*nx*Power(ny,2)*nz - 72.*nx*Power(nz,3) + 2.*nx*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,3)*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(ny,2)*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*nx*Power(nz,3)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j22=(nz*((0.5*Power(nx,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nx,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j23=(nx*nz*((0.5*ny*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*ny*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j24=(nx*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j31=(0. - 4.*ny*nz - 72.*Power(nx,2)*ny*nz - 72.*Power(ny,3)*nz - 72.*ny*Power(nz,3) + 2.*ny*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(nx,2)*ny*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*Power(ny,3)*nz*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 24.*ny*Power(nz,3)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j32=(ny*nz*((0.5*nx*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*nx*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j33=(nz*((0.5*Power(ny,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(ny,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j34=(ny*((0.5*Power(nz,2)*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - 1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))) + 0.5*(Power(nx,2) + Power(ny,2) + Power(nz,2))*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2);
	  j41=(0. + 28.*Power(nx,2) + 4.*Power(nx,4) + 28.*Power(ny,2) + 8.*Power(nx,2)*Power(ny,2) + 4.*Power(ny,4) + 24.*Power(nz,2) - 64.*Power(nx,2)*Power(nz,2) - 64.*Power(ny,2)*Power(nz,2) - 68.*Power(nz,4) + 13.*Power(nx,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(nx,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 13.*Power(ny,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 8.*Power(nx,2)*Power(ny,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 4.*Power(ny,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) + 15.*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(nx,2)*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 16.*Power(ny,2)*Power(nz,2)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))) - 20.*Power(nz,4)*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2));
	  j42=nx*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nz,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
	  j43=ny*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nz,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2));
	  j44=nz*((-41. + 12.*Power(nx,2) + 12.*Power(ny,2) + 12.*Power(nz,2) - 20.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) + (0.5*Power(nz,2)*(246. - 72.*Power(nx,2) - 72.*Power(ny,2) - 72.*Power(nz,2) + 120.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))))/((Power(nx,2) + Power(ny,2) + Power(nz,2))*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2)))*Power(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))),2)) - (1.*Power(nz,2)*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/Power(Power(nx,2) + Power(ny,2) + Power(nz,2),2) + (1.*(-1. + (3.*(3. + 4.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))/(5. + 2.*Sqrt(4. - 3.*(Power(nx,2) + Power(ny,2) + Power(nz,2))))))/(Power(nx,2) + Power(ny,2) + Power(nz,2)));

	  a=1;
	  b=-j22 - j33 - j44;
	  c=-(j23*j32) - j41 - j24*j42 - j34*j43 + j33*j44 + j22*(j33 + j44);
	  d=-(j21*j42) - j23*j34*j42 + j33*(j41 + j24*j42) - j31*j43 - j24*j32*j43 + j23*j32*j44 + j22*(j41 + j34*j43 - j33*j44);
	  e=j23*j32*j41 - j22*j33*j41 - j23*j31*j42 + j21*j33*j42 + j22*j31*j43 - j21*j32*j43;
	}

      //attemp to solve analytically:  
      int ret=gsl_poly_complex_solve_quartic(b/a,c/a,d/a,e/a,&z1,&z2,&z3,&z4);
  
      //if didn't work - solve numerically
      if(isnan(GSL_REAL(z1)) || isnan(GSL_REAL(z1)) || isnan(GSL_REAL(z1)) || isnan(GSL_REAL(z1)))
	{
            
	  double coef[5] = { (double)e,(double)d,(double)c,(double)b,(double)a };  
	  double z[8];      
	  gsl_poly_complex_workspace * w 
	    = gsl_poly_complex_workspace_alloc (5);      
	  int result=gsl_poly_complex_solve (coef, 5, w, z);      
	  gsl_poly_complex_workspace_free (w);
	  if(result!=GSL_SUCCESS)
	    {
	      printf("padaka w rozwiazywaniu czwormianu\n");
	      getchar();
	    }

	  //      return my_max(my_max(fabs(z[0]),fabs(z[2])),my_max(fabs(z[4]),fabs(z[6])));    
	  e1=z[0];
	  e2=z[2];
	  e3=z[4];
	  e4=z[6];
	}
      else
	{
	  e1=GSL_REAL(z1);
	  e2=GSL_REAL(z2);
	  e3=GSL_REAL(z3);
	  e4=GSL_REAL(z4);
	}

      ldouble al,ar;
      al=0.;
      if(e1<al) al=e1;
      if(e2<al) al=e2;
      if(e3<al) al=e3;
      if(e4<al) al=e4;
  
      ar=0.;
      if(e1>ar) ar=e1;
      if(e2>ar) ar=e2;
      if(e3>ar) ar=e3;
      if(e4>ar) ar=e4;

      avals[idim*2+1]=ar;
      avals[idim*2+0]=al;

    }

  return 0;
}

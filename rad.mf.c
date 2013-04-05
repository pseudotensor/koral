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

  redistribute_radfluids(pp,uu,&geom);
  //redistribute_radfluids_along_axes(pp,uu,&geom);

  u2p_rad(uu,pp,&geom,&iv);

  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz,uu[iv]);
      set_u(p,iv,ix,iy,iz,pp[iv]);
    }
	
  return 0;
}

//***********************************************************************************
//******* redistributes radiation fluids ***********************************************
//***********************************************************************************
int
redistribute_radfluids(ldouble *pp, ldouble *uu0, void* ggg)
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
//***********************************************************************************
int
redistribute_radfluids_new(ldouble *pp, ldouble *uu0, void* ggg)
{
#ifdef MULTIRADFLUID

  ldouble power=MFPOWER;
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

  ldouble MINVEL=1.e-3;

  for(irf=0;irf<NRF;irf++)
    {
      //in aval[NRF][6] one has rad.char.wavespeeds for [irf] fluid
      //aval[irf][0] - min left going in x
      //aval[irf][1] - max right going in x
      //aval[irf][2] - min left going in y 
      //etc...

      ldouble vxl,vxr,vyl,vyr,vzr,vzl,sumvel;
      vxl=my_min(aval[irf][0],-MINVEL);
      vxr=my_max(aval[irf][1],MINVEL);
      vyl=my_min(aval[irf][2],-MINVEL);
      vyr=my_max(aval[irf][3],MINVEL);
      vzl=my_min(aval[irf][4],-MINVEL);
      vzr=my_max(aval[irf][5],MINVEL);
    
      //wedges (x,y,z)
      //0 (x+)
      //1 (x-)
      //2 (y+)
      //3 (y-)
      //4 (z+)
      //5 (z-)

      vxl=pow(fabs(vxl),power);
      vxr=pow(fabs(vxr),power);
      vyl=pow(fabs(vyl),power);
      vyr=pow(fabs(vyr),power);
      vzl=pow(fabs(vzl),power);
      vzr=pow(fabs(vzr),power);

      sumvel=vxl+vxr+vyl+vyr+vzl+vzr;

      //test
      sumvel=vxl+vxr+vzl+vzr;

      //arbitrary power
      A[irf][0]=vxr/sumvel;
      A[irf][1]=vxl/sumvel;
      A[irf][2]=vyr/sumvel;
      A[irf][3]=vyl/sumvel;
      //A[irf][4]=vzr/sumvel;
      //A[irf][5]=vzl/sumvel;

      //test 2d, NRF =4, to work with along_axes
      A[irf][2]=vzr/sumvel;
      A[irf][3]=vzl/sumvel;

      /*
      //test - old approach
      ldouble vxr0=vxr;
      ldouble vzr0=vzr;
      ldouble vxl0=vxl;
      ldouble vzl0=vzl;
      
      vxr=sqrt(vxr0*vxr0+vzl0*vzl0);
      vzr=sqrt(vxr0*vxr0+vzr0*vzr0);
      vxl=sqrt(vxl0*vxl0+vzr0*vzr0);
      vzl=sqrt(vxl0*vxl0+vzl0*vzl0);
      A[irf][0]=vxr/(vxl+vxr)*vzr/(vzl+vzr);
      A[irf][1]=vxl/(vxl+vxr)*vzr/(vzl+vzr);
      A[irf][2]=vxl/(vxl+vxr)*vzl/(vzl+vzr);
      A[irf][3]=vxr/(vxl+vxr)*vzl/(vzl+vzr);
      */
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


//KORAL - rad.c
//radiation-related routines

#include "ko.h"

//*********************************************************************
//******* calculates total opacity over dx[] ***************************
//**********************************************************************
ldouble
calc_chi(ldouble *pp, ldouble *xx)
{
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  ldouble kappa=calc_kappa(rho,T,xx[1],xx[2],xx[3]);
  ldouble chi=kappa+calc_kappaes(rho,T,xx[1],xx[2],xx[3]);

  return chi;
}

//*********************************************************************
//******* calculates total opacity over dx[] ***************************
//**********************************************************************
int
calc_tautot(ldouble *pp, ldouble *xx, ldouble *dx, ldouble *tautot)
{
  ldouble rho=pp[RHO];
  ldouble u=pp[1];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  ldouble kappa=calc_kappa(rho,T,xx[1],xx[2],xx[3]);
  ldouble chi=kappa+calc_kappaes(rho,T,xx[1],xx[2],xx[3]);

  tautot[0]=chi*dx[0];
  tautot[1]=chi*dx[1];
  tautot[2]=chi*dx[2];

  return 0;
}

//**********************************************************************
//******* calculates abs opacity over dx[] ***************************
//**********************************************************************
int
calc_tauabs(ldouble *pp, ldouble *xx, ldouble *dx, ldouble *tauabs)
{
  ldouble rho=pp[RHO];
  ldouble u=pp[1];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;

  //xx[0] holds time
  ldouble kappa=calc_kappa(rho,T,xx[1],xx[2],xx[3]);

  tauabs[0]=kappa*dx[0];
  tauabs[1]=kappa*dx[1];
  tauabs[2]=kappa*dx[2];

  return 0;
}

      
//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame, returns ultimate deltas ***********************
//******* the fiducial approach *****************************************
//**********************************************************************

int f_implicit_lab_4dcon(ldouble *uu0,ldouble *uu,ldouble *pp0,ldouble dt,void* ggg,ldouble *f)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  ldouble pp[NV];
  int iv,i1,i2;
  for(iv=0;iv<NV;iv++)
    pp[iv]=pp0[iv];

  //opposite changes in gas quantities
  uu[1] = uu0[1] - (uu[6]-uu0[6]);
  uu[2] = uu0[2] - (uu[7]-uu0[7]);
  uu[3] = uu0[3] - (uu[8]-uu0[8]);
  uu[4] = uu0[4] - (uu[9]-uu0[9]);

  //calculating primitives  
  int corr[2],fixup[2],u2pret;
  //  print_Nvector(uu,NV);
  //  u2p_rad(uu,pp,ggg,&corr[1]);
  //  printf("%d \n",corr[1]);
  //print_Nvector(pp,NV);getchar();

  u2pret=u2p(uu,pp,ggg,corr,fixup);
  //printf("%d %d\n",corr[0],corr[1]);
  //print_Nvector(pp,NV);getchar();

  //if(corr[0]!=0 || corr[1]!=0) 
  //printf("corr: %d %d\n",corr[0],corr[1]);

  if(u2pret<-1) 
    {
      printf("implicit sub-sub-step failed\n");
      return -1; //allows for entropy but does not update conserved 
    }

  //radiative four-force
  ldouble Gi[4];
  calc_Gi(pp,ggg,Gi); 
  indices_21(Gi,Gi,gg);

  f[0] = uu[6] - uu0[6] + dt * gdetu * Gi[0];
  f[1] = uu[7] - uu0[7] + dt * gdetu * Gi[1];
  f[2] = uu[8] - uu0[8] + dt * gdetu * Gi[2];
  f[3] = uu[9] - uu0[9] + dt * gdetu * Gi[3];

  //fluid frame version for testing

  //zero - state (precalculate!)
  ldouble ucon0[4]={0.,pp0[VX],pp0[VY],pp0[VZ]},ucov0[4];
  conv_vels(ucon0,ucon0,VELPRIM,VEL4,gg,GG);
  indices_21(ucon0,ucov0,gg);
  ldouble Rij0[4][4],Rtt0;
  calc_Rij(pp0,ggg,Rij0);
  indices_2221(Rij0,Rij0,gg);
  Rtt0=0.;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      Rtt0+=-Rij0[i1][i2]*ucon0[i2]*ucov0[i1];

  /*
  boost22_lab2ff(Rij0,Rij0,pp0,gg,GG);
  trans22_cc2on(Rij0,Rij0,geom->tup);
  Rtt0=Rij0[0][0];
  */

  //new state
  ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]},ucov[4];
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);
  ldouble Rij[4][4],Rtt;
  calc_Rij(pp,ggg,Rij);
  indices_2221(Rij,Rij,gg);
  Rtt=0.;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      Rtt+=-Rij[i1][i2]*ucon[i2]*ucov[i1];

  /*
  boost22_lab2ff(Rij,Rij,pp,gg,GG);
  trans22_cc2on(Rij,Rij,geom->tup);
  Rtt=Rij[0][0];
  */

  ldouble T=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Ehat = -Rtt;
  ldouble dtau=dt/ucon[0];
  ldouble kappaabs=calc_kappa(pp[RHO],T,geom->xx,geom->yy,geom->zz);

  //printf("%.20e %.20e %.20e \n",uu[6],Gi[0],uu[6] - uu0[6] + dt * gdetu * Gi[0]);
  //printf("%.20e %.20e %.20e %.20e %.20e %.20e \n",uu[6],Rij[0][0],pp[6],Rtt,Ehat-4.*Pi*B,Rtt - Rtt0 - kappaabs*(Ehat-4.*Pi*B)*dtau);
  f[0]=Rtt - Rtt0 - kappaabs*(Ehat-4.*Pi*B)*dtau;
  
  return 0;
} 

int
print_state_implicit_lab_4dcon (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .13e % .13e % .13e % .3e "
	  "f(x) = % .13e % .13e % .13e % .13e\n",
	  iter,
	  //	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
	  x[0],x[1],x[2],x[3],f[0],f[1],f[2],f[3]);
  return 0;
}

int
solve_implicit_lab_4dcon(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose)
{
  int i1,i2,i3,iv,i,j;
  ldouble J[4][4],iJ[4][4];
  ldouble pp[NV],pp0[NV],uu[NV],uu0[NV],uu00[NV],uup[NV]; 
  ldouble f1[4],f2[4],f3hd[4],f3rad[4],xxx[4];

  ldouble (*gg)[5],(*GG)[5];

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz); //primitives corresponding to zero-state  
      uu[iv]=get_u(u,iv,ix,iy,iz);  
      uu00[iv]=uu[iv]; //total zero state
      uu0[iv]=uu[iv]; //zero state for substepping
      pp0[iv]=pp[iv]; //primitives
   }
  
  int corr[2],fixup[2];
  u2p(uu0,pp0,&geom,corr,fixup);
  p2u(pp0,uu0,&geom);

  ldouble EPS = 1.e-8;
  ldouble CONV = 1.e-6; 
  ldouble DAMP = 0.5;

  ldouble frdt = 1.0;
  ldouble dttot = 0.;

  int iter=0;
  int failed=0;

  //loop in dt
  if(verbose) 
    {
      ldouble xx,yy,zz;
      xx=get_x(ix,0);
      yy=get_x(iy,1);
      zz=get_x(iz,2);
      printf("=== i: %d %d %d\n=== x: %e %e %e\n",ix,iy,iz,xx,yy,zz);
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);
      print_metric(gg);
    }

  do 
    {

      if(verbose) 
	{
	  printf("====\n===\n Trying imp lab with frdt | dttot : %f | %f\n",frdt,dttot);
	}

      failed=0;
      iter=0;

      do
	{	 
	  iter++;
      
	  for(i=0;i<NV;i++)
	    {
	      uup[i]=uu[i];
	    }	

	  if(verbose>0)
	    {
	      for(i=0;i<4;i++)
		{
		  xxx[i]=uup[i+6];
		}  
	      print_Nvector(uu0,NV);
	      print_Nvector(uu,NV);
	      print_Nvector(pp0,NV);

	      int ret=f_implicit_lab_4dcon(uu0,uu,pp0,frdt*(1.-dttot)*dt,&geom,f1);
	      print_state_implicit_lab_4dcon (iter-1,xxx,f1); 
	      printf("f_lab_4dcon ret: %d\n",ret);
	    }


	  //values at base state
	  if(f_implicit_lab_4dcon(uu0,uu,pp0,frdt*(1.-dttot)*dt,&geom,f1)<0) 
	    {
	      failed=1;
		  
	      break;
	    }
	  
	  //calculating approximate Jacobian
	  for(j=0;j<4;j++)
	    {
	      ldouble del;

	      del=EPS*uup[6]; 

	      uu[j+6]=uup[j+6]-del;
	      
	      int fret=f_implicit_lab_4dcon(uu0,uu,pp0,frdt*(1.-dttot)*dt,&geom,f2);  

	      if(verbose>0)
		{
		  for(i=0;i<4;i++)
		    {
		      xxx[i]=uu[i+6];
		    }
		  print_state_implicit_lab_4dcon (iter-1,xxx,f2); 
		  printf("sub (%d) f_lab_4dcon ret: %d\n",j,fret);
		}

	      if(fret<0) 
		{
		  failed=1;	     		 
		}
  
	      //Jacobian matrix component
	      for(i=0;i<4;i++)
    	        {
		  J[i][j]=(f2[i] - f1[i])/(uu[j+6]-uup[j+6]);
		}

	      uu[j+6]=uup[j+6];

	      if(failed!=0) break;
	    }

	  if(failed!=0) break;
	  
	  //if(verbose)	    print_tensor(J);
	  
	  //inversion
	  if(inverse_44matrix(J,iJ)<0)
	    {
	      failed=1;
	      if(verbose) 
		{
		  print_tensor(J);
		  printf("Jacobian inversion failed\n");
		}
	      break;
	    }

	  //if(verbose)	    print_tensor(iJ);

	  //updating x
	  for(i=0;i<4;i++)
	    {
	      xxx[i]=uup[i+6];
	    }

	  for(i=0;i<4;i++)
	    {
	      for(j=0;j<4;j++)
		{
		  xxx[i]-=iJ[i][j]*f1[j];
		}
	    }

	  if(verbose>0)    print_state_implicit_lab_4dcon (iter,xxx,f1); 

	  for(i=0;i<4;i++)
	    {
	      uu[i+6]=xxx[i];
	    }
  
	   //opposite changes in gas quantities
	  uu[1] = uu0[1] - (uu[6]-uu0[6]);
	  uu[2] = uu0[2] - (uu[7]-uu0[7]);
	  uu[3] = uu0[3] - (uu[8]-uu0[8]);
	  uu[4] = uu0[4] - (uu[9]-uu0[9]);

	  //test convergence
	  for(i=0;i<4;i++)
	    {
	      f3rad[i]=(uu[i+6]-uup[i+6]);
	      f3hd[i]=(uu[i]-uup[i]);
	      
	      f3rad[i]=fabs(f3rad[i]/my_max(fabs(uup[6]),fabs(uup[i])));
	      f3hd[i]=fabs(f3hd[i]/my_max(fabs(uup[1]),fabs(uup[0])));
	    }

	  if(f3rad[0]<CONV && f3rad[1]<CONV && f3rad[2]<CONV && f3rad[3]<CONV)
	    {
	      if(verbose) printf("success ===\n");
	      break;
	    }

	  if(iter>10)
	    {
	      //return -1;
	      if(verbose) 
		{
		  printf("iter exceeded in solve_implicit_lab_4dcon() for frdt=%f | %f\n",frdt,dttot);	  
		}
	      failed=1;
	      break;
	      
	    }
     
	}
      while(1);

      if(failed==0) 
	{
	  if(verbose)
	    {
	      printf("output:\n");
	      print_Nvector(uu0,NV);
	      print_Nvector(uu,NV);
	      //	      getchar();
	    }

	  //opposite changes in gas quantities
	  uu[1] = uu0[1] - (uu[6]-uu0[6]);
	  uu[2] = uu0[2] - (uu[7]-uu0[7]);
	  uu[3] = uu0[3] - (uu[8]-uu0[8]);
	  uu[4] = uu0[4] - (uu[9]-uu0[9]);

	  //saving basis for next iteration if necessary
	  for(iv=0;iv<NV;iv++)
	    {
	      uu0[iv]=uu[iv];
	    }

	  if(frdt>=1.) 
	    {
	      if(verbose) {
		printf("worked!frdt=%f | %f===\n",frdt,dttot);
		if(dttot>0.) getchar();}
	      break;
	    }
	  else
	    {
	      if(verbose) 
		{
		  printf("worked but not the end! frdt=%f | %f===\n",frdt,dttot);		  
		}
	      dttot+=frdt*(1.-dttot);
	      frdt*=2.; 
	      if(frdt>1.) frdt=1.;
	    }

	  //update primitives corresponding to uu0
	  u2p(uu0,pp0,&geom,corr,fixup);
	  p2u(pp0,uu0,&geom);
	  continue;
	}
      
      //didn't work - decreasing time step
      frdt *= DAMP;

      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=uu0[iv];
	}
      
      if(frdt<0.00001 || 1)  //avoid substepping?
	{
	  if(verbose) 
	    {
	      printf("time step too small - aborting implicit_lab_4dcon() ===\n");
	    }
	  return -1;
	}
    }
  while(1);
  
  //  if(verbose) getchar();

  deltas[0]=uu[6]-uu00[6];
  deltas[1]=uu[7]-uu00[7];
  deltas[2]=uu[8]-uu00[8];
  deltas[3]=uu[9]-uu00[9];
  
  return 0;
}

//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame  working on primitives    ***********************
//******* rad or hydro (whichprim) **************************************
//**********************************************************************

int f_implicit_lab_4dprim(ldouble *pp,ldouble *uu0,ldouble *pp0,ldouble dt,void* ggg,ldouble *f,int *params)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int whichprim=params[0];
  int whicheq=params[1];

  ldouble uu[NV];
  int corr[2]={0,0},fixup[2]={0,0},u2pret,i1,i2;

  //total inversion, but only whichprim part matters
  p2u(pp,uu,geom);
 
  //opposite changes in the other quantities and inversion
  if(whichprim==1)
    {
      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);  

      u2pret=u2p(uu,pp,geom,corr,fixup); //total inversion (I should separate hydro from rad)
    }
  if(whichprim==0)
    {
      uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
      uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
      uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
      uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);

      u2pret=u2p_rad(uu,pp,geom,corr);
    }     

  //print_Nvector(uu,NV);getchar();

  //  if(corr[0]!=0 || corr[1]!=0) 
  //printf("corr: %d %d\n",corr[0],corr[1]);

  if(u2pret<-1) 
    {
      printf("implicit sub-sub-step failed\n");
      return -1; //allows for entropy but does not update conserved 
    }

  //radiative four-force
  ldouble Gi[4];
  calc_Gi(pp,ggg,Gi); 
  indices_21(Gi,Gi,gg);

  if(whichprim==1) //rad-primitives
    {
      f[0] = uu[6] - uu0[6] + dt * gdetu * Gi[0];
      f[1] = uu[7] - uu0[7] + dt * gdetu * Gi[1];
      f[2] = uu[8] - uu0[8] + dt * gdetu * Gi[2];
      f[3] = uu[9] - uu0[9] + dt * gdetu * Gi[3];
    }
  if(whichprim==0) //hydro-primitives
    {
      f[0] = uu[1] - uu0[1] - dt * gdetu * Gi[0];
      f[1] = uu[2] - uu0[2] - dt * gdetu * Gi[1];
      f[2] = uu[3] - uu0[3] - dt * gdetu * Gi[2];
      f[3] = uu[4] - uu0[4] - dt * gdetu * Gi[3];
    }

  //fluid frame version below
  if(whicheq==1)
    {
      
      //zero - state 
      ldouble ucon0[4]={0.,pp0[VX],pp0[VY],pp0[VZ]},ucov0[4];
      conv_vels(ucon0,ucon0,VELPRIM,VEL4,gg,GG);
      indices_21(ucon0,ucov0,gg);
      ldouble Rij0[4][4],Rtt0;
      calc_Rij(pp0,ggg,Rij0);
      indices_2221(Rij0,Rij0,gg);
      Rtt0=0.;
      for(i1=0;i1<4;i1++)
	for(i2=0;i2<4;i2++)
	  Rtt0+=-Rij0[i1][i2]*ucon0[i2]*ucov0[i1];

      //new state
      ldouble ucon[4]={0.,pp[VX],pp[VY],pp[VZ]},ucov[4];
      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
      indices_21(ucon,ucov,gg);
      ldouble Rij[4][4],Rtt;
      calc_Rij(pp,ggg,Rij);
      indices_2221(Rij,Rij,gg);
      Rtt=0.;
      for(i1=0;i1<4;i1++)
	for(i2=0;i2<4;i2++)
	  Rtt+=-Rij[i1][i2]*ucon[i2]*ucov[i1];

      ldouble T=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
      ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
      ldouble Ehat = -Rtt;
      ldouble dtau=dt/ucon[0];
      ldouble kappaabs=calc_kappa(pp[RHO],T,geom->xx,geom->yy,geom->zz);

      //fluid frame energy equation:
      f[0]=Rtt - Rtt0 - kappaabs*(Ehat-4.*Pi*B)*dtau;
    }
  
  return 0;
} 

int
print_state_implicit_lab_4dprim (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .13e % .13e % .13e % .3e "
	  "f(x) = % .13e % .13e % .13e % .13e\n",
	  iter,
	  //	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
	  x[0],x[1],x[2],x[3],f[0],f[1],f[2],f[3]);
  return 0;
}

int
solve_implicit_lab_4dprim(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int *params,int verbose)
{
  int i1,i2,i3,iv,i,j;
  ldouble J[4][4],iJ[4][4];
  ldouble pp[NV],pp0[NV],pp00[NV],ppp[NV],uu[NV],uu0[NV],uu00[NV],uup[NV]; 
  ldouble f1[4],f2[4],f3[4],xxx[4];

  ldouble (*gg)[5],(*GG)[5];

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  int whichprim=params[0];
  int whicheq=params[1];
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  for(iv=0;iv<NV;iv++)
    {
      uu00[iv]=get_u(u,iv,ix,iy,iz); //total zero-state    
      pp00[iv]=get_u(p,iv,ix,iy,iz); //only initial guess for u2p
    }
  
  int corr[2],fixup[2];
  u2p(uu00,pp00,&geom,corr,fixup);
  p2u(pp00,uu00,&geom);

  for(iv=0;iv<NV;iv++)
    {
      uu0[iv]=uu00[iv]; //zero state for substepping
      pp0[iv]=pp00[iv]; 
      uu[iv]=uu0[iv]; 
      pp[iv]=pp0[iv];     
    }
 
  ldouble EPS = 1.e-6;
  ldouble CONV = 1.e-6; 
  ldouble DAMP = 0.5;
  int sh;
  if(whichprim==0) 
    sh=UU; //solving in hydro primitives
  else
    sh=EE0; //solving in rad primitives

  ldouble frdt = 1.0;

  int iter=0;
  int failed=0;

  if(verbose) 
    {
      ldouble xx,yy,zz;
      xx=get_x(ix,0);
      yy=get_x(iy,1);
      zz=get_x(iz,2);
      printf("=== i: %d %d %d\n=== x: %e %e %e\n",ix,iy,iz,xx,yy,zz);
      print_Nvector(pp,NV);
      print_Nvector(uu,NV);
      print_metric(gg);
    }

  if(verbose) 
    {
      printf("====\n===\n Trying imp lab 4d prim with dt : %f \n",dt);
    }

  failed=0;
  iter=0;

  do //main solver loop
    {	 
      iter++;
      
      for(i=0;i<NV;i++)
	{
	  ppp[i]=pp[i];
	}	

      if(verbose>0)
	{
	  for(i=0;i<4;i++)
	    {
	      xxx[i]=ppp[i+sh];
	    }  
	  //print_Nvector(uu0,NV);
	  //print_Nvector(pp0,NV);
	  //print_Nvector(pp,NV);

	  int ret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,&geom,f1,params);
	  print_state_implicit_lab_4dprim (iter-1,xxx,f1); 
	  printf("f_lab_4dprim ret: %d\n",ret);
	}


      //values at base state
      if(f_implicit_lab_4dprim(pp,uu0,pp0,dt,&geom,f1,params)<0) 
	{
	  failed=1;
		  
	  break;
	}

      //calculating approximate Jacobian
      for(j=0;j<4;j++)
	{
	  ldouble del;

	  if(j==0) 
	    del=EPS*ppp[sh]; //eps in energy density
	  else
	    del=EPS; //eps in velocities

	  pp[j+sh]=ppp[j+sh]-del;

	      
	  int fret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,&geom,f2,params);  

	  if(verbose>0 && 0)
	    {
	      for(i=0;i<4;i++)
		{
		  xxx[i]=pp[i+sh];
		}
	      print_state_implicit_lab_4dprim (iter-1,xxx,f2); 
	      printf("sub (%d) f_lab_4dprim ret: %d\n",j,fret);
	    }

	  if(fret<0) 
	    {
	      return -1;
	    }
  
	  //Jacobian matrix component
	  for(i=0;i<4;i++)
	    {
	      J[i][j]=(f2[i] - f1[i])/(pp[j+sh]-ppp[j+sh]);
	    }

	  pp[j+sh]=ppp[j+sh];
	}
	  
      //inversion
      if(inverse_44matrix(J,iJ)<0)
	{
	  failed=1;
	  if(verbose) 
	    {
	      print_tensor(J);
	      printf("Jacobian inversion failed\n");
	    }
	  break;
	}

      //updating x
      for(i=0;i<4;i++)
	{
	  xxx[i]=ppp[i+sh];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      xxx[i]-=iJ[i][j]*f1[j];
	    }
	}

      if(verbose>0)    print_state_implicit_lab_4dprim (iter,xxx,f1); 

      for(i=0;i<4;i++)
	{
	  pp[i+sh]=xxx[i];
	}

      //implement overshooting check here
      //should not matter that before convergence check
	  
      //total inversion, but only whichprim part matters
      p2u(pp,uu,&geom);
      //opposite changes in the other quantities and inversion
      if(whichprim==1)
	{
	  uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
	  uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
	  uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
	  uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);
	  u2p(uu,pp,&geom,corr,fixup); //total inversion (I should separate hydro from rad)
	}
      if(whichprim==0)
	{
	  uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
	  uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
	  uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
	  uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);
	  u2p_rad(uu,pp,&geom,corr);
	}     

      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+sh]-ppp[i+sh]);
	  if(i==0)
	    f3[i]=fabs(f3[i]/ppp[sh]);
	  else
	    f3[i]=fabs(f3[i]/my_max(EPS,fabs(pp[i+sh])));
	}
	  
      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	{
	  if(verbose) printf("success ===\n");
	  break;
	}
	  
      if(iter>40)
	{
	  if(verbose) 
	    {
	      printf("iter exceeded in solve_implicit_lab_4dprim() for frdt=%f \n",dt);	  
	    }
	  return -1;
	}
    }
  while(1); //main solver loop

  //returns corrections to radiative primitives
  deltas[0]=uu[EE0]-uu00[EE0];
  deltas[1]=uu[FX0]-uu00[FX0];
  deltas[2]=uu[FY0]-uu00[FY0];
  deltas[3]=uu[FZ0]-uu00[FZ0];
  
  return 0;
}

//**********************************************************************
//* wrapper ************************************************************
//**********************************************************************

int
solve_implicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose)
{
  int RADPRIM=1;
  int HDPRIM=0;
  int LABEQ=0;
  int FFEQ=1;

  int params[2] = {RADPRIM, LABEQ};

  return solve_implicit_lab_4dprim(ix,iy,iz,dt,deltas,params,verbose);
  
  //return solve_implicit_lab_4dcon(ix,iy,iz,dt,deltas,verbose);
}


//**********************************************************************
//******* solves implicitly gas - radiation interaction ****************
//******* in the fluid frame, returns vector of deltas *****************
//******* the explicit-implicit approximate approach *******************
//******* used as the fail-safe backup method **************************
//**********************************************************************
int
solve_implicit_ff(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
{
  int i1,i2,i3,iv;
  ldouble pp[NV];

  ldouble (*gg)[5],(*GG)[5];

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  //temporary using local arrays
  gg=geom.gg;
  GG=geom.GG;

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
    }

  //transforming radiative primitives to ortonormal fluid frame
  prad_lab2ff(pp,pp,&geom);
  
  //four-force in the fluid frame
  ldouble Gi[4];
  calc_Gi_ff(pp,Gi);
  
  //implicit flux:
  ldouble rho=pp[RHO];
  ldouble u=pp[1];  
  ldouble E=pp[6];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble xx=get_x(ix,0);
  ldouble yy=get_x(iy,1);
  ldouble zz=get_x(iz,2);
  ldouble kappa=calc_kappa(rho,T,xx,yy,zz);
  ldouble chi=kappa+calc_kappaes(rho,T,xx,yy,zz);  
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;

  ldouble Fold[3]={pp[7],pp[8],pp[9]};
  ldouble Fnew[3];
  Fnew[0]=Fold[0]/(1.+dt*chi); 
  Fnew[1]=Fold[1]/(1.+dt*chi);
  Fnew[2]=Fold[2]/(1.+dt*chi);
 
  deltas[1]=Fnew[0]-Fold[0];
  deltas[2]=Fnew[1]-Fold[1];
  deltas[3]=Fnew[2]-Fold[2];

  //solving in parallel for E and u
  if(calc_LTE_ff(rho,&u,&E,dt,0)<0) 
    return -1;

  deltas[0]=E-pp[6];

  return 0;

}

//**********************************************************************
//******* solves explicitly gas - radiation interaction  *******************
//******* in the lab frame, returns vector of deltas **********************
//**********************************************************************
int
solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom.gg;
  GG=geom.GG;
  gdet=geom.gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int i1,i2,i3,iv;
  ldouble pp[NV];
  ldouble eup[4][4],elo[4][4];
  pick_T(emuup,ix,iy,iz,eup);
  pick_T(emulo,ix,iy,iz,elo);

  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz);      
    }

  ldouble Gi[4];
  calc_Gi(pp,&geom,Gi);
  indices_21(Gi,Gi,geom.gg);
  
  deltas[0]=-Gi[0]*dt*gdetu;
  deltas[1]=-Gi[1]*dt*gdetu;
  deltas[2]=-Gi[2]*dt*gdetu;
  deltas[3]=-Gi[3]*dt*gdetu;

  return 0;

}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical LTE solver used only by solve_implicit_ff()
//TODO: to be replaced with something faster or abandoned
//**********************************************************************
//**********************************************************************
//**********************************************************************

struct calc_LTE_ff_parameters
{
  ldouble rho,E,u,kappa,dt;
  int verbose;
};

double
f_calc_LTE_ff (double u, void *params)
{
  struct calc_LTE_ff_parameters *p 
    = (struct calc_LTE_ff_parameters *) params;
  
  ldouble rho=p->rho;
  ldouble Eold=p->E;
  ldouble uold=p->u;
  ldouble kappa=p->kappa;
  ldouble dt=p->dt;
  int verbose=p->verbose;

  ldouble pr= (GAMMA-1.)*(ldouble)u;
  ldouble T = pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Enew=(Eold +4.*Pi*kappa*B*dt)/(1.+kappa*dt);

  return u - uold - dt*(-4.*Pi*kappa*B + kappa*Enew);
}
     
double
df_calc_LTE_ff (double u, void *params)
{  
  struct calc_LTE_ff_parameters *p 
    = (struct calc_LTE_ff_parameters *) params;

  ldouble rho=p->rho;
  ldouble Eold=p->E;
  ldouble uold=p->u;
  ldouble kappa=p->kappa;
  ldouble dt=p->dt;
  int verbose=p->verbose;

  ldouble pr= (GAMMA-1.)*(ldouble)u;
  ldouble T = pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Enew=(Eold +4.*Pi*kappa*B*dt)/(1.+kappa*dt);

  ldouble dBdu=4.*B/(ldouble)u;
  ldouble dEdu=4.*Pi*kappa*dt*dBdu/(1.+dt*kappa);

  ldouble dfdu=1.+dt*kappa*4.*Pi*dBdu-dt*kappa*dEdu;

  return (double)dfdu;
}
     
void
fdf_calc_LTE_ff (double u, void *params, 
	       double *y, double *dy)
{ 
   struct calc_LTE_ff_parameters *p 
    = (struct calc_LTE_ff_parameters *) params;

  ldouble rho=p->rho;
  ldouble Eold=p->E;
  ldouble uold=p->u;
  ldouble kappa=p->kappa;
  ldouble dt=p->dt;
  int verbose=p->verbose;

  ldouble pr= (GAMMA-1.)*(ldouble)u;
  ldouble T = pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Enew=(Eold +4.*Pi*kappa*B*dt)/(1.+kappa*dt);

  ldouble dBdu=4.*B/(ldouble)u;
  ldouble dEdu=4.*Pi*kappa*dt*dBdu/(1.+dt*kappa);

  ldouble dfdu=1.+dt*kappa*4.*Pi*dBdu-dt*kappa*dEdu;
  ldouble f=u - uold - dt*(-4.*Pi*kappa*B + kappa*Enew);
 
  *y = f;
  *dy = dfdu;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
calc_LTE_ff(ldouble rho,ldouble *uint, ldouble *E,ldouble dt, int verbose)
{
  struct calc_LTE_ff_parameters cltep;
  cltep.rho=rho;
  cltep.u=*uint;
  cltep.E=*E;
  cltep.verbose=verbose;
  
  if(cltep.E<EEFLOOR && 0)
    {
      printf("imposing EEFLOOR 0\n");
      cltep.E=EEFLOOR;
    }
  
  ldouble p=(GAMMA-1.)*cltep.u;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  cltep.kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);

  ldouble tlte = 1./cltep.kappa;

  cltep.dt=(double)dt;

  //solving for E
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fdfsolver_type *TT;
  gsl_root_fdfsolver *s;
  double x0, x = cltep.u;
  gsl_function_fdf FDF;
     
  FDF.f = &f_calc_LTE_ff;
  FDF.df = &df_calc_LTE_ff;
  FDF.fdf = &fdf_calc_LTE_ff;
  FDF.params = &cltep;
  
  TT = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc (TT);
  gsl_root_fdfsolver_set (s, &FDF, x);

  do
    {
      iter++;
      status = gsl_root_fdfsolver_iterate (s);
      x0 = x;
      x = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (x, x0, 0, U2PRADPREC);
    }
  while (status == GSL_CONTINUE && iter < max_iter);
     
  if(iter>=max_iter) 
    {
      if(verbose) printf("lte error in calc_LTE_ff: %e %e %e %e\n",rho,*uint,*E,dt);
      return -1;
    }

  gsl_root_fdfsolver_free (s);

  *uint=x;
  ldouble pp1 = (GAMMA-1.)*(x);
  ldouble Ttu = pp1*MU_GAS*M_PROTON/K_BOLTZ/cltep.rho;
  ldouble Bp1 = SIGMA_RAD*pow(Ttu,4.)/Pi;
  *E=(cltep.E+4.*Pi*cltep.kappa*Bp1*dt)/(1.+cltep.kappa*dt);
  
  if(*uint<0 && 0)
    {
      *uint=1.e-30;
      printf("imposing ufloor in lte\n");
      getch();
    }
  return 0;
}
//**********************************************************************
//end of the numerical LTE solver
//**********************************************************************


//**********************************************************************
//******* opacities ****************************************************
//**********************************************************************
//absorption
ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
{
#include PR_KAPPA //PROBLEMS/XXX/kappa.c
}

//scattering
ldouble calc_kappaes(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)
{  
#include PR_KAPPAES //PROBLEMS/XXX/kappaes.c
}

//**********************************************************************
//****** takes radiative stress tensor and gas primitives **************
//****** and calculates contravariant four-force ***********************
//**********************************************************************
int
calc_Gi(ldouble *pp, void *ggg, ldouble Gi[4])
{
  int i,j,k;
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

#ifndef EDDINGTON_APR //M1 here
  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);

  //the four-velocity of fluid in lab frame
  ldouble ucon[4],ucov[4],vpr[3];
  ucon[1]=pp[2];
  ucon[2]=pp[3];
  ucon[3]=pp[4];
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);

  //covariant four-velocity
  indices_21(ucon,ucov,gg);  

  //gas properties
  ldouble rho=pp[RHO];
  ldouble u=pp[1];
  ldouble p= (GAMMA-1.)*(ldouble)u;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(Tgas,4.)/Pi;
  ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
  ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);
  ldouble chi=kappa+kappaes;

  //contravariant four-force in the lab frame

  //R^ab u_a u_b
  ldouble Ruu=0.;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Ruu+=Rij[i][j]*ucov[i]*ucov[j];

  ldouble Ru;
  for(i=0;i<4;i++)
    {
      Ru=0.;
      for(j=0;j<4;j++)
	Ru+=Rij[i][j]*ucov[j];
      Gi[i]=-chi*Ru - (kappaes*Ruu + kappa*4.*Pi*B)*ucon[i];
    }

#else //Eddington apr. here

  ldouble rho=pp[RHO];
  ldouble u=pp[1];
  ldouble ucov[4],ucon[4]={0,pp[2],pp[3],pp[4]};
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);
  ldouble EE=pp[6];
  ldouble Fcon[4]={0.,pp[7],pp[8],pp[9]};
  Fcon[0]=-1./ucov[0]*(Fcon[1]*ucov[1]+Fcon[2]*ucov[2]+Fcon[3]*ucov[3]); //F^0 u_0 = - F^i u_i
 
  ldouble p= (GAMMA-1.)*u;
  ldouble T = p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
  ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);
  ldouble chi=kappa+kappaes;

  for(i=0;i<4;i++)
    Gi[i]=kappa*(EE-4.*Pi*B)*ucon[i] + chi * Fcon[i];
#endif  

  //print_4vector(Gi);getchar();

  return 0;
}

//**********************************************************************
//******* takes fluid frame E,F^i in place of radiative primitives   *******
//******* and calculates force G^\mu in the fluid frame ****************
//**********************************************************************
int
calc_Gi_ff(ldouble *pp, ldouble Gi[4])
{
  ldouble rho=pp[RHO];
  ldouble u=pp[1];
  ldouble E=pp[6];
  ldouble F[3]={pp[7],pp[8],pp[9]};

  ldouble p= (GAMMA-1.)*(ldouble)u;
  ldouble T = p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
  ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);
  ldouble chi=kappa+kappaes;

  Gi[0]=kappa*(E-4.*Pi*B);
  Gi[1]=chi*F[0];
  Gi[2]=chi*F[1];
  Gi[3]=chi*F[2];

  return 0;
}


//***********************************************************************************
//******* takes primitives and closes R^ij in arbitrary frame ****************************
//***********************************************************************************
int
calc_Rij(ldouble *pp0, void *ggg, ldouble Rij[][4])
{
#ifdef RADIATION
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  ldouble (*tup)[4],(*tlo)[4];
  tup=geom->tup;
  tlo=geom->tlo;

  ldouble pp[NV],Erf;
  int verbose=0;
  int i,j;
 //relative velocity
  ldouble urfcon[4];
  //covariant formulation

#ifdef LABRADFLUXES
  //artificially puts pp=uu and converts them to urf and Erf using the regular converter
  u2p_rad_urf(pp0,pp,ggg,&i);
#else
  for(i=0;i<NV;i++)
    pp[i]=pp0[i];
#endif //LABRADFLUXES

#ifndef EDDINGTON_APR //M1 here
  //radiative energy density in the radiation rest frame
  Erf=pp[6];

  urfcon[0]=0.;
  urfcon[1]=pp[7];
  urfcon[2]=pp[8];
  urfcon[3]=pp[9];
  //converting to lab four-velocity
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
  //lab frame stress energy tensor:
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]=4./3.*Erf*urfcon[i]*urfcon[j]+1./3.*Erf*GG[i][j];

#else //Eddington

  ldouble h[4][4];
  ldouble ucov[4],ucon[4]={0,pp[2],pp[3],pp[4]};
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);
  ldouble EE=pp[6];
  ldouble Fcon[4]={0.,pp[7],pp[8],pp[9]};
  Fcon[0]=-1./ucov[0]*(Fcon[1]*ucov[1]+Fcon[2]*ucov[2]+Fcon[3]*ucov[3]); //F^0 u_0 = - F^i u_i
  //projection tensor
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      h[i][j]=GG[i][j] + ucon[i]*ucon[j];
  //Fragile's formula
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]=EE*ucon[i]*ucon[j] + Fcon[i]*ucon[j] + Fcon[j]*ucon[i] + 1./3.*EE*delta(i,j)*h[i][j];
#endif


#if (RADVISCOSITY!=NOVISCOSITY)
  ldouble Tvisc[4][4];
  calc_visc_Rij(pp,ggg,Tvisc,Rij);
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Rij[i][j]+=Tvisc[i][j];
#endif  

#endif //RADIATION

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
int
calc_visc_Rij(ldouble *pp, void* ggg, ldouble Tvisc[][4], ldouble Rij[][4])
{
  int i,j;
  ldouble Rij0[4][4];
  struct geometry *geom
   = (struct geometry *) ggg;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Tvisc[i][j]=0.;
	Rij0[i][j]=Rij[i][j];
      }

#if (RADVISCOSITY==SHEARVISCOSITY)
  ldouble shear[4][4];
  ldouble nu,vdiff2;
  calc_rad_shearviscosity(pp,geom,shear,&nu,&vdiff2);
  ldouble Erf=pp[EE0];

  //multiply by viscosity to get viscous tensor
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Tvisc[i][j]= -2. * nu * Erf * shear[i][j];
      }
#endif 

  return 0;

}

//**********************************************************************
//******* takes E and F^i from primitives (artificial) **********************
//******* and calculates radiation stress ******************************
//******* tensor R^ij in fluid frame using M1 closure scheme ***********
//**********************************************************************
int
calc_Rij_ff(ldouble *pp, ldouble Rij[][4])
{
  int irf=0;
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
 
  Rij[0][0]=E;
  Rij[0][1]=Rij[1][0]=F[0];
  Rij[0][2]=Rij[2][0]=F[1];
  Rij[0][3]=Rij[3][0]=F[2];

  Rij[1][1]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nx*nx);
  Rij[1][2]=E*(.5*(3.*f - 1.)*nx*ny);
  Rij[1][3]=E*(.5*(3.*f - 1.)*nx*nz);

  Rij[2][1]=E*(.5*(3.*f - 1.)*ny*nx);
  Rij[2][2]=E*(.5*(1.-f) + .5*(3.*f - 1.)*ny*ny);
  Rij[2][3]=E*(.5*(3.*f - 1.)*ny*nz);

  Rij[3][1]=E*(.5*(3.*f - 1.)*nz*nx);
  Rij[3][2]=E*(.5*(3.*f - 1.)*nz*ny);
  Rij[3][3]=E*(.5*(1.-f) + .5*(3.*f - 1.)*nz*nz);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns rad primitives for an atmosphere
int
set_radatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype)
{
#ifdef RADIATION  
  if(atmtype==0) //fixed Erf, urf of normal observer
    {
      pp[6]=ERADATMMIN; 
      ldouble ucon[4];
      calc_normalobs_4vel(GG,ucon);
      conv_vels(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);
      pp[7]=ucon[1]; 
      pp[8]=ucon[2];
      pp[9]=ucon[3];

    }
  if(atmtype==1) //fixed Erf, urf 0 in lab frame
    {
      ldouble ucon[4];
      ldouble xx2[4];
      ldouble GGBL[4][5];

      // BL coords
      coco_N(xx,xx2,MYCOORDS,BLCOORDS);
      calc_G_arb(xx2,GGBL,BLCOORDS);

      // normal observer in BL = stationary observer
      calc_normalobs_4vel(GGBL,ucon);
     
      // to MYCOORDS
      trans2_coco(xx2,ucon,ucon,BLCOORDS,MYCOORDS);
     
      // to VELPRIMRAD
      conv_vels(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);
     
      pp[7]=ucon[1];
      pp[8]=ucon[2];
      pp[9]=ucon[3];

      //    print_4vector(ucon); getchar();
      pp[6]=ERADATMMIN; 
     }
  if(atmtype==2) //optically thin atmosphere, scalings from numerical solution of radiall influx
    {
      ldouble gammamax=10.;
      //ldouble xxout[4]={0.,get_x(NX-1,0),0.,0.};
      //      coco_N(xxout,xxout,MYCOORDS,BLCOORDS);
      //      ldouble rout=xxout[1];

      ldouble rout=2.; //normalize at r_BL=2

      ldouble xxBL[4];
      coco_N(xx,xxBL,MYCOORDS,BLCOORDS);
      ldouble r=xxBL[1];
     
      pp[6]=ERADATMMIN*(rout/r)*(rout/r)*(rout/r)*(rout/r);

      ldouble ut[4]={0.,-gammamax*pow(r/rout,1.),0.,0.};

      ldouble ggBL[4][5],GGBL[4][5];
      calc_g_arb(xxBL,ggBL,KERRCOORDS);
      calc_G_arb(xxBL,GGBL,KERRCOORDS);

      conv_vels(ut,ut,VELR,VEL4,ggBL,GGBL);

      trans2_coco(xxBL,ut,ut,KERRCOORDS,MYCOORDS);

      conv_vels(ut,ut,VEL4,VELPRIM,gg,GG);
      
      pp[7]=ut[1];      
      pp[8]=ut[2];      
      pp[9]=ut[3];

    }
#endif
  return 0;
}

//**********************************************************************
//suplementary routines for conversions
//**********************************************************************
ldouble calc_PEQ_ufromTrho(ldouble T,ldouble rho)
{
  ldouble p=K_BOLTZ*rho*T/MU_GAS/M_PROTON;	      
  ldouble u=p/(GAMMA-1.);
  return u;
}

ldouble calc_PEQ_Tfromurho(ldouble u,ldouble rho)
{
  ldouble p=u*(GAMMA-1.);
  ldouble T=p/(K_BOLTZ*rho/MU_GAS/M_PROTON);
  return T;
}

ldouble calc_LTE_EfromT(ldouble T)
{
  return 4.*SIGMA_RAD*T*T*T*T;
}

ldouble calc_LTE_TfromE(ldouble E )
{
  return sqrt(sqrt((E/4./SIGMA_RAD)));
}


ldouble calc_LTE_Efromurho(ldouble u,ldouble rho)
{
  ldouble p=(GAMMA-1.)*(u);
  ldouble T=p*MU_GAS*M_PROTON/K_BOLTZ/rho;

  return calc_LTE_EfromT(T);
}

/************************************************************************/
/******* calculates wavespeeds in the lab frame takin 1/@3 in ************/
/******* radiative rest frame and boosting it to lab frame *****************/
/******* using the HARM algorithm - with taul limiter ********************/
/************************************************************************/
int
calc_rad_wavespeeds(ldouble *pp,void *ggg,ldouble tautot[3],ldouble *aval,int verbose)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  int i,j;

#ifdef LABRADFLUXES
  //artificially puts pp=uu and converts them to urf and Erf using the regular converter
  u2p_rad_urf(pp,pp,ggg,&i);
#endif
  
  //relative four-velocity
  ldouble urfcon[4];

#ifdef EDDINGTON_APR //in Eddington I take the fluid velocity - think over, maybe use max(cs2,1/3)?
  urfcon[0]=0.;
  urfcon[1]=pp[2];
  urfcon[2]=pp[3];
  urfcon[3]=pp[4];
  //converting to lab four-velocity
  conv_vels(urfcon,urfcon,VELPRIM,VEL4,gg,GG);
#else 
  urfcon[0]=0.;
  urfcon[1]=pp[7];
  urfcon[2]=pp[8];
  urfcon[3]=pp[9];
  //converting to lab four-velocity
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
#endif

  

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
      
#ifdef SKIPRADWAVESPEEDLIMITER
      rv2=rv2rad;
#endif

#ifdef FULLRADFRAMEWAVESPEED
      rv2=1.;
#endif

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

      aval[dim*2+0]=axl;
      aval[dim*2+1]=axr;
    }

  return 0;
}


/************************************************************************/
/******* aplying radiative four velocity-related changes in conserved ******/
/************************************************************************/
int
apply_rad_source_del4(int ix,int iy,int iz,ldouble *del4)
{
  ldouble delapl[NV];

  delapl[0]=0.; //zeros go to density and entropy so don't bother about the gdet there
  delapl[1]=-del4[0];
  delapl[2]=-del4[1];
  delapl[3]=-del4[2];
  delapl[4]=-del4[3];
  delapl[5]=0.;
  delapl[6]=del4[0];
  delapl[7]=del4[1];
  delapl[8]=del4[2];
  delapl[9]=del4[3];

  int iv;
  for(iv=0;iv<NV;iv++)
    {
      set_u(u,iv,ix,iy,iz, get_u(u,iv,ix,iy,iz)+delapl[iv] );
    }

  return 0;
}


/************************************************************************/
/******* explicit radiative source term  ***********************************/
/************************************************************************/
int explicit_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5])
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEEXPLICIT);   

  ldouble del4[4],delapl[NV];
  int iv;

  //new primitives before the source operator
  //skipped - there is one common call in finite.c
  //calc_primitives(ix,iy,iz);

  //applied explicitly directly in lab frame
  solve_explicit_lab(ix,iy,iz,dt,del4);
  indices_21(del4,del4,gg);

  apply_rad_source_del4(ix,iy,iz,del4);

  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}

/************************************************************************/
/******* implicit radiative source term in fluid frame and transported to lab  - backup method */
/************************************************************************/
int implicit_ff_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5],ldouble tlo[][4], ldouble tup[][4],ldouble *pp)
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITFF); 

  ldouble del4[4],delapl[NV];
  int iv;

  //skipped - there is one common call in finite.c
  //calc_primitives(ix,iy,iz);
  
  if(solve_implicit_ff(ix,iy,iz,dt,del4)<0) 
    {
      //failure, keeping u[] intact, reporting
      set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,-1); 
      return -1;
    }

  trans2_on2cc(del4,del4,tlo);
  boost2_ff2lab(del4,del4,pp,gg,GG);
  indices_21(del4,del4,gg);

  apply_rad_source_del4(ix,iy,iz,del4);

  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}

/************************************************************************/
/******* calculated Gi, compares with conserved and determines **********/
/******* if implicit step necessary ****************************************/
/************************************************************************/
int test_if_rad_implicit(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5], ldouble *del4)
{
  int iv;
  ldouble delapl[NV],uu[NV],pp[NV],uu0[NV];
 
  int method=0;
  
  if(method==0) //checks if inversion succesful and then max of du/u < DULIMIT
    {
      //calculating radforce
      //new primitives
      calc_primitives(ix,iy,iz);
      //rad-for-force
      solve_explicit_lab(ix,iy,iz,dt,del4);
      //del4[] will be passed up
      indices_21(del4,del4,gg); 
      //changes to conserved
      delapl[0]=0.;
      delapl[1]=-del4[0];
      delapl[2]=-del4[1];
      delapl[3]=-del4[2];
      delapl[4]=-del4[3];
      delapl[5]=0.;
      delapl[6]=del4[0];
      delapl[7]=del4[1];
      delapl[8]=del4[2];
      delapl[9]=del4[3];

      //gettin' pp & uu
      for(iv=0;iv<NV;iv++)
	{
	  uu0[iv]=get_u(u,iv,ix,iy,iz);
	  pp[iv]=get_u(p,iv,ix,iy,iz);
	  uu[iv]=uu0[iv]+delapl[iv];
	}

      struct geometry geom;
      fill_geometry(ix,iy,iz,&geom);
  
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=get_u(u,iv,ix,iy,iz);
	  pp[iv]=get_u(p,iv,ix,iy,iz);
	}

      //converting to primitives
      int corrected, fixups[2];
      u2p(uu,pp,&geom,&corrected,fixups);

      if(corrected!=0) return 1; //must do implicit
      
      //comparing with conserved to get the largest change
      ldouble maxdu=-1., uval;
      ldouble DULIMIT=0.1;

      for(iv=1;iv<NV;iv++)
	{
 	  if(iv==5) continue; //skip entropy
	  if(iv==3 && NY==1) continue; //skip y-momentum
	  if(iv==8 && NY==1) continue; //skip y-momentum
	  if(iv==4 && NZ==1) continue; //skip z-momentum
	  if(iv==9 && NZ==1) continue; //skip z-momentum

	  uval=uu[iv];

	  if(fabs(uval)<SMALL) //to avoid dividing by 0
	    maxdu=my_max(maxdu,1./SMALL);
	  else
	    {
	      maxdu=my_max(maxdu,fabs(delapl[iv]/uval));
	    }
	}

      if(maxdu<DULIMIT)
	return 0; //can do explicit
      else
	return 1; //must do implicit

    }

  if(method==1) //basing on maximal tautot
    {
      ldouble TAULIMIT=0.1;
      ldouble dx[3],xx[4];
      get_xx(ix,iy,iz,xx);
      dx[0]=get_size_x(ix,0)*sqrt(gg[1][1]);
      dx[1]=get_size_x(iy,1)*sqrt(gg[2][2]);
      dx[2]=get_size_x(iz,2)*sqrt(gg[3][3]);
      ldouble tautot[3],taumax;
      calc_tautot(pp,xx,dx,tautot);
      taumax=my_max(tautot[0],my_max(tautot[1],tautot[2]));
      if(taumax<TAULIMIT) 
	return 0; //can do explicit
      else
	return 1; //must do implicitxs
    }

  if(method==2) //based on comparison of Gi[i] with uu[i]
    {
      //calculating radforce
      //new primitives
      calc_primitives(ix,iy,iz);
      //rad-for-force
      solve_explicit_lab(ix,iy,iz,dt,del4);
      //del4[] will be passed up
      indices_21(del4,del4,gg); 

      ldouble DULIMIT=0.1;

      //gettin' pp & uu
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=get_u(u,iv,ix,iy,iz);
	  pp[iv]=get_u(p,iv,ix,iy,iz);
	}

      //changes to conserved
      delapl[0]=0.;
      delapl[1]=-del4[0];
      delapl[2]=-del4[1];
      delapl[3]=-del4[2];
      delapl[4]=-del4[3];
      delapl[5]=0.;
      delapl[6]=del4[0];
      delapl[7]=del4[1];
      delapl[8]=del4[2];
      delapl[9]=del4[3];
  
      //comparing with conserved to get the largest change
      ldouble maxdu=-1., uval;

      for(iv=1;iv<NV;iv++)
	{
 
	  if(iv==5) continue; //skip entropy
	  if(iv==3 && NY==1) continue; //skip y-momentum
	  if(iv==8 && NY==1) continue; //skip y-momentum
	  if(iv==4 && NZ==1) continue; //skip z-momentum
	  if(iv==9 && NZ==1) continue; //skip z-momentum

	  uval=uu[iv];

	  if(fabs(uval)<SMALL) //to avoid dividing by 0
	    maxdu=my_max(maxdu,1./SMALL);
	  else
	    {
	      maxdu=my_max(maxdu,fabs(delapl[iv]/uval));
	    }
	}

      //debug
      //printf("%d %d %d: if implicit = %e | %d\n",ix,iy,iz,maxdu,maxdu < DULIMIT ? 0 : 1); getchar();

      if(maxdu<DULIMIT)
	return 0; //can do explicit
      else
	return 1; //must do implicit

    }

  return 1;

}
  

/************************************************************************/
/******* explicit radiative source term with sub-stepping - inefficient *****/
/************************************************************************/
int explicit_substep_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5])
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEEXPLICITSUBSTEP); 

  int iv;
  ldouble del4[4],delapl[NV],uu[NV],pp[NV];
  double fdt, fdta, maxfu=-1., fu, uval,futau;
  //calculating reference time step basing on maximal tautot
  ldouble dx[3],xx[4];
  get_xx(ix,iy,iz,xx);
  dx[0]=get_size_x(ix,0)*sqrt(gg[1][1]);
  dx[1]=get_size_x(iy,1)*sqrt(gg[2][2]);
  dx[2]=get_size_x(iz,2)*sqrt(gg[3][3]);
  ldouble tautot[3],taumax;
  calc_tautot(pp,xx,dx,tautot);
  taumax=my_max(tautot[0],my_max(tautot[1],tautot[2]));
  futau=1./taumax;
  //reference time step only approximate
  fdt=fdta=0.;
  do
    {
      //new primitives
      //TODO: check here and there if worked
      calc_primitives(ix,iy,iz);
      //conserved
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=get_u(u,iv,ix,iy,iz);
	  pp[iv]=get_u(p,iv,ix,iy,iz);
	}

      //vector of changes of conserved assuming original dt which only multiplies source terms
      solve_explicit_lab(ix,iy,iz,dt,del4);
      indices_21(del4,del4,gg);
      //changes to conserved
      delapl[0]=0.;
      delapl[1]=-del4[0];
      delapl[2]=-del4[1];
      delapl[3]=-del4[2];
      delapl[4]=-del4[3];
      delapl[5]=0.;
      delapl[6]=del4[0];
      delapl[7]=del4[1];
      delapl[8]=del4[2];
      delapl[9]=del4[3];
		
#if(0) //my old dtsub esitmation based on single dimension
      //comparing with conserved to get the largest change
      maxfu=-1.;

      //fluxes can be zero and their relative change large
      //so far considering only energy densities
      for(iv=1;iv<NV;iv++)
	{
	  if(iv==5) continue;
	  uval=get_u(u,iv,ix,iy,iz);
	  if(fabs(uval)<SMALL)  //to avoid dividing by 0
	    fu=futau;
	  else
	    fu=fabs(delapl[iv]/uval);

	  if(ix==NX/2 ) printf("> %d %e %e %e\n",iv,fu,uval,delapl[iv]);
	  if(fu>maxfu) maxfu=fu;
	}
		  
      if(maxfu<MAXEXPLICITSUBSTEPCHANGE)
	fdt=1.;
      else
	fdt=MAXEXPLICITSUBSTEPCHANGE/maxfu;

      if(ix==NX/2 ) printf("----\n%e\n",maxfu);
		  
#else //Jon's spacetime
      //substep
      ldouble Umhd,Urad,Gtot,iUmhd,iUrad,idtsub,dtsub;
      Umhd=Urad=Gtot=0.;
      for(iv=0;iv<4;iv++)
	{
	  Umhd+=uu[1+iv]*uu[1+iv]*GG[iv][iv]; //GG?
	  Urad+=uu[6+iv]*uu[6+iv]*GG[iv][iv]; //GG?
	  Gtot+=del4[iv]*del4[iv]*GG[iv][iv]; //GG?
	}

      iUmhd=1.0/(fabs(Umhd)+SMALL);
      iUrad=1.0/(fabs(Urad)+SMALL);
      idtsub=SMALL+fabs(Gtot*my_max(iUmhd,iUrad));
      dtsub=1./idtsub;
      /*
	if(ix==42 && pp[7]!=0 && 0)
	{
	print_Nvector(pp,NV);
	print_4vector(del4);
	printf("----\n%e %e %e > %e\n",Gtot,Umhd,Urad,dtsub);getchar();
	}
      */
      if(dtsub>1.)
	fdt=1.;
      else
	fdt=dtsub;

      if(fdt<1.e-6 && 0) 
	{
	  printf("----\n%e %e %e > %e\n",Gtot,Umhd,Urad,dtsub);
	  printf("%d %d %d %f %f\n",ix,iy,iz,fdt,fdta);
	  getchar();
	}
#endif		 
		  
      if(fdta+fdt>1.) fdt=1.-fdta;
	
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]+=delapl[iv]*fdt;
	  set_u(u,iv,ix,iy,iz, uu[iv]);
	}

      fdta+=fdt;

    }
  while(fdta<1.);
	    
  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}

/************************************************************************/
/******* implicit radiative source term in lab frame *********************/
/******* solves numerically in 4D **************************************/
/************************************************************************/
int implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt, ldouble gg[][5], ldouble GG[][5],ldouble tlo[][4], ldouble tup[][4],ldouble *pp)
{
  ldouble del4[4],delapl[NV];
  int iv;
  int verbose=1;
  //  if(ix==63) verbose=1;

  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITLAB); 

  //test
  //if(ix<5) return 0;

  if(solve_implicit_lab(ix,iy,iz,dt,del4,0)<0)
      {
      set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,-1);
      //numerical implicit in 4D did not work
      if(verbose) 
	{
	  printf("===\nimp_lab didn't work at %d %d %d (%f %f %f)\ntrying imp_ff... ",ix,iy,iz,get_x(ix,0),get_x(iy,1),get_x(iz,1));
	  solve_implicit_lab(ix,iy,iz,dt,del4,1);
	  getchar();
	}
      //use the explicit-implicit backup method
      
      //test
      if(implicit_ff_rad_source_term(ix,iy,iz,dt,gg,GG,tlo,tup,pp)<0)
	{
	  if(verbose) printf("imp_ff didn't work either. requesting fixup.\n");
	  //this one failed too - failure
	  return -1;	  
	}
      else
	{
	  if(verbose) printf("worked.\n");
	  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITFF); 
	}	    
    }
  else
    {
      //success in lab frame
      set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 
      apply_rad_source_del4(ix,iy,iz,del4);
    }


  return 0;
}

/*************************************************************************/
/****** radiative lab M1-like primitives to lab Edd-like primitives*******/
/*************************************************************************/
int prad_m12edd(ldouble *pp1, ldouble *pp2, void* ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

#ifndef MULTIRADFLUID
  ldouble Rij[4][4],uufake[NV];
  int i,j;

  print_Nvector(pp1,NV);
  //to ortonormal fluid frame
  prad_lab2ff(pp1,pp1,ggg);

  print_Nvector(pp1,NV);
  //now set 1/3 on the diagonal and move back to lab frame and do u2p_rad()

  Rij[0][0]=pp1[6];
  Rij[0][1]=Rij[1][0]=pp1[7];
  Rij[0][2]=Rij[2][0]=pp1[8];
  Rij[0][3]=Rij[3][0]=pp1[9];
  Rij[1][1]=Rij[2][2]=Rij[3][3]=1./3.*Rij[0][0];
  Rij[1][2]=Rij[2][1]=Rij[1][3]=Rij[3][1]=Rij[2][3]=Rij[3][2]=0.;

  trans22_on2cc(Rij,Rij,tlo);  
  boost22_ff2lab(Rij,Rij,pp1,gg,GG); 
  indices_2221(Rij,Rij,gg);  

  uufake[6]=Rij[0][0];
  uufake[7]=Rij[0][1];
  uufake[8]=Rij[0][2];
  uufake[9]=Rij[0][3];

  //  print_Nvector(uufake,NV);
  //  getchar();

  u2p_rad(uufake,pp2,ggg,&i);
#else
  my_err("Edd does not work with multifluids yet\n");
#endif

  return 0;
} 

/***************************************/
/***************************************/
/***************************************/
int calc_rad_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,ldouble *vdiff2ret)
{  
#if(RADVISCOSITY==SHEARVISCOSITY) //full shear tensor
  struct geometry *geom
    = (struct geometry *) ggg;

  //calculating shear
  calc_shear_lab(geom->ix,geom->iy,geom->iz,shear,1);  
  indices_1122(shear,shear,geom->GG);

  //transforming to ortonormal
  ldouble shearon[4][4];
  trans22_cc2on(shear,shearon,geom->tup);

  //calculating the viscosity coefficient 

  //radiative energy density - so far in the radiative frame!
  ldouble Erf=pp[EE0];
  //mean free path
  ldouble chi=calc_chi(pp,geom->xxvec);
  ldouble mfp = 1./chi;

  //limiting in opt.thin region
  ldouble dx[3]={get_size_x(geom->ix,0)*sqrt(geom->gg[1][1]),   //here gg can be face or cell, get_size_x always refers to cell
		 get_size_x(geom->iy,1)*sqrt(geom->gg[2][2]),
		 get_size_x(geom->iz,2)*sqrt(geom->gg[3][3])};
  ldouble mindx;
  if(NY==1 && NZ==1) mindx = dx[0];
  else if(NZ==1) mindx = my_min(dx[0],dx[1]);
  else if(NY==1) mindx = my_min(dx[0],dx[2]);
  else mindx = my_min(dx[0],my_min(dx[1],dx[2]));
  if(mfp>mindx || chi<SMALL) mfp=mindx;

  ldouble ev[4],evmax,eta,nu,vdiff2;
  nu = ALPHARADVISC * 1./3. * mfp;
  
  //no longer necessary?
  if(PROBLEM==30 || PROBLEM==43 || PROBLEM==54) //RADNT to overcome huge gradients near fixed radiative atmosphere at r>rout
    if(geom->ix>=NX-2)
      nu = 0.; 

  //limiting basing on diffusive wavespeed
  ldouble MAXDIFFVEL=0.5; //max allowed vdiff
  ldouble MAXTOTVEL=0.75; //max allowed vdiff + vrad

  //limiting basing on maximal eigen value - slower and issues with tetrad  
  evmax=calc_eigen_4x4(shearon,ev);

  //limiting assuming maximal eigen value 1/dt
  //evmax=1./dt;

  //square of characteristic velocity for diffusion
  vdiff2=2.*nu*evmax;

  //checking if vdiff too large
  if(vdiff2 > MAXDIFFVEL*MAXDIFFVEL)
    {
      nu = MAXDIFFVEL*MAXDIFFVEL/2./evmax;
    }
  vdiff2=2.*nu*evmax;

  /*
  //checking if vdiff+vrad > 1
  if(vrad>MAXTOTVEL*MAXTOTVEL)
    {
      nu=0.;
      vdiff2=0.;
    }
  else if(sqrt(vrad)+sqrt(vdiff2)>MAXTOTVEL)
    {
      vdiff2=MAXTOTVEL-sqrt(vrad);
      vdiff2*=vdiff2;
      nu=vdiff2/2./evmax;
    }
  */
  
  *nuret=nu;
  *vdiff2ret=vdiff2;

  return 0;
#endif
  
}

 

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
  uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
  uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
  uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
  uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);

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

  f[0] = uu[EE0] - uu0[EE0] + dt * gdetu * Gi[0];
  f[1] = uu[FX0] - uu0[FX0] + dt * gdetu * Gi[1];
  f[2] = uu[FY0] - uu0[FY0] + dt * gdetu * Gi[2];
  f[3] = uu[FZ0] - uu0[FZ0] + dt * gdetu * Gi[3];

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

  //printf("%.20e %.20e %.20e \n",uu[EE0],Gi[0],uu[EE0] - uu0[EE0] + dt * gdetu * Gi[0]);
  //printf("%.20e %.20e %.20e %.20e %.20e %.20e \n",uu[EE0],Rij[0][0],pp[EE0],Rtt,Ehat-4.*Pi*B,Rtt - Rtt0 - kappaabs*(Ehat-4.*Pi*B)*dtau);
  //f[0]=Rtt - Rtt0 - kappaabs*(Ehat-4.*Pi*B)*dtau;
  
  return 0;
} 

int
print_state_implicit_lab_4dcon(int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .13e % .13e % .13e % .3e "
	  "f(x) = % .13e % .13e % .13e % .13e\n",
	  iter,
	  //	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
	  x[0],x[1],x[2],x[3],f[0],f[1],f[2],f[3]);
  return 0;
}

int
solve_implicit_lab_4dcon(ldouble *uu00,ldouble *pp00,void *ggg,ldouble dt,ldouble* deltas,int verbose)
{
  int i1,i2,i3,iv,i,j;
  ldouble J[4][4],iJ[4][4];
  ldouble pp[NV],pp0[NV],uu[NV],uu0[NV],uup[NV]; 
  ldouble f1[4],f2[4],f3hd[4],f3rad[4],xxx[4];

  ldouble (*gg)[5],(*GG)[5];

  struct geometry *geom
    = (struct geometry *) ggg;

  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;

  //temporary using local arrays
  gg=geom->gg;
  GG=geom->GG;

  int corr[2],fixup[2];
  u2p(uu00,pp00,geom,corr,fixup);
  p2u(pp00,uu00,geom);

  for(iv=0;iv<NV;iv++)
    {
      uu0[iv]=uu00[iv]; 
      pp0[iv]=pp00[iv]; 
      uu[iv]=uu0[iv]; 
      pp[iv]=pp0[iv];     
    }

  ldouble EPS = 1.e-8;
  ldouble CONV = 1.e-6; 
  int MAXITER=50;
  ldouble DAMP = 0.5;

  ldouble frdt = 1.0;
  ldouble dttot = 0.;

  int iter=0;
  int failed=0;

 
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
		  xxx[i]=uup[i+EE0];
		}  
	      print_Nvector(uu0,NV);
	      print_Nvector(uu,NV);
	      print_Nvector(pp0,NV);

	      int ret=f_implicit_lab_4dcon(uu0,uu,pp0,frdt*(1.-dttot)*dt,geom,f1);
	      print_state_implicit_lab_4dcon (iter-1,xxx,f1); 
	      printf("f_lab_4dcon ret: %d\n",ret);
	    }


	  //values at base state
	  if(f_implicit_lab_4dcon(uu0,uu,pp0,frdt*(1.-dttot)*dt,geom,f1)<0) 
	    {
	      failed=1;
		  
	      break;
	    }
	  
	  //calculating approximate Jacobian
	  for(j=0;j<4;j++)
	    {
	      ldouble del;

	      del=EPS*uup[EE0]; 

	      uu[j+EE0]=uup[j+EE0]-del;
	      
	      int fret=f_implicit_lab_4dcon(uu0,uu,pp0,frdt*(1.-dttot)*dt,geom,f2);  

	      if(verbose>0)
		{
		  for(i=0;i<4;i++)
		    {
		      xxx[i]=uu[i+EE0];
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
		  J[i][j]=(f2[i] - f1[i])/(uu[j+EE0]-uup[j+EE0]);
		}

	      uu[j+EE0]=uup[j+EE0];

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
	      xxx[i]=uup[i+EE0];
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
	      uu[i+EE0]=xxx[i];
	    }
  
	   //opposite changes in gas quantities
	  uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
	  uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
	  uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
	  uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);

	  //test convergence
	  for(i=0;i<4;i++)
	    {
	      f3rad[i]=(uu[i+EE0]-uup[i+EE0]);
	      f3hd[i]=(uu[i]-uup[i]);
	      
	      f3rad[i]=fabs(f3rad[i]/my_max(fabs(uup[EE0]),fabs(uup[i])));
	      f3hd[i]=fabs(f3hd[i]/my_max(fabs(uup[1]),fabs(uup[0])));
	    }

	  if(f3rad[0]<CONV && f3rad[1]<CONV && f3rad[2]<CONV && f3rad[3]<CONV)
	    {
	      if(verbose) printf("success ===\n");
	      break;
	    }

	  if(iter>MAXITER)
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
	  uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
	  uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
	  uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
	  uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);

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
	  u2p(uu0,pp0,geom,corr,fixup);
	  p2u(pp0,uu0,geom);
	  continue;
	}
      
      //didn't work - decreasing time step
      frdt *= DAMP;

      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=uu0[iv];
	}
      
      if(frdt<0.00001)  //avoid substepping?
	{
	  if(verbose) 
	    {
	      printf("time step too small - aborting implicit_lab_4dcon() ===\n");
	    }
	  return -1;
	}
    }
  while(1);
  
  deltas[0]=uu[EE0]-uu00[EE0];
  deltas[1]=uu[FX0]-uu00[FX0];
  deltas[2]=uu[FY0]-uu00[FY0];
  deltas[3]=uu[FZ0]-uu00[FZ0];

  if(verbose) 
    {
      u2p(uu,pp0,geom,corr,fixup);
      print_4vector(deltas);
      print_NVvector(uu);
      ldouble T=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO]);
      printf("Tgas: %e\n",T);
    }
  
  return 0;
}




/*************************************************************************/
/******* 1D solver in energy densities ***********************************/
/******* totenergy is the total energy density in ff *********************/
/******* ratio is the radio of energy denisities between ff and rad.rest frame */
/*************************************************************************/
ldouble
calc_implicit_1dprim_err(ldouble xxx,ldouble *uu0,ldouble *pp0,ldouble dtau,void *ggg,int *params,ldouble totenergy,ldouble ratio)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  int whichprim=params[0];

  // printf("input: %e %e %e %e\n",xxx,pp0[ENTR],dt,totenergy);
  
  ldouble ugas,Ehat,rho,Rtt,ucon[4],backup,err;

  rho=pp0[RHO];

  if(whichprim==RAD)
    {
      Ehat=xxx*ratio;
      ugas=totenergy-Ehat;
    }

  if(whichprim==MHD)
    {
      ugas=xxx;
      Ehat=totenergy-ugas;      
    }
       
  ldouble T=calc_PEQ_Tfromurho(ugas,rho);
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
  ldouble kappaabs=calc_kappa(rho,T,geom->xx,geom->yy,geom->zz);
  ldouble entr=calc_Sfromu(rho,ugas);
  ldouble entr0=pp0[ENTR];

  err = entr - entr0 - kappaabs*(Ehat-4.*Pi*B)*dtau;

  //  printf("enden in: %e %e %e\n",Ehat,ugas,totenergy);
  //  printf("output: %e %e %e %e %e\n",entr,entr0,Ehat,4.*Pi*B,dtau);getchar();

  return err;
}

int
solve_implicit_lab_1dprim(ldouble *uu0,ldouble *pp0,void *ggg,ldouble dt,ldouble* deltas,int verbose,ldouble *ppout)
{
  int i1,i2,i3,iv,i,j,iter;
  ldouble pp[NV],uu[NV]; 
  ldouble err[3],xxx[4],xxx0[4],err4d[4];
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  struct geometry *geom
    = (struct geometry *) ggg;

  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;

  //temporary using local arrays
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  //ff rad energy density
  ldouble Rtt0,Ehat,ugas0[4];
  calc_ff_Rtt(pp0,&Rtt0,ugas0,geom);
  Ehat=-Rtt0;

  //time step in fluid frame
  ldouble dtau;
  dtau=dt/ugas0[0];

  //make sure entropy in pp0 is consistent
  pp0[ENTR]=calc_Sfromu(pp0[RHO],pp0[UU]);

  //total ff energy
  ldouble totenergy;
  totenergy=Ehat + pp0[UU];

  if(verbose) printf("starting 1D \n");
  if(verbose) print_NVvector(pp0);
  if(verbose) printf("enden: %e %e %e\n",Ehat,pp0[UU],totenergy);
  
  //residual function parameters
  int params[3],whichprim;
  if(-Rtt0<1.e-3*pp0[UU]) //hydro preffered
    whichprim=RAD;
  else
    whichprim=MHD; 

  //override
  //whichprim=MHD;

  //solve in rad or mhd temperature?
  params[0]=whichprim;
  //energy or entropy equation to solve
  params[1]=RADIMPLICIT_ENTROPYEQ;
  //frame for energy/entropy equation to solve
  params[2]=RADIMPLICIT_FFEQ;

  //local vectors
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=uu0[iv]; 
      pp[iv]=pp0[iv];     
    }

  //initiate iterated quantity
  if(whichprim==MHD)
    {
      xxx[0]=pp0[UU];      
      xxx[1]=pp0[UU];
      xxx[2]=Ehat;

      //test
      xxx[2]=pp0[UU];

      if(verbose) printf("working on MHD\n");
    }

  ldouble ratio;
  if(whichprim==RAD)
    {
      ratio = Ehat/pp0[EE0]; //ratio of energy densities between the ff and rad rest frame
      xxx[0]=pp0[EE0];    
      xxx[1]=pp0[EE0]; 
      xxx[2]=pp0[UU]/ratio;

      //test
      xxx[2]=pp0[EE0];

      if(verbose) printf("working on RAD\n");
    }

  //order the lower-upper brackets 
  if(xxx[1]>xxx[2]) 
    {
      err[0]=xxx[1];
      xxx[1]=xxx[2];
      xxx[2]=err[0];
    }

  xxx0[0]=xxx[0];
  xxx0[1]=xxx[1];
  xxx0[2]=xxx[2];
  
  //test signs
  err[0]=calc_implicit_1dprim_err(xxx[0],uu0,pp0,dt,geom,params,totenergy,ratio);
  err[1]=calc_implicit_1dprim_err(xxx[1],uu0,pp0,dt,geom,params,totenergy,ratio);
  err[2]=calc_implicit_1dprim_err(xxx[2],uu0,pp0,dt,geom,params,totenergy,ratio);

  if(verbose) printf("%d >>> [%e : %e : %e] > ",0,xxx[1],xxx[0],xxx[2]);
  if(verbose) printf("[%e : %e : %e]\n",err[1],err[0],err[2]);  

  if(verbose) printf("searching for valid brackets\n");

  //first try low the lower bracket
  int MAXITER=30;
  ldouble FAC=2.;
  ldouble CONV=1.e-5; //cannot be too low!
  
  iter=0.;
  while(err[1]*err[2]>0.)
    {
      iter++;
      xxx[1]/=FAC;
      err[1]=calc_implicit_1dprim_err(xxx[1],uu0,pp0,dt,geom,params,totenergy,ratio);

      /*
      if(whichprim==RAD)
	xxx[2]=xxx[2]+(totenergy/ratio-xxx[2])/2.;
      if(whichprim==MHD)
	xxx[2]=xxx[2]+(totenergy-xxx[2])/2.;
      */

      //make sure this works all the time
      xxx[2]*=FAC;
      if(whichprim==RAD && xxx[2]>=totenergy/ratio)
	xxx[2]=(1.-1.e-10)*totenergy/ratio;
      if(whichprim==MHD && xxx[2]>=totenergy)
	xxx[2]=(1.-1.e-10)*totenergy;
      
      err[2]=calc_implicit_1dprim_err(xxx[2],uu0,pp0,dt,geom,params,totenergy,ratio);

      if(verbose) printf("%d (%d) >>> [%e : %e : %e] > ",0,iter,xxx[1],xxx[0],xxx[2]);
      if(verbose) printf("[%e : %e : %e]\n",err[1],err[0],err[2]);  
      if(isnan(err[1]) || isnan(err[2])) iter=MAXITER;
      if(iter>=MAXITER) break;
    }

  if(iter>=MAXITER)
    {
      if(verbose || 1) {printf("brackets not found!\n");}
      getchar();
      PLOOP(i) ppout[i]=pp0[i];
      return -1;
    }      
    
  if(verbose) printf("brackets found!\n");

  //call more sophisticated solver here
  do
    {
      iter++;
      xxx[0]=0.5*(xxx[1]+xxx[2]); //new estimate
      err[0]=calc_implicit_1dprim_err(xxx[0],uu0,pp0,dt,geom,params,totenergy,ratio); //new error
      

      if(err[0]*err[1]>0.) //same sign as the lower bracket
	{
	  xxx[1]=xxx[0];
	  err[1]=err[0];
	}
      else
	{
	  xxx[2]=xxx[0];
	  err[2]=err[0];
	}
      if(verbose) printf("%d >>> [%e : %e : %e] > ",iter,xxx[1],xxx[0],xxx[2]);
      if(verbose) printf("[%e : %e : %e]\n",err[1],err[0],err[2]);  

      if(isnan(xxx[0]) || isnan(err[0])) {my_err("nan in 1dprim\n");}
    }
  while(fabs((xxx[2]-xxx[1])/xxx[0])>CONV);

  if(verbose) printf("solution found: %e\n",xxx[0]);

  //converting to new set of primitives
  PLOOP(i) pp[i]=pp0[i];
  if(whichprim==RAD)
    pp[EE0]=xxx[0];
  else
    pp[UU]=xxx[0];
  
  //is what is below required?

  //total inversion, but only whichprim part matters
  p2u(pp,uu,geom);
   
  //opposite changes in the other quantities and inversion
  int corr[2],fixup[2];
  if(whichprim==RAD)
    {
      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);  

      u2p(uu,pp,geom,corr,fixup); //total inversion (I should separate hydro from rad)
    }
  if(whichprim==MHD)
    {
      uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
      uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
      uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
      uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);

      u2p_rad(uu,pp,geom,corr);
    }     

  //updating entropy
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

  //returning the new set of primitives
  PLOOP(i) ppout[i]=pp[i];

  if(verbose) print_NVvector(ppout);
  if(verbose) getchar();

  return 0; 
}

//**********************************************************************
//******* solves implicitidly four-force source terms *********************
//******* in the lab frame  working on primitives    ***********************
//******* rad or hydro (whichprim) **************************************
//**********************************************************************

int f_implicit_lab_4dprim(ldouble *ppin,ldouble *uu0,ldouble *pp0,ldouble dt,void* ggg,ldouble *f,int *params,ldouble *err0)
{
  int ret=0,i;
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
  int whichframe=params[2];

  //  printf("%d %d %d\n",params[0],params[1],params[2]);getchar();

  ldouble uu[NV],pp[NV],err[4];
  int corr[2]={0,0},fixup[2]={0,0},u2pret,i1,i2;

  for(i=0;i<NV;i++) pp[i]=ppin[i];
  
  //rho may be inconsistent on input if iterating MHD primitives
  ldouble ucon[4];
  //correcting rho for MHD prims
  if(whichprim==MHD)
    {
      ucon[1]=pp[2];
      ucon[2]=pp[3];
      ucon[3]=pp[4];
      ucon[0]=0.;
      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
      ldouble rho = uu0[RHO]/gdetu/ucon[0];
      pp[RHO]=rho;
    }

  //total inversion, but only whichprim part matters
  p2u(pp,uu,geom);
 
  //opposite changes in the other quantities and inversion
  if(whichprim==RAD)
    {
      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);  

      u2pret=u2p(uu,pp,geom,corr,fixup); //total inversion (I should separate hydro from rad)

      //calculate ff four-velocity
      ucon[1]=pp[2];
      ucon[2]=pp[3];
      ucon[3]=pp[4];
      ucon[0]=0.;
      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
    }
  if(whichprim==MHD)
    {
      uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
      uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
      uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
      uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);

      u2pret=u2p_rad(uu,pp,geom,corr);
    }     

  //updating entropy
  pp0[ENTR]=calc_Sfromu(pp0[RHO],pp0[UU]);
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);
  uu[ENTR]=pp[ENTR]*gdetu*ucon[0];

  //print_Nvector(uu,NV);getchar();

  if(corr[0]!=0 || corr[1]!=0) 
    ret=1;
  //printf("corr: %d %d\n",corr[0],corr[1]);

  if(u2pret<-1 || u2pret<0) 
    {
      printf("implicit sub-sub-step failed\n");
      return -1; //allows for entropy but does not update conserved 
    }

  //radiative four-force
  ldouble Gi[4];
  calc_Gi(pp,ggg,Gi); 
  indices_21(Gi,Gi,gg);

  //errors in momenta - always in lab frame
  if(whichprim==MHD) //mhd-primitives
    {
      f[1] = uu[2] - uu0[2] - dt * gdetu * Gi[1];
      f[2] = uu[3] - uu0[3] - dt * gdetu * Gi[2];
      f[3] = uu[4] - uu0[4] - dt * gdetu * Gi[3];

      if(fabs(f[1])>SMALL) err[1]=fabs(f[1])/(fabs(uu[2])+fabs(uu0[2])+fabs(dt*gdetu*Gi[1])); else err[1]=0.;
      if(fabs(f[2])>SMALL) err[2]=fabs(f[2])/(fabs(uu[3])+fabs(uu0[3])+fabs(dt*gdetu*Gi[2])); else err[2]=0.;
      if(fabs(f[3])>SMALL) err[3]=fabs(f[3])/(fabs(uu[4])+fabs(uu0[4])+fabs(dt*gdetu*Gi[3])); else err[3]=0.;
    }
  if(whichprim==RAD) //rad-primitives
    {
      f[1] = uu[FX0] - uu0[FX0] + dt * gdetu * Gi[1];
      f[2] = uu[FY0] - uu0[FY0] + dt * gdetu * Gi[2];
      f[3] = uu[FZ0] - uu0[FZ0] + dt * gdetu * Gi[3];

      if(fabs(f[1])>SMALL) err[1]=fabs(f[1])/(fabs(uu[FX0])+fabs(uu0[FZ0])+fabs(dt*gdetu*Gi[1])); else err[1]=0.;
      if(fabs(f[2])>SMALL) err[2]=fabs(f[2])/(fabs(uu[FY0])+fabs(uu0[FY0])+fabs(dt*gdetu*Gi[2])); else err[2]=0.;
      if(fabs(f[3])>SMALL) err[3]=fabs(f[3])/(fabs(uu[FZ0])+fabs(uu0[FZ0])+fabs(dt*gdetu*Gi[3])); else err[3]=0.;
    }

  /***** LAB FRAME ENERGY/ENTROPY EQS *****/
  if(whichframe==RADIMPLICIT_LABEQ)
    {      
      if(whichprim==RAD) //rad-primitives
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0] = uu[EE0] - uu0[EE0] + dt * gdetu * Gi[0];
	      if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[EE0]) + fabs(uu0[EE0]) + fabs(dt * gdetu * Gi[0])); else err[0]=0.;
	    }
	  else if(whicheq==RADIMPLICIT_ENTROPYEQ)
	    {
	      f[0] = uu[ENTR] - uu0[ENTR] + dt * gdetu * Gi[0]; //but this works on hydro entropy and may fail!
	      if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[ENTR]) + fabs(uu0[ENTR]) + fabs(dt * gdetu * Gi[0])); else err[0]=0.;
	    }
	  else
	    my_err("not implemented 3\n");
	}
      if(whichprim==MHD) //hydro-primitives
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0] = uu[UU] - uu0[UU] - dt * gdetu * Gi[0];
	      if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[UU]) + fabs(uu0[UU]) + fabs(dt * gdetu * Gi[0])); else err[0]=0.;
	    }
	  else if(whicheq==RADIMPLICIT_ENTROPYEQ)
	    {	  
	      f[0] = uu[ENTR] - uu0[ENTR] - dt * gdetu * Gi[0];
	      if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[ENTR]) + fabs(uu0[ENTR]) + fabs(dt * gdetu * Gi[0])); else err[0]=0.;
	    }      
	  else
	    my_err("not implemented 4\n");
	}
    }

  /***** FF FRAME ENERGY/ENTROPY EQS *****/
  if(whichframe==RADIMPLICIT_FFEQ)
    {

      ldouble ucon[4],Rtt0,Rtt;
      //zero - state 
      calc_ff_Rtt(pp0,&Rtt0,ucon,geom);
      //new state
      calc_ff_Rtt(pp,&Rtt,ucon,geom);
      
      ldouble T=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
      ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
      ldouble Ehat = -Rtt;
      ldouble Ehat0 = -Rtt0;
      ldouble dtau=dt/ucon[0];
      ldouble kappaabs=calc_kappa(pp[RHO],T,geom->xx,geom->yy,geom->zz);

      //fluid frame energy equation:
      if(whichprim==RAD) //rad-primitives
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0]=Ehat - Ehat0 + kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(Ehat) + fabs(Ehat0) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
	    }
	  else
	    {
	      pp[ENTR]= calc_Sfromu(pp[RHO],pp[UU]);
	      f[0]=pp[ENTR] - pp0[ENTR] - kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(pp[ENTR]) + fabs(pp0[ENTR]) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
	    }
	    
	}
      if(whichprim==MHD) //mhd-
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0]=pp[UU] - pp0[UU] - kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(pp[UU]) + fabs(pp0[UU]) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
	    }
	  else if(whicheq==RADIMPLICIT_ENTROPYEQ)
	    {
	      pp[ENTR]= calc_Sfromu(pp[RHO],pp[UU]);
	      f[0]=pp[ENTR] - pp0[ENTR] - kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(pp[ENTR]) + fabs(pp0[ENTR]) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
	    }
	  else
	    my_err("not implemented 2\n");	  
	}
    }

  *err0=my_max(my_max(err[0],err[1]),my_max(err[2],err[3]));
  
  return ret;
} 

int
print_state_implicit_lab_4dprim (int iter, ldouble *x, ldouble *f,ldouble err)
{
  printf ("iter = %3d x = % .13e % .13e % .13e % .13e "
	  "e = %.10e f(x) = % .10e % .10e % .10e % .10e\n",
	  iter,
	  //	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
	  x[0],x[1],x[2],x[3],err,f[0],f[1],f[2],f[3]);
  return 0;
}

/*****************************************/
/******* 4D solver working on primitives */
/*****************************************/
int
solve_implicit_lab_4dprim(ldouble *uu00,ldouble *pp00,void *ggg,ldouble dt,ldouble* deltas,int verbose,int *params)
{
  int i1,i2,i3,iv,i,j;
  ldouble J[4][4],iJ[4][4];
  ldouble pp[NV],pp0[NV],ppp[NV],uu[NV],uu0[NV],uup[NV]; 
  ldouble f1[4],f2[4],f3[4],xxx[4],err;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  struct geometry *geom
    = (struct geometry *) ggg;

  int ix=geom->ix;
  int iy=geom->iy;
  int iz=geom->iz;

  //temporary using local arrays
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif
  
  if(verbose) //dump the problematic case to file
    {
      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
      printf("\n@@@@@@@@ IMPLICIT @@@@@@@@@@@@");
      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");

      /*
      print_Nvector(uu00,NV);
      print_Nvector(pp00,NV);
      print_metric(gg);
      print_metric(GG);
      printf("%e %e %e\n",dt,geom->alpha,geom->gdet);
      */
      FILE *out = fopen("imp.problem.dat","w");
      for (i1=0;i1<NV;i1++)
	fprintf(out,"%.15e ",uu00[i1]);
      for (i1=0;i1<NV;i1++)
	fprintf(out,"%.15e ",pp00[i1]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  fprintf(out,"%.15e ",gg[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  fprintf(out,"%.15e ",GG[i1][i2]);
      fprintf(out,"%.15e \n",dt);
      fprintf(out,"%.15e \n",geom->alpha);
      fprintf(out,"%.15e \n",geom->gdet);
      fprintf(out,"\n");
      fclose(out);
      printf("dumped problematic case to imp.problem.dat\n");
    }

  //comparing energy densities
  ldouble ugas00[4],Rtt00,Tgas00,Trad00;
  int dominates;
  calc_ff_Rtt(pp00,&Rtt00,ugas00,geom);
   if(-Rtt00>pp00[UU]) 
    dominates = RAD;
  else
    dominates = MHD;
  Tgas00=calc_PEQ_Tfromurho(pp00[UU],pp00[RHO]);
  Trad00=calc_LTE_TfromE(-Rtt00);

  ldouble ppLTE[NV],uuLTE[NV];
  //  calc_LTE_state(pp00,ppLTE,geom);
  //  ldouble TLTE=calc_PEQ_Tfromurho(ppLTE[UU],ppLTE[RHO]);

  //nolonger keep TLTE - whats below can loop up
  //ldouble TLTE=calc_LTE_temp(pp00,geom);
  ldouble TLTE=0.;


  ldouble kappa=calc_kappa(pp00[RHO],Tgas00,0.,0.,0.);
  ldouble chi=kappa+calc_kappaes(pp00[RHO],Tgas00,0.,0.,0.);
  ldouble xi1=kappa*dt*(1.+16.*SIGMA_RAD*pow(Tgas00,4.)/pp00[UU]);
  ldouble xi2=chi*dt*(1.+(-Rtt00)/(pp00[RHO]+GAMMA*pp00[UU]));
  ldouble ucon[4];

  if(verbose) 
    {
      printf("\n===========\nxi1: %e xi2: %e\n===========\n\n",xi1,xi2);
      ldouble Rtt;
      calc_ff_Rtt(pp00,&Rtt,ucon,geom);
      printf("gamma gas: %e\n\n",ucon[0]);
      ldouble Gi00[4],Gihat00[4];
      //      print_Nvector(pp00,NV);
      calc_Gi(pp00, geom,Gi00);
      //print_4vector(Gi00);
      //getchar();
      
      indices_21(Gi00,Gi00,geom->gg);

      boost2_lab2ff(Gi00,Gihat00,pp00,geom->gg,geom->GG);
      for(iv=0;iv<4;iv++)
	{
	  Gi00[iv]*=dt*gdetu;
	  Gihat00[iv]*=dt*gdetu;
	}
      //print_4vector(Gi00);
      //print_4vector(Gihat00);
      //getchar();
 

      ldouble Trad=calc_LTE_TfromE(-Rtt);
      ldouble Tgas=calc_PEQ_Tfromurho(pp00[UU],pp00[RHO]);
      ldouble rho=pp00[RHO];
      ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
      ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);  
      printf("\n===========\nkappa: %e chi: %e\n===========\n",kappa,kappa+kappaes);
      printf("\nxi1: %e xi2: %e\n",xi1,xi2);
      
      printf("gas temp: %e\n",Tgas00);      
      printf("rad temp: %e\n",Trad00); 
      printf("LTE temp: %e\n\n",TLTE) ;  
    }
  
  /******************************************/
  /******************************************/
  /******************************************/
  //choice of primitives to evolve
  int whichprim;
  if(-Rtt00<1.e-3*pp00[UU]) //hydro preffered
    whichprim=RAD;
  else
    whichprim=MHD; 
  
  params[0]=whichprim;

  //override the given parameters
  //params[0]=MHD;
  //energy or entropy equation to solve
  //params[1]=RADIMPLICIT_ENTROPYEQ;
  //params[1]=RADIMPLICIT_ENERGYEQ;
  //frame for energy/entropy equation to solve
  //params[2]=RADIMPLICIT_LABEQ;
  /******************************************/
  /******************************************/
  /******************************************/

  //check if one can compare gas & rad velocities
  if(VELPRIM!=VELPRIMRAD) 
    my_err("implicit solver assumes VELPRIM == VELPRIMRAD\n");

  for(iv=0;iv<NV;iv++)
    {
      uu0[iv]=uu00[iv]; //zero state for substepping
      pp0[iv]=pp00[iv]; 
      uu[iv]=uu0[iv]; 
      pp[iv]=pp0[iv];     
    }
 
  ldouble EPS = 1.e-8;
  ldouble CONV = 1.e-6;
  ldouble MAXITER = 50;
  int corr[2],fixup[2];

  int sh;

  if(whichprim==MHD) 
    sh=UU; //solving in hydro primitives
  else if(whichprim==RAD) 
    sh=EE0; //solving in rad primitives

  ldouble frdt = 1.0;

  int iter=0;
  int failed=0;

  if(verbose) 
    {
      //printf("=== i: %d %d %d\n\n",ix,iy,iz);
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);
      //print_metric(gg);
    }

  if(verbose) 
    {
      printf("\n===\n Trying imp lab 4d prim with dt : %f \n",dt);
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

	  int ret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f1,params,&err);
	  print_state_implicit_lab_4dprim (iter-1,xxx,f1,err); 
	  if(ret<0) printf("f_lab_4dprim ret: %d\n",ret);
	}


      //values at base state
      if(f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f1,params,&err)<0) 
	{
	  return -1;	  
	}

	  
      if(err<CONV)
	{
	  if(verbose) printf("success ===\n");
	  break;
	}


      //calculating approximate Jacobian
      for(j=0;j<4;j++)
	{
	  ldouble del;

	  //one-way derivatives
	  if(j==0)
	    {
	      //uses EPS of the dominating quantity
	      if(dominates==RAD)
		del=EPS*ppp[EE0];
	      else
		del=EPS*ppp[UU];	   
	      
	      //EPS of the iterated quantity
	      del=EPS*ppp[sh];

	      //EPS of the geometrical mean
	      del=EPS*sqrt(ppp[EE0]*ppp[UU]);
	    }
	  else //decreasing velocity
	    {
	      if(ppp[j+sh]>=0.)
		del=-EPS; 
	      else
		del=EPS;
	    }

	  pp[j+sh]=ppp[j+sh]+del;
	      
	  int fret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f2,params,&err);  

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

      //applying corrections
      ldouble xiapp=1.;
      do //check whether energy density positive and the error function returns 0
	{	    
	  //original x
	  for(i=0;i<4;i++)
	    {
	      xxx[i]=ppp[i+sh];
	    }

	  //updating x
	  for(i=0;i<4;i++)
	    {
	      for(j=0;j<4;j++)
		{
		  xxx[i]-=xiapp*iJ[i][j]*f1[j];
		}
	    }

	  //update primitives
	  for(i=0;i<4;i++)
	    {
	      pp[i+sh]=xxx[i];
	    }

	  if(whichprim==MHD)
	    {
	      //correct rho to follow new velocity (only for MHD primitives)
	      ucon[1]=pp[2];ucon[2]=pp[3];ucon[3]=pp[4]; ucon[0]=0.;
	      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);  
	      pp[RHO]=uu00[RHO]/gdetu/ucon[0];
	    }

	  //updating the other set of quantities
	  p2u(pp,uu,geom);

	  int u2pret;
	  //opposite changes in the other quantities and inversion
	  if(whichprim==RAD)
	    {
	      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
	      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
	      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
	      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);
	      u2pret=u2p(uu,pp,geom,corr,fixup); //total inversion (I should separate hydro from rad)
	      ucon[1]=pp[2]; ucon[2]=pp[3]; ucon[3]=pp[4]; ucon[0]=0.;
	      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);  
	    }
	  if(whichprim==MHD)
	    {
	      uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
	      uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
	      uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
	      uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);
	      u2pret=u2p_rad(uu,pp,geom,corr);
	    }     

	  if(xxx[0]>0. && u2pret==0) break;

	  //decrease the applied fraction
	  xiapp/=2.; //to be generous

	  if(xiapp<1.e-7) 
	    {
	      printf("damped unsuccesfully in implicit_4dprim\n");
	      return -1;	      
	    }
	}
      while(1); 

      /*
      int overshoot=0,overcnt=0;	      
      ldouble xi[4]={1.,1.,1.,1.}; //fraction of the Jacobian-implied step to apply
      ldouble xiapp=1.,fneg;
      do //overshooting
	{
	  fneg=1.;

	  do //negative energy density
	    {	    
	      //original x
	      for(i=0;i<4;i++)
		{
		  xxx[i]=ppp[i+sh];
		}

	      //updating x
	      for(i=0;i<4;i++)
		{
		  for(j=0;j<4;j++)
		    {
		      xxx[i]-=fneg*xiapp*iJ[i][j]*f1[j];
		    }
		}

	      if(xxx[0]>0.) break;

	      //negative en.density 
	      fneg = ppp[sh] / (ppp[sh]-xxx[0]);
	      fneg/=1.5; //to be generous
	    }
	  while(1); //neg.energy

	  if(verbose>0) print_state_implicit_lab_4dprim (iter,xxx,f1); 

	  overshoot=0;

	  for(i=0;i<4;i++)
	    {
	      pp[i+sh]=xxx[i];
	    }

	  if(whichprim==MHD)
	    {
	      //correct rho to follow new velocity (only for MHD primitives)
	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];
	      ucon[0]=0.;
	      //converting to 4-velocity
	      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);  
	      ldouble rho = uu00[RHO]/gdetu/ucon[0];
	      pp[RHO]=rho;

	      if(verbose) printf("rho: %e\n",rho);
	    }

	  //updating the other set of quantities
	  //total inversion, but only whichprim part matters
	  p2u(pp,uu,geom);
	  //opposite changes in the other quantities and inversion
	  if(whichprim==RAD)
	    {
	      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
	      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
	      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
	      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);
	      u2p(uu,pp,geom,corr,fixup); //total inversion (I should separate hydro from rad)
	      //calc the ff 4-velocity
	      ucon[1]=pp[2];
	      ucon[2]=pp[3];
	      ucon[3]=pp[4];
	      ucon[0]=0.;
	      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);  
	    }
	  if(whichprim==MHD)
	    {
	      uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
	      uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
	      uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
	      uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);
	      u2p_rad(uu,pp,geom,corr);
	    }     

	  //overshooting check
	  ldouble ucon[4],Rtt,Tgas,Trad;

	  //new temperatures
	  calc_ff_Rtt(pp,&Rtt,ucon,geom);
	  Tgas=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
	  Trad=calc_LTE_TfromE(-Rtt);
	  //calc_LTE_state(pp,ppLTE,geom); 
	  //TLTE=calc_PEQ_Tfromurho(ppLTE[UU],ppLTE[RHO]);
	  //TLTE=calc_LTE_temp(pp,geom);
	  TLTE=0.;

	  if(verbose) 
	    {
	      printf("\nTgas: %.20e Trad: %.20e TLTE: %.20e\n",Tgas,Trad,TLTE);
	      printf("ug1:  %.20e ur1:  %.20e\n",pp[VX],pp[FX0]);
	      printf("ug2:  %.20e ur2:  %.20e\n",pp[VY],pp[FY0]);
	      printf("ug3:  %.20e ur3:  %.20e\n\n",pp[VZ],pp[FZ0]);
	      print_NVvector(uu);
	    }
      */
		  
	  //checking if overshooted significantly
	  /*
	  if(((Tgas-Trad)*(Tgas00-Trad00)<0. && fabs(my_max(Tgas-Tgas00,Trad-Trad00)/my_min(Tgas00,Trad00))>1.e-4) ||
	     ((pp[VX]-pp[FX0])*(pp00[VX]-pp00[FX0])<0. && 
	      fabs(my_max(fabs(pp[VX]-pp00[VX]),fabs(pp[FX0]-pp00[FX0]))/my_max(1.e-8,my_min(pp00[VX],pp00[FX0])))>1.e-4) ||
	     ((pp[VY]-pp[FY0])*(pp00[VY]-pp00[FY0])<0. && 
	      fabs(my_max(fabs(pp[VY]-pp00[VY]),fabs(pp[FY0]-pp00[FY0]))/my_max(1.e-8,my_min(pp00[VY],pp00[FY0])))>1.e-4) ||
	     ((pp[VZ]-pp[FZ0])*(pp00[VZ]-pp00[FZ0])<0. && 
	      fabs(my_max(fabs(pp[VZ]-pp00[VZ]),fabs(pp[FZ0]-pp00[FZ0]))/my_max(1.e-8,my_min(pp00[VZ],pp00[FZ0])))>1.e-4)
	     )
	    overshoot=1;
	  */
	  //temperature only

	  /* //order only
	  if(((Tgas-Trad)*(Tgas00-Trad00)<0. && fabs(my_max(Tgas-Tgas00,Trad-Trad00)/my_min(Tgas00,Trad00))>1.e-4))
	    {
	      overshoot=1;
	      //return -1;
	    }
	  */

	  /*
	  //control over initial and LTE values
  ldouble OSEPS=1.e-5;

	  if(Tgas00 > Trad00)
	    if(Tgas>(1.+OSEPS)*Tgas00 || Tgas<(1.-OSEPS)*TLTE || Trad<(1.-OSEPS)*Trad00 || Trad>(1.+OSEPS)*TLTE)
	      overshoot=1;
	  if(Tgas00 < Trad00)
	    if(Trad>(1.+OSEPS)*Trad00 || Trad<(1.-OSEPS)*TLTE || Tgas<(1.-OSEPS)*Tgas00 || Tgas>(1.+OSEPS)*TLTE)
	      overshoot=1;
	  */
	  /*
	  //control over LTE temperature
	  ldouble OSEPS=1.e-6;

	  if(Tgas00 > Trad00)
	    if(Tgas<(1.-OSEPS)*TLTE || Trad>(1.+OSEPS)*TLTE)
	      overshoot=1;
	  if(Tgas00 < Trad00)
	    if(Trad<(1.-OSEPS)*TLTE || Tgas>(1.+OSEPS)*TLTE)
	      overshoot=1;
	  
	  //to dump a case which worked but overshot
	  //if(overshoot==1) return -1;
	  
	  //**************
	  //override
	  overshoot=0;
	  //**************
	  
	  if(overshoot==1)
	    {
	      ldouble Tgasnew=Tgas;
	      ldouble Tradnew=Trad;
	      //making crossing temperature LTE or give it bndr value
	      if(Tgas00 > Trad00)
		{
		  //if(Tgas>(1.+OSEPS)*Tgas00)
		  //Tgasnew=Tgas00;
		  if(Tgas<(1.-OSEPS)*TLTE)
		    Tgasnew=TLTE;
		  //if(Trad<(1.-OSEPS)*Trad00)
		  //Tradnew=Trad00;
		  if(Trad>(1.+OSEPS)*TLTE)
		    Tradnew=TLTE;
		}
	      else
		{
		  //if(Trad>(1.+OSEPS)*Trad00)
		  //Tradnew=Trad00;
		  if(Trad<(1.-OSEPS)*TLTE)
		    Tradnew=TLTE;
		  //if(Tgas<(1.-OSEPS)*Tgas00)
		  //Tgasnew=Trad00;
		  if(Tgas>(1.+OSEPS)*TLTE)
		    Tgasnew=TLTE;
		}
	      
	      if(verbose)
		{
		  printf("overshooted: gas: %.20e -> %.20e -> %.20e\n",Tgas00,Tgas,Tgasnew);
		  printf("overshooted: rad: %.20e -> %.20e -> %.20e\n\n",Trad00,Trad,Tradnew);
		      
		  //getchar();
		}


	      //correct primitives:
	      pp[UU] = calc_PEQ_ufromTrho(Tgasnew,pp[RHO]);
	      pp[EE0] = calc_LTE_EfromT(Tradnew);
	      overshoot=0; //to go directly to the new step
	     
	      
	    
	  */
	      /*
	      //damping the correction approach

	      xi[0]=xi[1]=xi[2]=xi[3]=1.;
	      if((Tgas-Trad)*(Tgas00-Trad00)<0.)
		xi[0] = fabs(Trad00-Tgas00)/(fabs(Trad-Trad00)+fabs(Tgas-Tgas00));
	      for(i1=0;i1<3;i1++)
		if((pp[VX+i1]-pp[FX0+i1])*(pp00[VX+i1]-pp00[FX0+i1])<0.)
		  xi[1+i1] = fabs(pp00[FX0+i1]-pp00[VX+i1])/(fabs(pp[FX0+i1]-pp00[FX0+i1])+fabs(pp[VX+i1]-pp00[VX+i1]));

	      //ldouble minxi = my_min(my_min(xi[0],xi[1]),my_min(xi[2],xi[3]));
	      ldouble minxi=xi[0]; //temperature only

	      ldouble maxxi=0.99;
	      if(minxi>maxxi) //skip if correction small or if fraction already small
		overshoot=0;
	      else
		{
		  //verbose=1;
	      
		  if(verbose)
		    {
		      printf("cnt: %d\n",overcnt);
		      printf("overshooted: gas: %.20e -> %.20e vs rad: %.20e -> %.20e\n",Tgas00,Tgas,Trad00,Trad);
		      printf("overshooted: v1: %.20e -> %.20e vs rad: %.20e -> %.20e\n",pp00[VX],pp[VX],pp00[FX0],pp[FX0]);
		      printf("overshooted: v2: %.20e -> %.20e vs rad: %.20e -> %.20e\n",pp00[VY],pp[VY],pp00[FY0],pp[FY0]);
		      printf("overshooted: v3: %.20e -> %.20e vs rad: %.20e -> %.20e\n",pp00[VZ],pp[VZ],pp00[FZ0],pp[FZ0]);
		      print_4vector(xi);
		      getchar();
		    }
		  
		  minxi/=1.5; //to be generous when damping step

		  //minxi=0.5; //override step
		  xiapp*=minxi;
		}
	      */
      /*
	      overcnt++;

    }
    }
      while(overshoot==1); //overshooting loop

      //if(overcnt>10) printf("over cnt: %d\n",overcnt);	    
      */


      /*
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+sh]-ppp[i+sh]);
	  if(i==0)
	    if(dominates==RAD)
	      f3[i]=fabs(f3[i]/ppp[EE0]);
	    else
	      f3[i]=fabs(f3[i]/ppp[UU]);
	  else
	    f3[i]=fabs(f3[i]/my_max(EPS,fabs(ppp[i+sh])));		  
	}
      //override - convergence with respect to the iterated quantity
      f3[0]=fabs((pp[sh]-ppp[sh])/ppp[sh]);
	  
      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	{
	  if(verbose) printf("success ===\n");
	  break;
	}
      */
	  
      if(iter>MAXITER)
	{
	  if(verbose) 
	    {
	      printf("iter exceeded in solve_implicit_lab_4dprim() for frdt=%f \n",dt);	  
	    }
	  return -1;
	}
    }
  while(1); //main solver loop

  //to print number of iterations
  //if(geom->iy==0 && geom->ix==NX/3) printf("%d\n",iter);

  //returns corrections to radiative primitives
  deltas[0]=uu[EE0]-uu00[EE0];
  deltas[1]=uu[FX0]-uu00[FX0];
  deltas[2]=uu[FY0]-uu00[FY0];
  deltas[3]=uu[FZ0]-uu00[FZ0];

  if(verbose) print_4vector(deltas);

  return 0;
}

//**********************************************************************
//* wrapper ************************************************************
//**********************************************************************

int
solve_implicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  int iv;ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz); //conserved after advection and geometry
      pp[iv]=get_u(p,iv,ix,iy,iz); //some reasonable estimate of primitives 
   }
  
  //inversion to get the right pp[]
  //(u2p checks against proper entropy evolution and uses entropy inversion if necessary

  int corr[2],fixup[2],params[3],ret;
  u2p(uu,pp,&geom,corr,fixup);
  p2u(pp,uu,&geom);

  //1d solver in temperatures only
  ret=solve_implicit_lab_1dprim(uu,pp,&geom,dt,deltas,0,pp);
  //if(ret<0) solve_implicit_lab_1dprim(uu,pp,&geom,dt,deltas,1,pp);
  
  //4d solver starting from the solution satisfying above

  params[1]=RADIMPLICIT_ENERGYEQ;
  params[2]=RADIMPLICIT_LABEQ;

  //deltas[0]=deltas[1]=deltas[2]=deltas[3]=0.;return 0;

  ret=solve_implicit_lab_4dprim(uu,pp,&geom,dt,deltas,verbose,params);

  if(ret<0)
    {
      params[1]=RADIMPLICIT_ENTROPYEQ;
      params[2]=RADIMPLICIT_FFEQ;
      ret=solve_implicit_lab_4dprim(uu,pp,&geom,dt,deltas,verbose,params);
      if(ret<0) return -1;
    }
    
  return 0;

  //return solve_implicit_lab_4dcon(uu,pp,&geom,dt,deltas,verbose);
}


///**********************************************************************
//* test wrapper ************************************************************
//**********************************************************************

int
test_solve_implicit_lab()
{
  FILE *in = fopen("imp.problem.0","r");
  int i1,i2,iv;
  ldouble uu[NV],pp[NV],dt;
  struct geometry geom;

  for (i1=0;i1<NV;i1++)
    iv=fscanf(in,"%lf",&uu[i1]);
  for (i1=0;i1<NV;i1++)
    iv=fscanf(in,"%lf",&pp[i1]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<5;i2++)
      iv=fscanf(in,"%lf ",&geom.gg[i1][i2]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<5;i2++)
      iv=fscanf(in,"%lf ",&geom.GG[i1][i2]);
  iv=fscanf(in,"%lf ",&dt);
  iv=fscanf(in,"%lf ",&geom.alpha);
  iv=fscanf(in,"%lf ",&geom.gdet);
  fclose(in);

  /*
  print_Nvector(uu,NV);
  print_Nvector(pp,NV);
  print_metric(geom.gg);
  print_metric(geom.GG);
  printf("%e %e %e\n",dt,geom.alpha,geom.gdet);
  getchar();
  */
 
  ldouble deltas[4];
  int verbose=1;
 
  int params[3];
  params[1]=RADIMPLICIT_ENERGYEQ;
  params[2]=RADIMPLICIT_LABEQ;
  return solve_implicit_lab_4dprim(uu,pp,&geom,dt,deltas,verbose,params);
}

///**********************************************************************
//* Jon test wrapper ************************************************************
//**********************************************************************

int
test_jon_solve_implicit_lab()
{
  //NOGDET please

  FILE *in = fopen("jon.problem.pre","r");
  int i1,i2,iv,ifile;
  ldouble uu[NV],pp[NV],pp0[NV],dt;
  struct geometry geom;

  for(ifile=1;ifile<=165;ifile++)
    {
      printf("\n            &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
      printf("            &&&&&&&&&&&& case %4d &&&&&&&&&&&&&\n",ifile);
      printf("            &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n\n");

      for (i1=0;i1<NV;i1++)
	iv=fscanf(in,"%lf",&uu[i1]);
      for (i1=0;i1<NV;i1++)
	iv=fscanf(in,"%lf",&pp[i1]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  iv=fscanf(in,"%lf ",&geom.gg[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  iv=fscanf(in,"%lf ",&geom.GG[i1][i2]);
      iv=fscanf(in,"%lf ",&dt);
      iv=fscanf(in,"%lf ",&geom.alpha);
      iv=fscanf(in,"%lf ",&geom.gdet);

      //      uu[EE0]/=1000.;

      //fill missing parts
      ldouble ucon[4];
      ucon[1]=pp[VX];
      ucon[2]=pp[VY];
      ucon[3]=pp[VZ];
      conv_vels(ucon,ucon,VEL4,VEL4,geom.gg,geom.GG);
      geom.alpha=sqrt(-1./geom.GG[0][0]);
      pp[5]=calc_Sfromu(pp[0],pp[1]);

      //destroy magn field
      //uu[B1]=uu[B2]=uu[B3]=pp[B1]=pp[B2]=pp[B3]=0.;

      printf("\n...........................\nJon's input:\n\n");
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);
      //print_metric(geom.gg);
      //print_metric(geom.GG);
      //printf("%e %e %e\n",dt,geom.alpha,geom.gdet);
      //printf("ut: %e\n",ucon[0]);
      int corr[2],fixup[2],u2pret,radcor;
     
      //test
      /*
      u2p_rad(uu,pp,&geom,&radcor);
      printf("radcor: %d\n",radcor);
      print_Nvector(pp,NV);
      p2u_rad(pp,uu,&geom);
      print_Nvector(uu,NV);
      getchar();
      */

      //printf("inverting...\n");
      u2pret=u2p_solver(uu,pp,&geom,U2P_HOT,0); //hd
      if(u2pret<0) printf("u2pret mhd: (%d)\n",u2pret);
      u2p_rad(uu,pp,&geom,&radcor); //rad
      if(radcor!=0) printf("u2pcor rad: (%d)\n",radcor);
      printf("\n..........................\nafter u2p_HOT:\n\n");
      print_Nvector(pp,NV);

      //compare change in entropy
      ucon[1]=pp[VX];
      ucon[2]=pp[VY];
      ucon[3]=pp[VZ];
      conv_vels(ucon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
      ldouble s1=exp(uu[ENTR]/ucon[0]/pp[RHO]);
      ldouble s2=exp(pp[ENTR]/pp[RHO]);

      printf("\n..........................\nchange in entropy:\n\n");
      printf("s(adv) | s(inv): %e | %e\n",s1,s2); 
     
      if(s2/s1 < 0.9 | u2pret<0.)
	{ 
	  printf("\n PROBLEM DETECTED IN ENTROPY OR U2P_HOT DID NOT SUCCEED!\n");
	  u2pret=u2p_solver(uu,pp,&geom,U2P_ENTROPY,0); //hd
	  if(u2pret<0) printf("u2pret mhd: (%d)\n",u2pret);
	  printf("\n..........................\nafter u2p_ENTROPY:\n\n");
	  print_Nvector(pp,NV);
	}      
      
      printf("\n..........................\nafter p2u:\n\n");
      p2u(pp,uu,&geom);
      for (i1=0;i1<NV;i1++)
	pp0[i1]=pp[i1];
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);


      getchar();
   
      ldouble deltas[4];
      int verbose=1;
      int params[3];
      
      solve_explicit_lab_core(uu,pp,&geom,dt,deltas,verbose);
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LABEQ;
      solve_implicit_lab_4dprim(uu,pp,&geom,dt,deltas,verbose,params);
      //solve_implicit_lab_4dcon(uu,pp,&geom,dt,deltas,verbose);

      getchar();
    }

  fclose(in);

  return 0;
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
  ldouble E=pp[EE0];  
  ldouble pr=(GAMMA-1.)*(u);
  ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble xx=get_x(ix,0);
  ldouble yy=get_x(iy,1);
  ldouble zz=get_x(iz,2);
  ldouble kappa=calc_kappa(rho,T,xx,yy,zz);
  ldouble chi=kappa+calc_kappaes(rho,T,xx,yy,zz);  
  ldouble B = SIGMA_RAD*pow(T,4.)/Pi;

  ldouble Fold[3]={pp[FX0],pp[FY0],pp[FZ0]};
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

  deltas[0]=E-pp[EE0];

  return 0;

}


//**********************************************************************
//******* solves explicitly gas - radiation interaction  *******************
//******* in the lab frame, returns vector of deltas **********************
//**********************************************************************
int
solve_explicit_lab_core(ldouble *uu,ldouble *pp,void* ggg,ldouble dt,ldouble* deltas,int verbose)
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

  ldouble Gi[4];
  calc_Gi(pp,geom,Gi);
  indices_21(Gi,Gi,geom->gg);
  
  deltas[0]=-Gi[0]*dt*gdetu;
  deltas[1]=-Gi[1]*dt*gdetu;
  deltas[2]=-Gi[2]*dt*gdetu;
  deltas[3]=-Gi[3]*dt*gdetu;

  if(verbose)
    {
      ldouble delapl[NV];
      ldouble uu0[NV],pp0[NV];
      
      int iv;
      for(iv=0;iv<NV;iv++)
	{
	  delapl[iv]=0.;
	  uu0[iv]=uu[iv];
	  pp0[iv]=pp[iv];
	}
      
      delapl[1]=-deltas[0];
      delapl[2]=-deltas[1];
      delapl[3]=-deltas[2];
      delapl[4]=-deltas[3];
      delapl[EE0]=deltas[0];
      delapl[FX0]=deltas[1];
      delapl[FY0]=deltas[2];
      delapl[FZ0]=deltas[3];

      for(iv=0;iv<NV;iv++)
	{
	  uu0[iv]+=delapl[iv];
	}

      int corr[2],fixup[2];

      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@");
      printf("\n@@@@@@@@ EXPLICIT @@@@@@@@@@@@");
      printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n");
  
      //print_metric(geom->gg);
      //print_metric(geom->GG);
      //printf("%e\n",geom->alpha);
      //printf("%e\n",geom->gdet);
      //print_Nvector(uu0,NV);
      //print_Nvector(pp,NV);
      
      //u2p(uu,pp,geom,corr,fixup);
      //printf("%d %d\n",corr[0],corr[1]);
      //print_Nvector(pp,NV);
      //calc_Gi(pp,geom,Gi);
      
      //print_4vector(Gi);
      //getchar();
      //indices_21(Gi,Gi,geom->gg);

      print_4vector(deltas);
      ldouble T=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO]);
      printf("Tgas: %e\n",T);
    }
  
  return 0;

}


//**********************************************************************
//******* solves explicitly gas - wrapper  *******************
//**********************************************************************
int
solve_explicit_lab(int ix,int iy,int iz,ldouble dt,ldouble* deltas,int verbose)
{
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);
  
  int iv,ret;
  ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=get_u(p,iv,ix,iy,iz); //primitives corresponding to zero-state  
      uu[iv]=get_u(u,iv,ix,iy,iz);  
    }

  ret= solve_explicit_lab_core(uu,pp,&geom,dt,deltas,verbose);

  return ret;
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

  /*
  //as in Ramesh's code
  ldouble Ehat=0.;
  indices_2221(Rij,Rij,gg);
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Ehat+=Rij[i][j]*ucov[i]*ucon[j];

  for(i=0;i<4;i++)
    {
      Gi[i]= - (kappaes*Ehat + kappa*4.*Pi*B)*ucon[i];
      for(j=0;j<4;j++)
	Gi[i]-=(chi)*Rij[i][j]*ucon[j];
    }
  */

#else //Eddington apr. here

  ldouble rho=pp[RHO];
  ldouble u=pp[1];
  ldouble ucov[4],ucon[4]={0,pp[2],pp[3],pp[4]};
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);
  ldouble EE=pp[EE0];
  ldouble Fcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};
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
  ldouble E=pp[EE0];
  ldouble F[3]={pp[FX0],pp[FY0],pp[FZ0]};

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
  Erf=pp[EE0];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
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
  ldouble EE=pp[EE0];
  ldouble Fcon[4]={0.,pp[FX0],pp[FY0],pp[FZ0]};
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
      pp[EE0]=ERADATMMIN; 
      ldouble ucon[4];
      calc_normalobs_4vel(GG,ucon);
      conv_vels(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);
      pp[FX0]=ucon[1]; 
      pp[FY0]=ucon[2];
      pp[FZ0]=ucon[3];

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
     
      pp[FX0]=ucon[1];
      pp[FY0]=ucon[2];
      pp[FZ0]=ucon[3];

      //    print_4vector(ucon); getchar();
      pp[EE0]=ERADATMMIN; 
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
     
      pp[EE0]=ERADATMMIN*(rout/r)*(rout/r)*(rout/r)*(rout/r);

      ldouble ut[4]={0.,-gammamax*pow(r/rout,1.),0.,0.};

      ldouble ggBL[4][5],GGBL[4][5];
      calc_g_arb(xxBL,ggBL,KERRCOORDS);
      calc_G_arb(xxBL,GGBL,KERRCOORDS);

      conv_vels(ut,ut,VELR,VEL4,ggBL,GGBL);

      trans2_coco(xxBL,ut,ut,KERRCOORDS,MYCOORDS);

      conv_vels(ut,ut,VEL4,VELPRIM,gg,GG);
      
      pp[FX0]=ut[1];      
      pp[FY0]=ut[2];      
      pp[FZ0]=ut[3];

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
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
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

  int iv;
  for(iv=0;iv<NV;iv++)
    delapl[iv]=0.;

  delapl[1]=-del4[0];
  delapl[2]=-del4[1];
  delapl[3]=-del4[2];
  delapl[4]=-del4[3];
  delapl[EE0]=del4[0];
  delapl[FX0]=del4[1];
  delapl[FY0]=del4[2];
  delapl[FZ0]=del4[3];

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
  solve_explicit_lab(ix,iy,iz,dt,del4,0);
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
  for(iv=0;iv<NV;iv++)
    delapl[iv]=0.;

  int method=0;
  
  if(method==0) //checks if inversion succesful and then max of du/u < DULIMIT
    {
      //calculating radforce
      //new primitives
      calc_primitives(ix,iy,iz);
      //rad-for-force
      solve_explicit_lab(ix,iy,iz,dt,del4,0);
      //del4[] will be passed up
      indices_21(del4,del4,gg); 
      //changes to conserved
      delapl[1]=-del4[0];
      delapl[2]=-del4[1];
      delapl[3]=-del4[2];
      delapl[4]=-del4[3];
      delapl[EE0]=del4[0];
      delapl[FX0]=del4[1];
      delapl[FY0]=del4[2];
      delapl[FZ0]=del4[3];

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
      solve_explicit_lab(ix,iy,iz,dt,del4,0);
      //del4[] will be passed up
      indices_21(del4,del4,gg); 

      ldouble DULIMIT=0.1;

      //gettin' pp & uu
      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=get_u(u,iv,ix,iy,iz);
	  pp[iv]=get_u(p,iv,ix,iy,iz);
	  delapl[iv]=0.;
	}

      //changes to conserved
      delapl[1]=-del4[0];
      delapl[2]=-del4[1];
      delapl[3]=-del4[2];
      delapl[4]=-del4[3];
      delapl[EE0]=del4[0];
      delapl[FX0]=del4[1];
      delapl[FY0]=del4[2];
      delapl[FZ0]=del4[3];
  
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

  Rij[0][0]=pp1[EE0];
  Rij[0][1]=Rij[1][0]=pp1[FX0];
  Rij[0][2]=Rij[2][0]=pp1[FY0];
  Rij[0][3]=Rij[3][0]=pp1[FZ0];
  Rij[1][1]=Rij[2][2]=Rij[3][3]=1./3.*Rij[0][0];
  Rij[1][2]=Rij[2][1]=Rij[1][3]=Rij[3][1]=Rij[2][3]=Rij[3][2]=0.;

  trans22_on2cc(Rij,Rij,tlo);  
  boost22_ff2lab(Rij,Rij,pp1,gg,GG); 
  indices_2221(Rij,Rij,gg);  

  uufake[EE0]=Rij[0][0];
  uufake[FX0]=Rij[0][1];
  uufake[FY0]=Rij[0][2];
  uufake[FZ0]=Rij[0][3];

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

 
//calculates fluid frame radiative energy density
//and lab-frame four-velocity of gas
int
calc_ff_Rtt(ldouble *pp,ldouble *Rttret, ldouble* ucon,void* ggg)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  //print_Nvector(pp,NV);getchar();
  ucon[0]=0.;
  ucon[1]=pp[VX];
  ucon[2]=pp[VY];
  ucon[3]=pp[VZ];
  ldouble ucov[4];
  conv_vels(ucon,ucon,VELPRIM,VEL4,geom->gg,geom->GG);
  indices_21(ucon,ucov,geom->gg);
  ldouble Rij[4][4],Rtt;
  calc_Rij(pp,ggg,Rij);
 
  indices_2221(Rij,Rij,geom->gg);
  Rtt=0.;
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      Rtt+=-Rij[i1][i2]*ucon[i2]*ucov[i1];

  *Rttret = Rtt;

  return 0;

}

//calculates LTE state
int
calc_LTE_state(ldouble *pp,ldouble *ppLTE,void *ggg)
{
  int i,j;
  for(i=0;i<NV;i++)
    ppLTE[i]=pp[i];

  ldouble Rtt,Eff,ucon[4],ugas;
  
  calc_ff_Rtt(pp,&Rtt,ucon,ggg);
  Eff=-Rtt; //en.density of radiation in the fluid frame
  ugas=pp[UU];
  
  ldouble C = -(Eff + ugas);
  ldouble kt = K_BOLTZ/MU_GAS/M_PROTON;
  ldouble A = 4.*SIGMA_RAD*pow(GAMMAM1/pp[RHO]/kt,4.);
  ldouble Trad=calc_LTE_TfromE(Eff);
  ldouble Tgas=calc_PEQ_Tfromurho(ugas,pp[RHO]);

  ldouble cbrtnaw=cbrt(9.*A + Sqrt(3.)*Sqrt(27.*Power(A,2.) - 256.*Power(A,3.)*Power(C,3.)));
  ldouble ugasLTE=-Sqrt((4*cbrt(2./3.)*C)/
		  cbrtnaw +  cbrtnaw/
		   (cbrt(2.*3.*3.)*A))/2. +
    Sqrt((-4*cbrt(2./3.)*C)/cbrtnaw - cbrtnaw/(cbrt(2.*3.*3.)*A) +
	 2./(A*Sqrt((4*cbrt(2./3.)*C)/cbrtnaw + cbrtnaw/(cbrt(2.*3.*3.)*A))))/2.;

  gsl_complex z0,z1,z2,z3;
  gsl_poly_complex_solve_quartic (0.,0.,1./A,C/A,
				  &z0,&z1,&z2,&z3);

  if(fabs(GSL_IMAG(z0)) < SMALL && GSL_REAL(z0)>0.) ugasLTE=GSL_REAL(z0);
  if(fabs(GSL_IMAG(z1)) < SMALL && GSL_REAL(z1)>0.) ugasLTE=GSL_REAL(z1);
  if(fabs(GSL_IMAG(z2)) < SMALL && GSL_REAL(z2)>0.) ugasLTE=GSL_REAL(z2);
  if(fabs(GSL_IMAG(z3)) < SMALL && GSL_REAL(z3)>0.) ugasLTE=GSL_REAL(z3);
 
  if(isnan(ugasLTE)) 
    {
      gsl_complex z0,z1,z2,z3;
      gsl_poly_complex_solve_quartic (0.,0.,1./A,C/A,
				      &z0,&z1,&z2,&z3);

      printf("%e + %e i\n%e + %e i\n%e + %e i\n%e + %e i\n",
	     GSL_REAL(z0),GSL_IMAG(z0),
	     GSL_REAL(z1),GSL_IMAG(z1),
	     GSL_REAL(z2),GSL_IMAG(z2),
	     GSL_REAL(z3),GSL_IMAG(z3));


      print_Nvector(pp,NV);
      printf("%e %e %e | %e %e | %e %e\n",Eff,ugas,pp[RHO],Tgas,Trad,A,C);
     
      my_err("calc_LTE_state() provided nan\n");
      return -1;
    }

  ldouble EffLTE = -C - ugasLTE;
  ldouble TradLTE=calc_LTE_TfromE(EffLTE);
  ldouble TgasLTE=calc_PEQ_Tfromurho(ugasLTE,pp[RHO]);

  /*
  printf("Trad: %e -> %e\nTgas: %e -> %e\nurad: %e -> %e\nugas: %e -> %e\n\n",
	 Trad,TradLTE,Tgas,TgasLTE,
	 Eff,EffLTE,ugas,ugasLTE)
    ;getchar();
  */  

  //updates only uint gas
  ppLTE[UU]=ugasLTE;

  return 0;

  
}

 
//calculates LTE temperature
ldouble
calc_LTE_temp(ldouble *pp,void *ggg)
{
  int i,j;
  ldouble Rtt,Ehat,ucon[4],ugas0,ugas,rho,Ehat0;
  
  calc_ff_Rtt(pp,&Rtt,ucon,ggg);
  Ehat=-Rtt; //en.density of radiation in the fluid frame
  ugas=pp[UU];
  rho=pp[RHO];
    
  ldouble C = (Ehat + ugas);
  ldouble kt = K_BOLTZ/MU_GAS/M_PROTON;
  ldouble arad = 4.*SIGMA_RAD;
  //ldouble A = 4.*SIGMA_RAD*pow(GAMMAM1/pp[RHO]/kt,4.);

  ldouble Trad=calc_LTE_TfromE(Ehat);
  ldouble Tgas=calc_PEQ_Tfromurho(ugas,pp[RHO]);

  ldouble ccc,TradLTE,TgasLTE;
  if(ugas<Ehat)
    {
      do
	{
	  ugas0=ugas;
	  ugas=(kt*rho/GAMMAM1)*sqrt(sqrt((C-ugas)/arad));
	  ccc=fabs((ugas-ugas0)/ugas0);
	}
      while(ccc>1.e-8);
      Ehat = C - ugas;
    }
  else
    {
      do
	{
	  Ehat0=Ehat;
	  Ehat=arad*pow(GAMMAM1*(C-Ehat)/rho/kt,4.);
	  ccc=fabs((Ehat-Ehat0)/Ehat0);
	}
      while(ccc>1.e-8);
      ugas = C - Ehat;
    }
  
  //TODO: when comparable solve quartic

  TradLTE=calc_LTE_TfromE(Ehat);
  TgasLTE=calc_PEQ_Tfromurho(ugas,rho);
  
  
  /*
    printf("Trad: %e -> %e\nTgas: %e -> %e\n\n",
	 Trad,TradLTE,Tgas,TgasLTE)
    ;getchar();
  */
  

  return TradLTE;  
}

 

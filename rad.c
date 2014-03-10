//KORAL - rad.c
//radiation-related routines

#include "ko.h"

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
//******* calculates total opacity over dx[] *****************************
//*********************************************************************
int
calc_tautot(ldouble *pp, ldouble *xx, ldouble *dx, ldouble *tautot)
{
  ldouble rho=pp[RHO];
  ldouble u=pp[UU];  
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

  ldouble uu[NV],pp[NV],err[5]={0.,0.,0.,0.,0.};
  int corr[2]={0,0},fixup[2]={0,0},u2pret,i1,i2;

  for(i=0;i<NV;i++) pp[i]=ppin[i];
  
  //rho may be inconsistent on input if iterating MHD primitives
  ldouble ucon[4],utcon[4],ucov[4];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  utcon[0]=0.;
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

  //correcting rho for MHD prims
  if(whichprim==MHD)
    {
      ldouble rho = (uu0[RHO])/gdetu/ucon[0];
      pp[RHO]=rho;
    }

  //updating entropy
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

  //total inversion, but only whichprim part matters
  p2u(pp,uu,geom);

  //corresponding change in entropy
  uu[ENTR] = uu0[ENTR] - (uu[EE0]-uu0[EE0]);

  //opposite changes in the other quantities and inversion
  u2pret=0;
  if(whichprim==RAD)
    {
      uu[RHO]=uu0[RHO];
      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);  

      int rettemp=0;
      rettemp=u2p_solver(uu,pp,geom,U2P_HOT,0); 
      if(rettemp<0)
	rettemp=u2p_solver(uu,pp,geom,U2P_ENTROPY,0); 

      if(rettemp<0) 
	  u2pret=-2; //to return error

      else u2pret=0;
    }
  if(whichprim==MHD)
    {
      uu[EE0] = uu0[EE0] - (uu[1]-uu0[1]);
      uu[FX0] = uu0[FX0] - (uu[2]-uu0[2]);
      uu[FY0] = uu0[FY0] - (uu[3]-uu0[3]);
      uu[FZ0] = uu0[FZ0] - (uu[4]-uu0[4]);

      u2pret=u2p_rad(uu,pp,geom,corr);
    }   

  if(corr[0]!=0 || corr[1]!=0) 
    ret=1;
  
  if(u2pret<-1) 
    {
      return -1; //allows for entropy and radiation ceilings but does not update conserved 
    }

#ifdef NCOMPTONIZATION
  ldouble nsource=calc_nsource(pp,ggg);
  //test
  //nsource=0.;
  f[4]=uu[NF0] - uu0[NF0] - dt * gdetu * nsource;
  if(fabs(f[4])>SMALL) err[4]=fabs(f[4])/(fabs(uu[NF0])+fabs(uu0[NF0])+fabs(dt*gdetu*nsource)); else err[4]=0.;
#endif

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

      if(fabs(f[1])>SMALL) err[1]=fabs(f[1])/(fabs(uu[FX0])+fabs(uu0[FX0])+fabs(dt*gdetu*Gi[1])); else err[1]=0.;
      if(fabs(f[2])>SMALL) err[2]=fabs(f[2])/(fabs(uu[FY0])+fabs(uu0[FY0])+fabs(dt*gdetu*Gi[2])); else err[2]=0.;
      if(fabs(f[3])>SMALL) err[3]=fabs(f[3])/(fabs(uu[FZ0])+fabs(uu0[FZ0])+fabs(dt*gdetu*Gi[3])); else err[3]=0.;
    }

  /***** LAB FRAME ENERGY/ENTROPY EQS *****/
  if(whichframe==RADIMPLICIT_LAB)
    {      
      if(whichprim==RAD) //rad-primitives
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0] = uu[EE0] - uu0[EE0] + dt * gdetu * Gi[0];
	      if(fabs(f[0])>SMALL) err[0] = fabs(f[0])/(fabs(uu[EE0]) + fabs(uu0[EE0]) + fabs(dt*gdetu*Gi[0])); else err[0]=0.;

	      ldouble bsq=0.;
	      ldouble bcon[4],bcov[4];
	    }
	  else if(whicheq==RADIMPLICIT_ENTROPYEQ)
	    {
	      f[0] = uu[ENTR] - uu0[ENTR] + dt * gdetu * Gi[0];//but this works on hydro entropy and may fail!
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
  if(whichframe==RADIMPLICIT_FF)
    {

      ldouble uconf[4],Rtt0,Rtt;
      //zero - state 
      calc_ff_Rtt(pp0,&Rtt0,uconf,geom);
      //new state
      calc_ff_Rtt(pp,&Rtt,uconf,geom);
      
      ldouble T=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
      ldouble B = SIGMA_RAD*pow(T,4.)/Pi;
      ldouble Ehat = -Rtt;
      ldouble Ehat0 = -Rtt0;
      ldouble dtau=dt/uconf[0];
      ldouble kappaabs=calc_kappa(pp[RHO],T,geom->xx,geom->yy,geom->zz);

      //fluid frame energy equation:
      if(whichprim==RAD) //rad-primitives
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0]=Ehat - Ehat0 + kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(Ehat) + fabs(Ehat0) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
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
      if(whichprim==MHD) //mhd-
	{
	  if(whicheq==RADIMPLICIT_ENERGYEQ)
	    {
	      f[0]=pp[UU] - pp0[UU] - kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(pp[UU]) + fabs(pp0[UU]) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
	    }
	  else if(whicheq==RADIMPLICIT_ENTROPYEQ)
	    {
	      pp0[ENTR]= calc_Sfromu(pp0[RHO],pp0[UU]);
	      pp[ENTR]= calc_Sfromu(pp[RHO],pp[UU]);
	      f[0]=pp[ENTR] - pp0[ENTR] - kappaabs*(Ehat-4.*Pi*B)*dtau;
	      err[0]=fabs(f[0])/(fabs(pp[ENTR]) + fabs(pp0[ENTR]) + fabs(kappaabs*(Ehat-4.*Pi*B)*dtau));
	    }
	  else
	    my_err("not implemented 2\n");	  
	}
    }

#ifndef NCOMPTONIZATION
  if(!isfinite(f[0]) || !isfinite(f[1]) || !isfinite(f[2]) || !isfinite(f[3]))
    return -1;
  *err0=my_max(my_max(err[0],err[1]),my_max(err[2],err[3]));
#else
  if(!isfinite(f[0]) || !isfinite(f[1]) || !isfinite(f[2]) || !isfinite(f[3]) || !isfinite(f[4]))
    return -1;
  *err0=my_max(my_max(my_max(err[0],err[1]),my_max(err[2],err[3])),err[4]);
#endif
  
  return ret;
} 

int
print_state_implicit_lab_4dprim (int iter, ldouble *x, ldouble *f,ldouble err)
{
#ifndef NCOMPTONIZATION
  printf ("iter = %3d x = % .13e % .13e % .13e % .13e "
	  "e = %.10e f(x) = % .10e % .10e % .10e % .10e\n",
	  iter,
	  //	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
	  x[0],x[1],x[2],x[3],err,f[0],f[1],f[2],f[3]);
#else 
 printf ("iter = %3d x = % .13e % .13e % .13e % .13e % .13e "
	  "e = %.10e f(x) = % .10e % .10e % .10e % .10e % .10e\n",
	  iter,
	  //	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
	  x[0],x[1],x[2],x[3],x[4],err,f[0],f[1],f[2],f[3],f[4]);
#endif
   return 0;
}

/*****************************************/
/******* 4D solver working on primitives */
/*****************************************/
int
solve_implicit_lab_4dprim(ldouble *uu00,ldouble *pp00,void *ggg,ldouble dt,ldouble* deltas,int verbose,int *params,ldouble *pp)
{
  int i1,i2,i3,iv,i,j;
  int mom_over_flag;
  ldouble J[NRADVAR][NRADVAR],iJ[NRADVAR][NRADVAR]; 
  ldouble pp0[NV],ppp[NV],uu[NV],uu0[NV],uup[NV]; 
  ldouble f1[NRADVAR],f2[NRADVAR],f3[NRADVAR],xxx[NRADVAR],xxxbest[NRADVAR],err,errbase,errbest;
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;

  for(iv=0;iv<NV;iv++)
    {
      uu0[iv]=uu00[iv]; //zero state for substepping
      pp0[iv]=pp00[iv]; 
      uu[iv]=uu0[iv]; 
      pp[iv]=pp0[iv];     
    }

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

     
      FILE *out = fopen("imp.problem.dat","w");
      for (i1=0;i1<NV;i1++)
	fprintf(out,"%.20e ",uu0[i1]);
      for (i1=0;i1<NV;i1++)
	fprintf(out,"%.20e ",pp0[i1]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  fprintf(out,"%.20e ",gg[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<5;i2++)
	  fprintf(out,"%.20e ",GG[i1][i2]);
      fprintf(out,"%.20e \n",dt);
      fprintf(out,"%.20e \n",geom->alpha);
      fprintf(out,"%.20e \n",geom->gdet);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<4;i2++)
	  fprintf(out,"%.20e ",geom->tup[i1][i2]);
      for (i1=0;i1<4;i1++)
	for (i2=0;i2<4;i2++)
	  fprintf(out,"%.20e ",geom->tlo[i1][i2]);
      fprintf(out,"\n");
      fclose(out);
      printf("dumped problematic case to imp.problem.dat\n");
    }

  //comparing energy densities
  ldouble ugas00[4],Rtt00,Tgas00,Trad00;
  int dominates;
  calc_ff_Rtt(pp0,&Rtt00,ugas00,geom);
   if(-Rtt00>pp0[UU]) 
    dominates = RAD;
  else
    dominates = MHD;
  Tgas00=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO]);
  Trad00=calc_LTE_TfromE(-Rtt00);

  ldouble TLTE=0.;

  ldouble kappa=calc_kappa(pp0[RHO],Tgas00,0.,0.,0.);
  ldouble chi=kappa+calc_kappaes(pp0[RHO],Tgas00,0.,0.,0.);
  ldouble xi1=kappa*dt*(1.+16.*SIGMA_RAD*pow(Tgas00,4.)/pp0[UU]);
  ldouble xi2=chi*dt*(1.+(-Rtt00)/(pp0[RHO]+GAMMA*pp0[UU]));
  ldouble ucon[4];

  if(verbose) 
    {
      printf("\nparams: %d %d %d %d\n\n",params[0],params[1],params[2],params[3]);

      if(params[1]==RADIMPLICIT_ENERGYEQ)
	printf("energy eq.\n");
      if(params[1]==RADIMPLICIT_ENTROPYEQ)
	printf("entropy eq.\n");
      if(params[1]==RADIMPLICIT_LTEEQ)
	printf("LTE eq. - no longer supported\n");
      if(params[2]==RADIMPLICIT_FF)
	printf("fluid\n");
      if(params[2]==RADIMPLICIT_LAB)
	printf("lab\n");

      printf("\n===========\nxi1: %e xi2: %e\n===========\n\n",xi1,xi2);
      ldouble Rtt;
      calc_ff_Rtt(pp0,&Rtt,ucon,geom);
      printf("Ehat: %e\n\n",-Rtt);

      ldouble qsq=0.;
      int i,j;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=pp0[UU+i]*pp0[UU+j]*geom->gg[i][j];
      ldouble gamma2=1.+qsq;
      printf("gamma gas: %e\n\n",sqrt(gamma2));

      ldouble bcon[4],bcov[4],bsq;
      calc_bcon_prim(pp0,bcon,geom);
      indices_21(bcon,bcov,geom->gg); 
      bsq = dot(bcon,bcov);
      printf("bsq: %e\n\n",bsq);

      qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=pp0[EE0+i]*pp0[EE0+j]*geom->gg[i][j];
      gamma2=1.+qsq;
      printf("gamma rad: %e\n\n",sqrt(gamma2));

      ldouble Gi00[4],Gihat00[4];
      calc_Gi(pp0, geom,Gi00);
       
      indices_21(Gi00,Gi00,geom->gg);

      boost2_lab2ff(Gi00,Gihat00,pp0,geom->gg,geom->GG);
      for(iv=0;iv<4;iv++)
	{
	  Gi00[iv]*=dt*gdetu;
	  Gihat00[iv]*=dt*gdetu;
	}
 
      ldouble Trad=calc_LTE_TfromE(-Rtt);
      ldouble Tgas=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO]);
      ldouble rho=pp0[RHO];
      ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
      ldouble kappaes=calc_kappaes(rho,Tgas,-1.,-1.,-1.);  
      printf("\n===========\nkappa: %e chi: %e\n===========\n",kappa,kappa+kappaes);
      printf("\nxi1: %e xi2: %e\n",xi1,xi2);
      
      printf("gas temp: %e\n",Tgas00);      
      printf("rad temp: %e\n",Trad00); 
      printf("LTE temp: %e\n\n",TLTE) ;  
    }

  //test
  //verbose=1;
  
  /******************************************/
  /******************************************/
  /******************************************/
  //choice of primitives to evolve
  int whichprim,whicheq,do_mom_over;
  if(params[0]==-1) //choose here
    {
      if(-Rtt00<1.e-3*pp0[UU]) //rad sub-dominant
	whichprim=RAD;
      else //gas dominant
	whichprim=MHD; 
      params[0]=whichprim;  
    }
  else
    {
      if(params[0]==RAD)
	whichprim=RAD;
      if(params[0]==MHD)
	whichprim=MHD;
    }
 
  whicheq=params[1];
  do_mom_over =params[3];
  
  if(verbose && whichprim==MHD) printf("Working on MHD\n\n");
  if(verbose && whichprim==RAD) printf("Working on RAD\n\n");

  /******************************************/
  /******************************************/
  /******************************************/

  //check if one can compare gas & rad velocities
  if(VELPRIM!=VELPRIMRAD) 
    my_err("implicit solver assumes VELPRIM == VELPRIMRAD\n");

 
  //4dprim
  ldouble EPS = 1.e-8;
  ldouble CONV = RADIMPCONV;  
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
      printf("=== i: %d %d %d\n\n",ix,iy,iz);
      print_conserved(uu);
      print_primitives(pp);
      print_metric(gg);
    }

  if(verbose) 
    {
      printf("\n===\n Trying imp lab 4d prim with dt : %e \n",dt);
    }

  failed=0;
  iter=0;
  errbest=BIG;

  int np;
#ifdef NCOMPTONIZATION
  np=5;
#else
  np=4;
#endif


  do //main solver loop
    {	 
      iter++;
      
      for(i=0;i<NV;i++)
	{
	  ppp[i]=pp[i];
	}	
      
      for(i=0;i<4;i++)
	{
	  xxx[i]=ppp[i+sh];
	}  
      #ifdef NCOMPTONIZATION
      xxx[4]=ppp[NF0];
      #endif
      
      if(verbose>0)
	{
	  int ret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f1,params,&err);
	  print_state_implicit_lab_4dprim (iter-1,xxx,f1,err); 
	  if(ret<0) printf("f_lab_4dprim ret: %d\n",ret);
	}


      //values at base state
      if(f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f1,params,&err)<0) 
	{
	  if(verbose>0) printf("base state\n");
	  return -1;	  
	}

      errbase=err;	  

      //criterion of convergence on the error
      if(err<CONV)
	{
	  //if(verbose) print_NVvector(pp0);
	  if(verbose) printf("\n === success (error) ===\n");
	  if(verbose==2) getchar();
	  break;
	}

      if(err<errbest)
	{
	  errbest=err;
	  for(j=0;j<np;j++)
	    xxxbest[j]=xxx[j];          
	}
	  

      ldouble del;
      //calculating approximate Jacobian
      for(j=0;j<np;j++)
	{
	  //one-way derivatives
	  //try both signs
	  ldouble sign=-1.;
	  for(;;)
	  {
	    if(j==0)
	      {
		//EEE
		
		//uses EPS of the dominating quantity
		/*
		if(dominates==RAD)
		  del=sign*EPS*ppp[EE0];
		else
		  del=sign*EPS*ppp[UU];	   
		*/
		
		//EPS of the iterated quantity
		del=sign*EPS*ppp[sh]; //minus avoids u2p_mhd errors when working on radiative

		pp[j+sh]=ppp[j+sh]+del;
	      }
	    else if(j==4) //number of photons, only for NCOMPTONIZATION
	      {
		del=sign*EPS*ppp[NF0]; 
		pp[NF0]=ppp[NF0]+del;
	      }	    
	    else //decreasing velocity
	      {
		ldouble veleps = EPS*my_sign(ppp[j+sh])*my_max(1.e-6/sqrt(geom->gg[j][j]),fabs(ppp[j+sh]));
		if(ppp[j+sh]>=0.)
		  del=sign*veleps; 
		else
		  del=-sign*veleps; 
		pp[j+sh]=ppp[j+sh]+del;
	      }	    
	   
	      
	    int fret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f2,params,&err);  
	    
	    if(fret<0) 
	      {
		if(sign>0.) //already switched signs
		  {	      
		    if(verbose) printf("Jac mid-state (%d) both signs lead to trouble\n",j);
		    failed=1;
		    break;
		  }
		else
		  {
		    if(verbose) printf("Jac mid-state (%d) trying the other one\n",j);
		    pp[j+sh]=ppp[j+sh];
		    sign*=-1.;
		    continue;
		  }
	      }
  
	    //Jacobian matrix component
	    for(i=0;i<np;i++)
	      {
		if(j<4)
		  J[i][j]=(f2[i] - f1[i])/(pp[j+sh]-ppp[j+sh]);
		else
		  J[i][j]=(f2[i] - f1[i])/(pp[NF0]-ppp[NF0]);		  
	      }

	    if(j<4) 
	      pp[j+sh]=ppp[j+sh];
	    else
	      pp[NF0]=ppp[NF0];
	    break;
	  }

	  if(failed==1)
	    break;
	}
      
      if(failed==1)
	break;
 
      //inversion
#ifndef NCOMPTONIZATION
      if(inverse_44matrix(J,iJ)<0)
	{
	  failed=1;
	  if(verbose) 
	    {
	      print_tensor(J);
	      print_tensor(iJ);

	      int k,l;
	      ldouble JJ[4][4];
	      for(i=0;i<4;i++)
		{
		  for(j=0;j<4;j++)
		    {
		      JJ[i][j]=0.;

		      for(k=0;k<4;k++)     
			JJ[i][j]+=J[i][k]*iJ[k][j];
		    }
		}

	      print_tensor(JJ);

	      printf("Jacobian inversion failed\n");getchar();
	    }
	  break;
	}
#else
      if(inverse_matrix(&J[0][0],&iJ[0][0],5)<0)
	{
	  failed=1;
	  if(verbose) 
	      printf("Jacobian 5x5 inversion failed\n");getchar();
	  break;
	}  
#endif     
	 
      if(verbose) 
	{
	  //	  if(verbose!=2) exit(0);
	  //getchar();
	}


      //applying corrections
      ldouble xiapp=1.;
      do //check whether energy density positive and the error function returns 0
	{	    
	  //original x
	  for(i=0;i<np;i++)
	    {
	      if(i<4)
		xxx[i]=ppp[i+sh];
	      else
		xxx[i]=ppp[NF0];		
	    }

	  //updating x
	  for(i=0;i<np;i++)
	    {
	      for(j=0;j<np;j++)
		{
		  xxx[i]-=xiapp*iJ[i][j]*f1[j];
		}
	    }

	  if(verbose)
	    {
	      printf("\nsub> trying with xi=%e\n",xiapp);
	      print_Nvector(xxx,np);
	      //print_Nvector(&J[0][0],25);
	    }

	  mom_over_flag=0;
	  //check if momenta overshoot
	  ldouble ALLOWANCE; 
	  if(params[3]==1) ALLOWANCE=100.;
	  if(params[3]==2) ALLOWANCE=1.1;	  
	  if(params[3]==3) ALLOWANCE=1.0001;	  
	  ldouble mommin,mommax,momsep;
	  if(do_mom_over)
	    {
	      for(i=1;i<4;i++)
		{
		  //allowed brackets 
		  //TODO: precalculate
		  mommin=my_min(pp0[EE0+i],pp0[UU+i]);
		  mommax=my_max(pp0[EE0+i],pp0[UU+i]);
		  momsep=fabs(mommin-mommax);

		  mommin=mommin - ALLOWANCE*momsep;
		  mommax=mommax + ALLOWANCE*momsep;

		  if(verbose) printf("mom min/max (%d) %e %e\n",i,mommin,mommax);

		  if(xxx[i]<mommin)
		    {
		      if(verbose) printf("overshoot %d-momentum type 1 (%e). resetting to %e\n",i,xxx[i],mommin);
		      xxx[i]=mommin;
		      mom_over_flag=1;
		    }
		  if(xxx[i]>mommax)
		    {
		      if(verbose) printf("overshoot %d-momentum type 2 (%e). resetting to %e\n",i,xxx[i],mommax);
		      xxx[i]=mommax;
		      mom_over_flag=1;
		    }
		}
	    }
	
	  //update primitives
	  for(i=0;i<np;i++)
	    {
	      if(i<4)
		pp[i+sh]=xxx[i];
	      else
		pp[NF0]=xxx[i];
	    }

	  if(whichprim==MHD)
	    {
	      //correct rho to follow new velocity (only for MHD primitives)
	      ucon[1]=pp[2];ucon[2]=pp[3];ucon[3]=pp[4]; ucon[0]=0.;
	      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);  
	      pp[RHO]=(uu0[RHO])/gdetu/ucon[0];
	    }

	  //updating entropy
	  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

	  //updating the other set of quantities
	  p2u(pp,uu,geom);

	  uu[ENTR] = uu0[ENTR] - (uu[EE0]-uu0[EE0]);
  
	  int u2pret;
	  //opposite changes in the other quantities and inversion
	  if(whichprim==RAD)
	    {
	      uu[RHO]=uu0[RHO];
	      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
	      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
	      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
	      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);

	      int rettemp=0;
	      rettemp=u2p_solver(uu,pp,geom,U2P_HOT,0); 
	      if(rettemp<0)
		rettemp=u2p_solver(uu,pp,geom,U2P_ENTROPY,0); 

	      if(rettemp<0) u2pret=-2; //to return error if even entropy inversion failed
	      else u2pret=0;	      

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
	      //report on ceilings
	      if(corr[0]>0)
		{
		  if(verbose) printf("corr: %d\n",corr[0]);
		  u2pret=-2; //not to allow hitting radiation ceiling in rad when doing iterations
		}
	    }    

	    
	  
	  //check if energy density positive and the inversion worked using U2P_HOT
	  if(xxx[0]>0. && u2pret>=-1) break;

	  
	  //if not decrease the applied fraction
	  if(xxx[0]<=0. && 1)
	    {
	      xiapp*=ppp[sh]/(ppp[sh]+fabs(xxx[0]));
	      xiapp*=sqrt(EPS); //not to land too close to zero, but sometimes prevents from finding proper solution
	    }
	  else //u2p error only
	    {
	      xiapp/=2.; 
	    }

	  if(xiapp<1.e-20) 
	    {
	      if(verbose) printf("damped unsuccesfully in implicit_4dprim\n");
	      failed=1;
	      break;
	    }
	}
      while(1); 


      //TODO:
      //this may fail but necessary to start of a.neq.0 runs      
#ifdef BHDISK_PROBLEMTYPE
      if(1 && failed==0 && global_time<100.)
	{
	  //criterion of convergence on relative change of quantities
	  f3[0]=fabs((pp[sh]-ppp[sh])/ppp[sh]);
	  for(i=1;i<np;i++)
	    {
	      if(i<4)
		{
		  f3[i]=pp[i+sh]-ppp[i+sh];
		  f3[i]=fabs(f3[i]/my_max(EPS,fabs(ppp[i+sh])));	
		}
	      if(i==5)
		{
		  f3[i]=pp[NF0]-ppp[NF0];
		  f3[i]=fabs(f3[i]/my_max(EPS,fabs(ppp[NF0])));	
		}
	    }

	  if(verbose) {
	    //print_4vector(&ppp[sh]);
	    //print_4vector(&pp[sh]);
	    //print_4vector(f3);
	  }
	  
	  ldouble CONVREL=EPS;
	  ldouble CONVRELERR=1.-EPS;

	  if(f3[0]<CONVREL && f3[1]<CONVREL && f3[2]<CONVREL && f3[3]<CONVREL && errbase<CONVRELERR)
	    {
	      if(verbose) printf("\n === success (rel.change) ===\n");
	      break;
	    }     
	}
#endif
      
      
     

      if(iter>MAXITER || failed==1)
	break;
    }
  while(1); //main solver loop

  if(iter>MAXITER || failed==1)
    {
      if(verbose)
	{
	  printf("iter (%d) or failed in solve_implicit_lab_4dprim() for frdt=%f (%e)\n",iter,dt,errbest);	  
	}
      
      /*
      ldouble CONVLOOSE=CONV*1.e0;
      if(errbest<CONVLOOSE)
	{
	  if(verbose) printf("\n === success (looser error) ===\n === coming back to errbest (%e): %e %e %e %e === \n",errbest,xxxbest[0],xxxbest[1],xxxbest[2],xxxbest[3]);
	  //update primitives
	  for(i=0;i<4;i++)
	    pp[i+sh]=xxxbest[i];
	}
      else
      */

      return -1;       
    }

  //updating conserved
  if(whichprim==MHD)
    {
      ucon[1]=pp[2];ucon[2]=pp[3];ucon[3]=pp[4]; ucon[0]=0.;
      conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);  
      pp[RHO]=(uu0[RHO])/gdetu/ucon[0];
    }
  //updating entropy
  pp[ENTR]=calc_Sfromu(pp[RHO],pp[UU]);

  p2u(pp,uu,geom);
  
  int u2pret;
  if(whichprim==RAD)
    {
      uu[RHO]=uu0[RHO];
      uu[1] = uu0[1] - (uu[EE0]-uu0[EE0]);
      uu[2] = uu0[2] - (uu[FX0]-uu0[FX0]);
      uu[3] = uu0[3] - (uu[FY0]-uu0[FY0]);
      uu[4] = uu0[4] - (uu[FZ0]-uu0[FZ0]);
      uu[ENTR] = uu0[ENTR];

      int rettemp=0;
      rettemp=u2p_solver(uu,pp,geom,U2P_HOT,0); 
      if(rettemp<0)
	rettemp=u2p_solver(uu,pp,geom,U2P_ENTROPY,0); 
      
      if(rettemp<0) //to return error if neither entropy or energy inversion succeeded
	u2pret=-2; 
      else 
	u2pret=0;
	      
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

      if(corr[0]>0) //if final solution hit ceiling, discard 
	{
	  if(verbose>1)
	    printf("final solution corrected in rad\n");
	  u2pret=-2;

	  //printf("discarding gammamax solution at %d %d\n",geom->ix,geom->iy);
	}
    }    
  if(u2pret<-1) 
    {
      if(verbose>1)
	printf("final solution rejected\n");
      return -1;
    }
  
  //returns corrections to radiative primitives
  deltas[0]=uu[EE0]-(uu0[EE0]);
  deltas[1]=uu[FX0]-(uu0[FX0]);
  deltas[2]=uu[FY0]-(uu0[FY0]);
  deltas[3]=uu[FZ0]-(uu0[FZ0]);

#ifdef NCOMPTONIZATION
  deltas[4]=uu[NF0]-(uu0[NF0]);
#endif
  
  if(verbose) print_4vector(deltas);
  
  if(verbose)
    {
      ldouble delapl[NV],uu[NV];

      int iv;
      for(iv=0;iv<NV;iv++)
	delapl[iv]=0.;

      delapl[0]=0.;
      delapl[1]=-deltas[0];
      delapl[ENTR]=-deltas[0];
      delapl[2]=-deltas[1];
      delapl[3]=-deltas[2];
      delapl[4]=-deltas[3];
#ifdef MAGNFIELD
      delapl[B1]=0.;
      delapl[B2]=0.;
      delapl[B3]=0.;
#endif
      delapl[EE0]=deltas[0];
      delapl[FX0]=deltas[1];
      delapl[FY0]=deltas[2];
      delapl[FZ0]=deltas[3];

      for(iv=0;iv<NV;iv++)
	{
	  uu[iv]=uu0[iv]+delapl[iv];
	}

      int corr[2],fixup[2];
  
      print_conserved(uu);
      print_primitives(pp);

      calc_ff_Rtt(pp,&Rtt00,ugas00,geom);
      Tgas00=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);
      Trad00=calc_LTE_TfromE(-Rtt00);

      printf("Tgas: %e\n",Tgas00);
      printf("Trad: %e\n",Trad00);
      
      int fret=f_implicit_lab_4dprim(pp,uu0,pp0,dt,geom,f2,params,&err);  

      printf("err %e del0: %e\n",err,deltas[0]);
     }

  //succeeded
  //primitives vector returned to pp[]

  //to calculate average number of succesful iterations

  if(whichprim==RAD && whicheq==RADIMPLICIT_ENERGYEQ)
    {
      //#pragma omp critical //to make the counting precise
      global_int_slot[GLOBALINTSLOT_ITERIMPENERRAD]+=iter;
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_NIMPENERRAD]+=1;
    }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENERGYEQ)
    {
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPENERMHD]+=iter;
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_NIMPENERMHD]+=1;
    }
  if(whichprim==RAD && whicheq==RADIMPLICIT_ENTROPYEQ)
    {
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPENTRRAD]+=iter;
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_NIMPENTRRAD]+=1;
    }
  if(whichprim==MHD && whicheq==RADIMPLICIT_ENTROPYEQ)
    {
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPENTRMHD]+=iter;
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_NIMPENTRMHD]+=1;
    }
  if(whicheq==RADIMPLICIT_LTEEQ)
    {
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_ITERIMPLTE]+=iter;
      //#pragma omp critical
      global_int_slot[GLOBALINTSLOT_NIMPLTE]+=1;
    }

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
  set_cflag(RADFIXUPFLAG,ix,iy,iz,0);

  
  int iv;ldouble pp[NV],uu[NV];
  for(iv=0;iv<NV;iv++)
    {
      uu[iv]=get_u(u,iv,ix,iy,iz); //conserved after advection and geometry
      pp[iv]=get_u(p,iv,ix,iy,iz); //some reasonable estimate of primitives 
   }
  
  //inversion to get the right pp[]
  //(u2p checks against proper entropy evolution and uses entropy inversion if necessary

  int corr[2],fixup[2],params[4],ret;
  u2p(uu,pp,&geom,corr,fixup,0);
  p2u(pp,uu,&geom);

  ldouble pp0[NV],pp00[NV],uu0[NV],uu00[NV];
  PLOOP(iv) 
  {
    pp0[iv]=pp[iv];
    uu0[iv]=uu[iv];
    pp00[iv]=pp[iv];
    uu00[iv]=uu[iv];
  }

  //mostly for LTE solver - can be moved further once debug is done
  ldouble ugas0[4],Rtt0,Tgas0,Trad0,Ehat;
  calc_ff_Rtt(pp0,&Rtt0,ugas0,&geom);
  Ehat=-Rtt0;
  Tgas0=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO]);
  Trad0=calc_LTE_TfromE(Ehat);
  ldouble kappa=calc_kappa(pp0[RHO],Tgas0,geom.xx,geom.yy,geom.zz);
  ldouble chi=kappa+calc_kappaes(pp0[RHO],Tgas0,geom.xx,geom.yy,geom.zz);
  ldouble xi1=kappa*dt*(1.+16.*SIGMA_RAD*pow(Tgas0,4.)/pp0[UU]);
  ldouble xi2=chi*dt*(1.+(-Rtt0)/(pp0[RHO]+GAMMA*pp0[UU]));

  ldouble qsq=0.;
  int i,j;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=pp0[UU+i]*pp0[UU+j]*geom.gg[i][j];
  ldouble gammagas2=1.+qsq;
  qsq=0.;
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      qsq+=pp0[EE0+i]*pp0[EE0+j]*geom.gg[i][j];
  ldouble gammarad2=1.+qsq;

  if(verbose) printf("gammas: %e %e\n\n",sqrt(gammagas2),sqrt(gammarad2));

  //**** 0th ****

  //in pp00[] initial guess for solvers
  params[3]=0; //no overshooting check

  //*********** 1th ************
  PLOOP(iv) 
  { pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
  params[1]=RADIMPLICIT_ENERGYEQ;
  params[2]=RADIMPLICIT_LAB;
  if(Ehat<1.e-2*pp0[UU]) params[0]=RAD; else params[0]=MHD;

  ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp); 
  if(ret!=0)
    { 
      //test
      //ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,1,params,pp); 
      //exit(0);    

      params[2]=RADIMPLICIT_FF;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp);
    }      
  if(ret==0) set_cflag(RADFIXUPFLAG,ix,iy,iz,0);
  
  //*********** 2th ************
  if(ret!=0) {
      PLOOP(iv) 
      {	pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LAB;
      if(params[0]==RAD) params[0]=MHD; else params[0]=RAD;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp);
      if(ret!=0)
	{
	  params[2]=RADIMPLICIT_FF;
	  ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp);
	}      
      if(ret==0) set_cflag(RADFIXUPFLAG,ix,iy,iz,0);
    }
 
  //*********** 3th ************
  if(ret!=0) {
      PLOOP(iv) 
      {	pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[1]=RADIMPLICIT_ENTROPYEQ;
      params[2]=RADIMPLICIT_FF;
      if(Ehat<1.e-2*pp0[UU]) params[0]=RAD; else params[0]=MHD;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp);
      if(ret==0) set_cflag(RADFIXUPFLAG,ix,iy,iz,0);
    }

  //*********** 4th ************
  if(ret!=0) {
      PLOOP(iv) 
      {	pp0[iv]=pp00[iv]; uu0[iv]=uu00[iv]; }
      params[1]=RADIMPLICIT_ENTROPYEQ;
      params[2]=RADIMPLICIT_FF;
      if(params[0]==RAD) params[0]=MHD; else params[0]=RAD;
      ret=solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp);
      if(ret==0) set_cflag(RADFIXUPFLAG,ix,iy,iz,0);
    }


  if(ret!=0)
    {
      //report failure, stop and rerun with verbose
      return -1;

      //****
      //nothing worked - ask for fixup

      //to regard as critical error only if far from starting with gammamax which already unphysical
      if(gammagas2<0.9*GAMMAMAXHD*GAMMAMAXHD && gammarad2<0.9*GAMMAMAXRAD*GAMMAMAXRAD)
	{
	  #ifndef MPI
	  fprintf(fout_fail,"%d > rad implicit > (%4d %4d %4d) (t=%.5e) (otpt=%d) > critical failure no. %d!\n",
		  PROCID,geom.ix,geom.iy,geom.iz,global_time,nfout1,global_int_slot[GLOBALINTSLOT_NCRITFAILURES]);
	  #endif

	  printf("%d > rad implicit > (%4d %4d %4d) (t=%.5e) (otpt=%d) > critical failure no %d!\n",
		 PROCID,geom.ix,geom.iy,geom.iz,global_time,nfout1,global_int_slot[GLOBALINTSLOT_NCRITFAILURES]);

	  //#pragma omp critical
	  global_int_slot[GLOBALINTSLOT_NCRITFAILURES]++;
	  //#pragma omp critical
	  global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES]++;
	  if(global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES]>1.e3)
	    {
	      printf("exceeded # of critical failures (%d) - exiting.\n",
		     global_int_slot[GLOBALINTSLOT_NTOTALCRITFAILURES]);
	      exit(-1);
	    }
	}
      else
	{
	  //#pragma omp critical
	  global_int_slot[GLOBALINTSLOT_NRADFIXUPS]++;
	}
  

      set_cflag(RADFIXUPFLAG,ix,iy,iz,-1);

      return 0; 
    }

  //succeeded!
  //solution given in pp[]

  //calculate deltas here
  p2u(pp,uu,&geom);
  PLOOP(iv)
  {
    set_u(p,iv,ix,iy,iz,pp[iv]);
    set_u(u,iv,ix,iy,iz,uu[iv]);
  }

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
      delapl[ENTR]=-deltas[0];
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
      print_Nvector(uu,NV);
      print_Nvector(pp,NV);

      u2p(uu0,pp0,geom,corr,fixup,0);
      printf("%d %d\n",corr[0],corr[1]);

      print_Nvector(uu0,NV);
      print_Nvector(pp0,NV);

      //print_Nvector(pp,NV);
      //calc_Gi(pp,geom,Gi);
      
      //print_4vector(Gi);
      //getchar();
      //indices_21(Gi,Gi,geom->gg);

      print_4vector(deltas);
      getchar();
      //ldouble T=calc_PEQ_Tfromurho(pp0[UU],pp0[RHO]);
      //printf("Tgas: %e\n",T);
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

  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);

  //the four-velocity of fluid in lab frame
  ldouble ucon[4],utcon[4],ucov[4],vpr[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,gg,GG);

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

  //R^ab u_a u_b = Erad in fluid frame
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

  //test
#ifdef COMPTONIZATION
  ldouble Ehatrad = Ruu;
  ldouble Thatrad;

  #ifdef NCOMPTONIZATION //number of photons conserved

  Thatrad = calc_ncompt_Thatrad(pp,ggg,Ehatrad);
  ldouble ThatradBB=calc_LTE_TfromE(Ehatrad);

  //test
  #ifdef MAXDIFFTRADS
  ldouble maxfac=MAXDIFFTRADS;
  if(Thatrad>maxfac*ThatradBB) 
    Thatrad=maxfac*ThatradBB;
  if(Thatrad<1./maxfac*ThatradBB) 
    Thatrad=1./maxfac*ThatradBB; 
  #endif

  #else //thermal comptonization
  
  Thatrad = calc_LTE_TfromE(Ehatrad);
  
  #endif

  
  for(i=0;i<4;i++)
    Gi[i]+=kappaes * Ehatrad * (4.*K_BOLTZ*(Thatrad - Tgas)/M_ELECTR) * (1. + 4.*K_BOLTZ*Tgas/M_ELECTR) * ucon[i]; 

#endif 

  return 0;
}

/****************************************************/
/***** radiative stress energy tensor ***************/
/***** M1 only here *********************************/
/****************************************************/
int
calc_Rij(ldouble *pp, void* ggg, ldouble Rij[][4])
{
#ifdef RADIATION
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  ldouble Erf;
  int verbose=0;
  int i,j;
  ldouble urfcon[4];

  //covariant formulation

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

  #endif
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
int
calc_Rij_visc(ldouble *pp, void* ggg, ldouble Rvisc[][4], int *derdir)
{
  int i,j;
  
  struct geometry *geom
   = (struct geometry *) ggg;
  
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rvisc[i][j]=0.;	
      }

#if (RADVISCOSITY==SHEARVISCOSITY)
  ldouble shear[4][4];
  ldouble nu;
  calc_rad_shearviscosity(pp,ggg,shear,&nu,derdir);

  //radiation rest frame energy density :
  ldouble Erad;
  Erad=pp[EE0]; 

  //multiply by viscosity to get viscous tensor
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rvisc[i][j]= -2. * nu * Erad * shear[i][j];
      }
#endif //SHEARVISCOSITY

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
  
  if(nlen>=1.)
	{
	  f=1.;
	}
  else
    f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
  
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
      calc_normalobs_relvel(GG,ucon);
      conv_vels(ucon,ucon,VELR,VELPRIMRAD,gg,GG);
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
      conv_vels_ut(ucon,ucon,VEL4,VELPRIMRAD,gg,GG);
     
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

      conv_vels_ut(ut,ut,VEL4,VELPRIM,gg,GG);
      
      pp[FX0]=ut[1];      
      pp[FY0]=ut[2];      
      pp[FZ0]=ut[3];

    }

#ifdef NCOMPTONIZATION
  pp[NF0]=calc_NFfromE(pp[EE0]);
#endif
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

ldouble calc_NFfromT(ldouble T)
{
  return A_RAD*T*T*T/2.70118/K_BOLTZ;
}

ldouble calc_NFfromE(ldouble E)
{
  ldouble temp=calc_LTE_TfromE(E);
  return calc_NFfromT(temp);
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
/******* calculates wavespeeds in the lab frame takin 1/@3 in ***********/
/******* radiative rest frame and boosting it to lab frame **************/
/******* using the HARM algorithm - with taul limiter *******************/
/******* or calculating numerically at cell centers *********************/
/************************************************************************/
int
calc_rad_wavespeeds(ldouble *pp,void *ggg,ldouble tautot[3],ldouble *aval,int verbose)
{
  verbose=0;
  //if(geom->ix==NX/2 && geom->iy==0) verbose = 1;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int i,j;
  
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
  //transforming 1/sqrt(3) from radiation rest frame + limiter based on optical depth
  //the regular, non viscous, approach
  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  ldouble urfcon[4]; 
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);

  //square of radiative wavespeed in radiative rest frame
  ldouble rv2rad = 1./3.;
  ldouble rv2,rv2tau;

  //algorithm from HARM to transform the fluid frame wavespeed into lab frame
  ldouble Acov[4],Acon[4],Bcov[4],Bcon[4],Asq,Bsq,Au,Bu,AB,Au2,Bu2,AuBu,A,B,discr,wspeed2;
  ldouble axl,axr,ayl,ayr,azl,azr;
  axl=axr=ayl=ayr=azl=azr=1.;
   
  int dim;
  for(dim=0;dim<3;dim++)
    {
      //characterisitic limiter based on the optical depth (Sadowski+13a)
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

#ifdef FULLRADWAVESPEEDS
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
  ldouble pp0[NV];

  for(iv=0;iv<NV;iv++)
    {
      delapl[iv]=0.;
      pp0[iv]=get_u(p,iv,ix,iy,iz);
    }

  delapl[RHO]=0.;
  delapl[UU]=-del4[0];
  delapl[VX]=-del4[1];
  delapl[VY]=-del4[2];
  delapl[VZ]=-del4[3];
  delapl[ENTR]=-del4[0];
#ifdef MAGNFIELD
  delapl[B1]=0.;
  delapl[B2]=0.;
  delapl[B3]=0.;
#endif
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
int explicit_rad_source_term(int ix,int iy, int iz,ldouble dt)
{
  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEEXPLICIT);   
  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

  ldouble del4[4],delapl[NV];
  int iv;

  //applied explicitly directly in lab frame
  solve_explicit_lab(ix,iy,iz,dt,del4,0);

  //indices_21(del4,del4,geom.gg);

  apply_rad_source_del4(ix,iy,iz,del4);

  set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 

  return 0;
}


/************************************************************************/
/******* implicit radiative source term in lab frame *********************/
/******* solves numerically in 4D **************************************/
/************************************************************************/
int implicit_lab_rad_source_term(int ix,int iy, int iz,ldouble dt)
{
  
  ldouble del4[NRADVAR],delapl[NV];
  int iv;
  int verbose=1;

  set_cflag(RADSOURCETYPEFLAG,ix,iy,iz,RADSOURCETYPEIMPLICITLAB); 

  if(solve_implicit_lab(ix,iy,iz,dt,del4,0)<0)
    {
      set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,-1);
      //numerical implicit in 4D did not work
      if(verbose) 
	{
	  printf("===\nimp_lab didn't work at %d %d %d (%f %f %f)\n",ix,iy,iz,get_x(ix,0),get_x(iy,1),get_x(iz,1));
	  solve_implicit_lab(ix,iy,iz,dt,del4,2);
	  exit(0);
	  getchar();
	}
      //use the explicit-implicit backup method
    }
  else
    {
      //success in lab frame
      set_cflag(RADSOURCEWORKEDFLAG,ix,iy,iz,0); 
      //apply_rad_source_del4(ix,iy,iz,del4);
    }


  return 0;
}


/***************************************/
/* rad viscosity and shear at cell centers */
/***************************************/
int calc_rad_shearviscosity(ldouble *pp,void* ggg,ldouble shear[][4],ldouble *nuret,int *derdir)
{  
  
#if(RADVISCOSITY==SHEARVISCOSITY) //full shear tensor
  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble mfp,nu,mindx;

  if(geom->ix<-1) 
    {
      printf("rad_shear too far at %d %d\n",geom->ix,geom->iy);
      getchar();
    }
  
  calc_shear_lab(pp,ggg,shear,RAD,derdir);  
  indices_1122(shear,shear,geom->GG);
 
  //calculating the mean free path
  calc_rad_visccoeff(pp,ggg,&nu,&mfp,&mindx);
  
  *nuret=nu;

#else //no rad.viscosity
  int i1,i2;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      shear[i1][i2]=0.;
  *nuret=0.;
#endif
  return 0;  
}

 
//calculates fluid frame radiative energy density (-R^t_t = Ehat)
//and lab-frame four-velocity of gas
int
calc_ff_Rtt(ldouble *pp,ldouble *Rttret, ldouble* ucon,void* ggg)
{
  struct geometry *geom
    = (struct geometry *) ggg;
  //print_Nvector(pp,NV);getchar();
  ldouble ucov[4],utcon[4];
  utcon[0]=0.;
  utcon[1]=pp[VX];
  utcon[2]=pp[VY];
  utcon[3]=pp[VZ];
  conv_vels_both(utcon,ucon,ucov,VELPRIM,VEL4,geom->gg,geom->GG);
  //conv_velscov(utcon,ucov,VELPRIM,VEL4,geom->gg,geom->GG);
  //indices_21(ucon,ucov,geom->gg);
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

//calculates normal observer's radiative energy density
//and returns his four-velocity
int
calc_normal_Rtt(ldouble *pp,ldouble *Rttret, ldouble* ncon,void* ggg)
{
  int i1,i2,i,j;
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble ncov[4];
  ldouble Rij[4][4],Rtt;

  calc_normalobs_4vel(geom->GG,ncon);

  //print_4vector(ncon);getchar();
  indices_21(ncon,ncov,geom->gg);

  calc_Rij(pp,ggg,Rij);
  indices_2221(Rij,Rij,geom->gg);
  Rtt=0.;
  for(i1=0;i1<4;i1++)
    for(i2=0;i2<4;i2++)
      Rtt+=-Rij[i1][i2]*ncon[i2]*ncov[i1];

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
  ldouble ugasLTE=0.;

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

  //updates only uint gas
  ppLTE[UU]=ugasLTE;

  return 0;

  
}

//calculates LTE state and return temperature
int
calc_LTE_state_temp(ldouble *pp,void *ggg)
{
  int i,j;
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
  //troublesome
  
  /*
  ldouble ugasLTE=-Sqrt((4*cbrt(2./3.)*C)/
		  cbrtnaw +  cbrtnaw/
		   (cbrt(2.*3.*3.)*A))/2. +
    Sqrt((-4*cbrt(2./3.)*C)/cbrtnaw - cbrtnaw/(cbrt(2.*3.*3.)*A) +
	 2./(A*Sqrt((4*cbrt(2./3.)*C)/cbrtnaw + cbrtnaw/(cbrt(2.*3.*3.)*A))))/2.;
  */
  
  ldouble ugasLTE=0.;

  gsl_complex z0,z1,z2,z3;
  gsl_poly_complex_solve_quartic (0.,0.,1./A,C/A,
				  &z0,&z1,&z2,&z3);

  if(fabs(GSL_IMAG(z0)) < SMALL && GSL_REAL(z0)>0. && GSL_REAL(z0)<-C)  ugasLTE=GSL_REAL(z0);
  if(fabs(GSL_IMAG(z1)) < SMALL && GSL_REAL(z1)>0. && GSL_REAL(z1)<-C) ugasLTE=GSL_REAL(z1);
  if(fabs(GSL_IMAG(z2)) < SMALL && GSL_REAL(z2)>0. && GSL_REAL(z2)<-C) ugasLTE=GSL_REAL(z2);
  if(fabs(GSL_IMAG(z3)) < SMALL && GSL_REAL(z3)>0. && GSL_REAL(z3)<-C) ugasLTE=GSL_REAL(z3);
 
  if(isnan(ugasLTE) || 1)
    {
      //      gsl_complex z0,z1,z2,z3;
      //      gsl_poly_complex_solve_quartic (0.,0.,1./A,C/A,
      //&z0,&z1,&z2,&z3);

      printf("%e\n%e + %e i\n%e + %e i\n%e + %e i\n%e + %e i\n",
	     ugasLTE,GSL_REAL(z0),GSL_IMAG(z0),
	     GSL_REAL(z1),GSL_IMAG(z1),
	     GSL_REAL(z2),GSL_IMAG(z2),
	     GSL_REAL(z3),GSL_IMAG(z3));


      print_Nvector(pp,NV);
      printf("%e %e %e | %e %e | %e %e\n",Eff,ugas,pp[RHO],Tgas,Trad,A,C);
     
      //_err("calc_LTE_state() provided nan\n");
      //return -1;
    }

  ugasLTE=GSL_REAL(z0);
  
  ldouble EffLTE = -C - ugasLTE;
  ldouble TradLTE=calc_LTE_TfromE(EffLTE);
  ldouble TgasLTE=calc_PEQ_Tfromurho(ugasLTE,pp[RHO]);

  
  printf("Trad: %e -> %e\nTgas: %e -> %e\nurad: %e -> %e\nugas: %e -> %e\n\n",
	 Trad,TradLTE,Tgas,TgasLTE,
	 Eff,EffLTE,ugas,ugasLTE)
    ;getchar();
  

  
  return TgasLTE;

  
}

 
 
//calculates LTE temperature
ldouble
calc_LTE_temp(ldouble *pp,void *ggg,int verbose)
{
  int i,j;
  ldouble Rtt,Ehat,ucon[4],ugas0,ugas,rho,Ehat0;
  

  calc_ff_Rtt(pp,&Rtt,ucon,ggg);
  Ehat=-Rtt; //en.density of radiation in the fluid frame
  ugas=pp[UU];
  rho=pp[RHO];

  ldouble Ehat00,ugas00;
  Ehat00=Ehat;
  ugas00=ugas;
    
  ldouble C = (Ehat + ugas);
  ldouble kt = K_BOLTZ/MU_GAS/M_PROTON;
  ldouble arad = 4.*SIGMA_RAD;
  //ldouble A = 4.*SIGMA_RAD*pow(GAMMAM1/pp[RHO]/kt,4.);

  ldouble Trad=calc_LTE_TfromE(Ehat);
  ldouble Tgas=calc_PEQ_Tfromurho(ugas,pp[RHO]);

  if(verbose)
    {
      printf("%e %e %e\n",C,Trad,Tgas);
      printf("= %e %e\n",Ehat,ugas);
    }
  
  ldouble ccc,TradLTE,TgasLTE;
  int iter=0;

  if(ugas00>Ehat00)
    {
      ugas=ugas00;
      do
	{
	  iter++;
	  ugas0=ugas;
	  ugas=(kt*rho/GAMMAM1)*sqrt(sqrt((C-ugas)/arad));
	  ccc=fabs((ugas-ugas0)/ugas0);
	  if(verbose) printf(">> %e %e\n",ugas,ccc);
	}
      while(ccc>1.e-8 && iter<50);
      Ehat = C - ugas;
    }
  
  
  if(ugas00<=Ehat00)
    {
      Ehat=Ehat00;
      do
	{
	  iter++;
	  Ehat0=Ehat;
	  //	  Ehat=arad*pow(GAMMAM1*(C-Ehat)/rho/kt,4.);
	  Ehat = C - sqrt(sqrt(Ehat/arad))*kt*rho/GAMMAM1;
	  ccc=fabs((Ehat-Ehat0)/Ehat0);
	  if(verbose) printf("> %e %e %e\n",Ehat0,Ehat,ccc);
	}
      while(ccc>1.e-8 && iter<50);
      ugas = C - Ehat;
    }
  

  if(!isfinite(ugas) || iter>=50)
    {
      printf("LTE_temp failed\n"); if(verbose) getchar();
      return -1.;
    }
      

  //TODO: when failed solve quartic

  TradLTE=calc_LTE_TfromE(Ehat);
  TgasLTE=calc_PEQ_Tfromurho(ugas,rho);
  
  
  return TradLTE;  
}



//***************************************
// calculates radiative tensor at faces or centers
// returns total, pure M1, visc ~ R^i_j
//***************************************
int f_flux_prime_rad_total(ldouble *pp, void *ggg,ldouble Rij[][4],ldouble RijM1[][4], ldouble Rijvisc[][4])
{    

#ifdef RADIATION
  int i,j;

  ldouble uu[NV];
  p2u(pp,uu,ggg);

  struct geometry *geom
    = (struct geometry *) ggg;
  ldouble Rvisc1[4][4],Rvisc2[4][4];
  int ix,iy,iz; int iix,iiy,iiz; 
  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  calc_Rij(pp,ggg,RijM1); //regular M1 R^ij
  indices_2221(RijM1,RijM1,gg); //R^i_j
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rij[i][j]=RijM1[i][j];
	Rijvisc[i][j]=0.;
      }

  
#if (RADVISCOSITY==SHEARVISCOSITY)
  //when face and shear viscosity put the face primitives at both cell centers and average the viscous stress tensor
  if(geom->ifacedim>-1)
    //face fluxes
    {      
       
      ix=geom->ix;
      iy=geom->iy;
      iz=geom->iz;
      struct geometry geomcent;
      
      iix=ix;iiy=iy;iiz=iz;

      //left
      if(geom->ifacedim==0)
	iix=ix-1;
      if(geom->ifacedim==1)
	iiy=iy-1;
      if(geom->ifacedim==2)
	iiz=iz-1;
      
      //using the face interpolated primitives
      //calc_Rij_visc(pp,&geomcent,Rvisc1);

      int derdir[3]={0,0,0}; //by default centered derivatives in calc_shear
      //using primitives from the cell centers of neighbours

      //left
      fill_geometry(iix,iiy,iiz,&geomcent);
      derdir[geom->ifacedim]=0; //right derivative
      calc_Rij_visc(&get_u(p,0,iix,iiy,iiz),&geomcent,Rvisc1,derdir);
      indices_2221(Rvisc1,Rvisc1,geomcent.gg); //R^i_j

      //right
      fill_geometry(ix,iy,iz,&geomcent);      
      derdir[geom->ifacedim]=0; //left derivative
      calc_Rij_visc(&get_u(p,0,ix,iy,iz),&geomcent,Rvisc2,derdir);
      indices_2221(Rvisc2,Rvisc2,geomcent.gg); //R^i_j
      
      //adding up to M1 tensor
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  {
	    Rijvisc[i][j]=.5*(Rvisc1[i][j]+Rvisc2[i][j]);
	  }      

    }
  else
    //cell centered fluxes for char. wavespeed evaluation
    {
      int derdir[3]={0,0,0};
      calc_Rij_visc(pp,ggg,Rijvisc,derdir);
      indices_2221(Rijvisc,Rijvisc,gg); //R^i_j
      //adding up to M1 tensor
    }  

  /**********************************/
  //damping the viscous term based on char. viscous velocity
  /**********************************/
#ifdef RADVISCMAXVELDAMP
  //basing on radiative Reynolds number Re = diffusiv flux of conserved quantity / conserved quantity
  ldouble dampfac=1.;
  int idim;
  ldouble vel[3]={0.,0.,0.},maxvel=-1.;
 
  //face centers, using given uu
  for(idim=0;idim<3;idim++)
    for(i=1;i<4;i++)
      {
	if(i==2 && NY==1) continue;
	if(i==3 && NZ==1) continue;
	if(fabs(uu[EE0+i])<1.e-10 * fabs(uu[EE0])) continue;
	vel[idim]=Rijvisc[idim+1][i]/(uu[EE0+i]/gdetu)*sqrt(gg[idim+1][idim+1]);
	//printf("%d %e\n",geom->ix,vel); if(geom->ix==0)getchar();
      }
  if(geom->ifacedim>-1)
    maxvel=fabs(vel[geom->ifacedim]);
  else
    maxvel=sqrt(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]);

  //adjust:
  if(maxvel>MAXRADVISCVEL)
    {
      dampfac=MAXRADVISCVEL/maxvel;
      //printf("damping at %d (%e)\n",geom->ix,maxvel);
    }
  else
    dampfac=1.;

  //todo: choose best prescription
  //dampfac= 1. / (1.+pow(maxvel/MAXRADVISCVEL,2.));
  //dampfac = 1.-step_function(-1.+maxvel/MAXRADVISCVEL,.5);
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rijvisc[i][j]*=dampfac;
      }
  
#endif
  //adding up to Rij
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Rij[i][j]+=Rijvisc[i][j];
      }

#endif
#endif

  return 0;
}

//***************************************
// calculates radiative fluxes at faces or centers
//***************************************
int f_flux_prime_rad( ldouble *pp, int idim, void *ggg,ldouble *ff)
{  
#ifdef RADIATION
  int i,j;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  ldouble Rij[4][4],RijM1[4][4],Rijvisc[4][4];
  f_flux_prime_rad_total(pp,ggg,Rij,RijM1,Rijvisc);

  //fluxes to ff[EE0+]

  ff[EE0]= gdetu*(Rij[idim+1][0]);
      
  ff[FX0]= gdetu*(Rij[idim+1][1]);
      
  ff[FY0]= gdetu*(Rij[idim+1][2]);
      
  ff[FZ0]= gdetu*(Rij[idim+1][3]);

#ifdef NCOMPTONIZATION
  ldouble nphrad=pp[NF0];
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,gg,GG);
  ff[NF0]= gdetu*nphrad*urfcon[idim+1];
#endif

#endif
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates shear tensor sigma_ij in the lab frame at cell centers only!
//hdorrad == MHD -> using gas velocity
//hdorrad == RAD -> using radiative velocity
//derdir[] determines the type of derivative in each dimension (left,right,centered)

int
calc_shear_lab(ldouble *pp0, void* ggg,ldouble S[][4],int hdorrad,int *derdir)
{
  int i,j,k,iv;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;
  
  int ix,iy,iz;
  ix=geom->ix;
  iy=geom->iy;
  iz=geom->iz;
  //let's start with derivatives
  ldouble du[4][4]; //du_i,j
  ldouble du2[4][4]; //du^i,j


  int istart,whichvel;
  if(hdorrad==MHD)
    {
      whichvel=VELPRIM;
      istart=VX;
    }
  else if(hdorrad==RAD)
    {
      whichvel=VELPRIMRAD;
      istart=FX(0);
    }

  //neglecting time derivatives
  /*
  ldouble ucontm1[4],ucovtm1[4],ucontm2[4],ucovtm2[4];
  ucontm1[0]=ucontm2[0]=0.; //time component will be calculated

  ucontm1[1]=get_u(ptm1,istart,ix,iy,iz);
  ucontm1[2]=get_u(ptm1,istart+1,ix,iy,iz);
  ucontm1[3]=get_u(ptm1,istart+2,ix,iy,iz);
  ucontm2[1]=get_u(ptm2,istart,ix,iy,iz);
  ucontm2[2]=get_u(ptm2,istart+1,ix,iy,iz);
  ucontm2[3]=get_u(ptm2,istart+2,ix,iy,iz);

  conv_vels(ucontm1,ucontm1,whichvel,VEL4,gg,GG);
  conv_vels(ucontm2,ucontm2,whichvel,VEL4,gg,GG);

  indices_21(ucontm1,ucovtm1,gg);
  indices_21(ucontm2,ucovtm2,gg);

  for(i=0;i<4;i++)
    {
#ifndef ZEROTIMEINSHEAR
      if(fabs(ttm1-ttm2) < SMALL)
	{
	  du[i][0] = 0.;
	  du2[i][0] = 0.;
	}
      else
	{
	  du[i][0]=(ucovtm1[i]-ucovtm2[i])/(ttm1-ttm2);
	  du2[i][0]=(ucontm1[i]-ucontm2[i])/(ttm1-ttm2);
	}
#else //force d/dt = 0 in shear
      du[i][0] = 0.;
      du2[i][0] = 0.;
#endif
    }
  */

  //instead:
  for(i=0;i<4;i++)
    {
      //force d/dt = 0 in shear
      du[i][0] = 0.;
      du2[i][0] = 0.;
    }

  ldouble ppm1[NV],ppp1[NV],pp[NV];
  ldouble ggm1[4][5],GGm1[4][5];
  ldouble ggp1[4][5],GGp1[4][5];
  ldouble xxvecm1[4],xxvec[4],xxvecp1[4];
  ldouble uconm1[4],uconp1[4],utconm1[4],utconp1[4],utcon[4],ucon[4];
  ldouble ucovm1[4],ucovp1[4],ucov[4];
  ldouble enl,enr;
  int idim;

  //four-velocity at cell basing on pp[]
  //xxvec=geom->xxvec;
  get_xx(ix,iy,iz,xxvec);
  for(iv=0;iv<NV;iv++)
    {
      pp[iv]=pp0[iv];
    }
  utcon[1]=pp[istart];  utcon[2]=pp[istart+1];  utcon[3]=pp[istart+2];
  conv_vels_both(utcon,ucon,ucov,whichvel,VEL4,gg,GG);  
   
  //derivatives
  for(idim=1;idim<4;idim++)
    {
      if(idim==1)
	{
	  get_xx(ix-1,iy,iz,xxvecm1);
	  get_xx(ix+1,iy,iz,xxvecp1);
	  if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix-1,iy,iz);
	      enr=get_u(u,EE0,ix+1,iy,iz);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix-1,iy,iz);
	      enr=get_u(u,UU,ix+1,iy,iz);
	    }
	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix-1,iy,iz);
	      ppp1[iv]=get_u(p,iv,ix+1,iy,iz);
	    }
	  pick_g(ix-1,iy,iz,ggm1);  pick_G(ix-1,iy,iz,GGm1);
	  pick_g(ix+1,iy,iz,ggp1);  pick_G(ix+1,iy,iz,GGp1);
	}

      if(idim==2)
	{
	  get_xx(ix,iy-1,iz,xxvecm1);
	  get_xx(ix,iy+1,iz,xxvecp1);	  
	  if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix,iy-1,iz);
	      enr=get_u(u,EE0,ix,iy+1,iz);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix,iy-1,iz);
	      enr=get_u(u,UU,ix,iy+1,iz);
	    }
	  for(iv=0;iv<NV;iv++)
	    {
	      ppm1[iv]=get_u(p,iv,ix,iy-1,iz);
	      ppp1[iv]=get_u(p,iv,ix,iy+1,iz);
	    }
	  pick_g(ix,iy-1,iz,ggm1);  pick_G(ix,iy-1,iz,GGm1);
	  pick_g(ix,iy+1,iz,ggp1);  pick_G(ix,iy+1,iz,GGp1);
	    
	}

     if(idim==3)
       {
	 get_xx(ix,iy,iz-1,xxvecm1);
	 get_xx(ix,iy,iz+1,xxvecp1);
	 if(hdorrad==RAD) 
	    {
	      enl=get_u(u,EE0,ix,iy,iz-1);
	      enr=get_u(u,EE0,ix,iy,iz+1);
	    }
	  else
	    {
	      enl=get_u(u,UU,ix,iy,iz-1);
	      enr=get_u(u,UU,ix,iy,iz+1);
	    }
	 for(iv=0;iv<NV;iv++)
	   {
	     ppm1[iv]=get_u(p,iv,ix,iy,iz-1);
	     ppp1[iv]=get_u(p,iv,ix,iy,iz+1);
	   }
	 pick_g(ix,iy,iz-1,ggm1);  pick_G(ix,iy,iz-1,GGm1);
	 pick_g(ix,iy,iz+1,ggp1);  pick_G(ix,iy,iz+1,GGp1);
       }

     //calculating four velocity

     utconm1[1]=ppm1[istart];  utconm1[2]=ppm1[istart+1];  utconm1[3]=ppm1[istart+2];
     utconp1[1]=ppp1[istart];  utconp1[2]=ppp1[istart+1];  utconp1[3]=ppp1[istart+2];

     conv_vels_both(utconm1,uconm1,ucovm1,whichvel,VEL4,ggm1,GGm1);
     conv_vels_both(utconp1,uconp1,ucovp1,whichvel,VEL4,ggp1,GGp1);

     ldouble dl,dr,dc;
     ldouble dl2,dr2,dc2;
     for(i=0;i<4;i++)
       {
	 dc=(ucovp1[i]-ucovm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
	 dr=(ucovp1[i]-ucov[i]) / (xxvecp1[idim] - xxvec[idim]);
	 dl=(ucov[i]-ucovm1[i]) / (xxvec[idim] - xxvecm1[idim]);
	 dc2=(uconp1[i]-uconm1[i]) / (xxvecp1[idim] - xxvecm1[idim]);
	 dr2=(uconp1[i]-ucon[i]) / (xxvecp1[idim] - xxvec[idim]);
	 dl2=(ucon[i]-uconm1[i]) / (xxvec[idim] - xxvecm1[idim]);

	 //to avoid corners
	 if((ix<0 && iy==0 && iz==0 && idim!=1) ||
	    (iy<0 && ix==0 && iz==0 && idim!=2) ||
	    (iz<0 && ix==0 && iy==0 && idim!=3))
	   {
	     du[i][idim]=dr;
	     du2[i][idim]=dr2;
	   }
	 else if((ix<0 && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy<0 && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz<0 && ix==NX-1 && iy==NY-1 && idim!=3))
	   {
	     du[i][idim]=dl;
	     du2[i][idim]=dl2;
	   }
	 else if((ix>=NX && iy==0 && iz==0 && idim!=1) ||
	    (iy>=NY && ix==0 && iz==0 && idim!=2) ||
	    (iz>=NZ && ix==0 && iy==0 && idim!=3))
	   {
	     du[i][idim]=dr;
	     du2[i][idim]=dr2;
	   }
	 else if((ix>=NX && iy==NY-1 && iz==NZ-1 && idim!=1) ||
	    (iy>=NY && ix==NX-1 && iz==NZ-1 && idim!=2) ||
	    (iz>=NZ && ix==NX-1 && iy==NY-1 && idim!=3))
	   {
	     du[i][idim]=dl;
	     du2[i][idim]=dl2;
	   }
	 else
	   {
	     //choice of 1st order derivative
	     #ifdef NUMRADWAVESPEEDS
	     if(fabs(enl)>fabs(enr))
	       {
		 du[i][idim]=dl;
		 du2[i][idim]=dl2;
	       }
	     else
	       {
		 du[i][idim]=dr;
		 du2[i][idim]=dr2;
	       }
	     #else

	     if(derdir[idim-1]==0)
	       {
		 du[i][idim]=dc;
		 du2[i][idim]=dc2;
	       }
	     if(derdir[idim-1]==1)
	       {
		 du[i][idim]=dl;
		 du2[i][idim]=dl2;
	       }
	     if(derdir[idim-1]==2)
	       {
		 du[i][idim]=dr;
		 du2[i][idim]=dr2;
	       }
	     #endif
	   }

	 if(isnan(du[i][idim])) {
	   printf("nan in shear_lab : %d %d %d %d\n",ix,iy,iz,idim);
	   print_4vector(ucovm1);
	   print_4vector(ucov);
	   print_4vector(ucovp1);
	   print_4vector(xxvecm1);
	   print_4vector(xxvec);	   
	   print_4vector(xxvecp1);
	   getchar();
	 }	 	 
       }       
    }

  //covariant derivative tensor du_i;j
  ldouble dcu[4][4];
  //covariant derivative tensor du^i;j - only for expansion
  ldouble dcu2[4][4];
  ldouble Krsum;

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Krsum=0.;
	  for(k=0;k<4;k++)
	    Krsum+=get_gKr(k,i,j,ix,iy,iz)*ucov[k];

	  dcu[i][j] = du[i][j] - Krsum;
	}
    
      //only diagonal terms for expansion
      Krsum=0.;
      for(k=0;k<4;k++)
	Krsum+=get_gKr(i,i,k,ix,iy,iz)*ucon[k];
    
      dcu2[i][i] = du2[i][i] + Krsum; 
    }

  //expansion
  ldouble theta=0.;
  for(i=0;i<4;i++)
    theta+=dcu2[i][i];

  //projection tensors P11=P_ij, P21=P^i_j
  ldouble P11[4][4],P21[4][4];
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	P11[i][j] = gg[i][j] + ucov[i]*ucov[j];
	P21[i][j] = delta(i,j) + ucon[i]*ucov[j];
      }

  //the shear tensor sigma_ij - only spatial components
  for(i=1;i<4;i++)
    for(j=1;j<4;j++)
      {
	ldouble sum1,sum2;
	sum1=sum2=0.;
	for(k=0;k<4;k++)
	  {
	    sum1+=dcu[i][k]*P21[k][j];
	    sum2+=dcu[j][k]*P21[k][i];
	  }
	S[i][j] = 0.5*(sum1+sum2) - 1./3.*theta*P11[i][j];

     }

  //filling the time component from u^mu sigma_munu = 0 
  //(zero time derivatives in the comoving frame - no need for separate comoving routine)
  for(i=1;i<4;i++)
    S[i][0]=S[0][i]=-1./ucon[0]*(ucon[1]*S[1][i]+ucon[2]*S[2][i]+ucon[3]*S[3][i]);

  S[0][0]=-1./ucon[0]*(ucon[1]*S[1][0]+ucon[2]*S[2][0]+ucon[3]*S[3][0]);


  return 0;
}
        
//calculates mean free path used in the viscosity coefficient
int 
calc_rad_visccoeff(ldouble *pp,void *ggg,ldouble *nuret,ldouble *mfpret,ldouble *mindxret)
{
#if (RADVISCOSITY==SHEARVISCOSITY)
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble chi=calc_chi(pp,geom->xxvec);
  ldouble mfp = 1./chi; // dr / dtau
  ldouble mindx,nu;

  //limiting in opt.thin region

  //minimal cell size
  ldouble dx[3]={get_size_x(geom->ix,0)*sqrt(geom->gg[1][1]),   
		 get_size_x(geom->iy,1)*sqrt(geom->gg[2][2]),
		 get_size_x(geom->iz,2)*sqrt(geom->gg[3][3])};

  if(NY==1 && NZ==1) mindx = dx[0];
  else if(NZ==1) mindx = my_min(dx[0],dx[1]);
  else if(NY==1) mindx = my_min(dx[0],dx[2]);
  else mindx = my_min(dx[0],my_min(dx[1],dx[2]));

//choice of mean free path limiter for tau<<1
#ifdef RADVISCMFPCONST

  if(mfp>RADVISCMFPCONST || chi<SMALL) mfp=RADVISCMFPCONST;

#elif defined(RADVISCMFPCYL)

  ldouble xxBL[4];
  coco_N(geom->xxvec,xxBL,MYCOORDS, BLCOORDS);
  if(mfp>mindx || chi<SMALL) mfp=xxBL[1]*sin(xxBL[2]); //Rcyl = Rsph * sin(th)
  if(mfp<0. || !isfinite(mfp)) mfp=0.;

#elif defined(RADVISCMFPSPH)

  ldouble xxBL[4];
  coco_N(geom->xxvec,xxBL,MYCOORDS, BLCOORDS);
  ldouble mfplim=xxBL[1];
  if(mfp>mfplim || chi<SMALL) mfp=mfplim; //Rcyl = Rsph
  if(mfp<0. || !isfinite(mfp)) mfp=0.;

#else

  if(mfp>mindx || chi<SMALL) mfp=mindx;

#endif

  nu = ALPHARADVISC * mfp;

  /**********************************/
  //damping the viscous term based on maximal diffusion coefficient for given dt ans cell size
  /**********************************/

#ifdef RADVISCNUDAMP
  ldouble nulimit = mindx*mindx / 2. / dt;
  ldouble fac=nu/nulimit;
  if(nu>nulimit)
    {
      //printf("famping at %d %d | %e | %e %e | %e\n",geom->ix,geom->iy,nu,mfp,mindx,dt); getchar();
      nu=nulimit;
    }
#endif


  *nuret=nu;
  *mfpret=mfp;
  *mindxret=mindx;
#endif
  return 0;
}

//**********************************************************************
//* calculates the source term for number of photons
//**********************************************************************
ldouble
calc_nsource(ldouble *pp, void* ggg)
{
  int i,j;
      
#ifdef NCOMPTONIZATION
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet; gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  //velocities of the frames
  ldouble ut[4];ut[1]=pp[VX];ut[2]=pp[VY];ut[3]=pp[VZ];
  ldouble uffcov[4],uffcon[4];
  conv_vels_both(ut,uffcon,uffcov,VELPRIM,VEL4,gg,GG);
  ldouble urfcov[4],urfcon[4];
  ut[1]=pp[FX0];ut[2]=pp[FY0];ut[3]=pp[FZ0];
  conv_vels_both(ut,urfcon,urfcov,VELPRIMRAD,VEL4,gg,GG);
  //radiative stress tensor in the lab frame
  ldouble Rij[4][4];
  calc_Rij(pp,ggg,Rij);

  //R^ab u_a u_b = Erad in fluid frame
  ldouble Ehatrad = 0.;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Ehatrad+=Rij[i][j]*uffcov[i]*uffcov[j];

  //radiation temperature
  ldouble Thatrad = calc_ncompt_Thatrad_4vel(pp,ggg,Ehatrad,urfcon,uffcov);
  //test
  #ifdef MAXDIFFTRADS
  ldouble ThatradBB=calc_LTE_TfromE(Ehatrad);
  ldouble maxfac=MAXDIFFTRADS;
  if(Thatrad>maxfac*ThatradBB) 
    Thatrad=maxfac*ThatradBB;
  if(Thatrad<1./maxfac*ThatradBB) 
    Thatrad=1./maxfac*ThatradBB; 
  #endif

  //gas properties
  ldouble rho=pp[RHO];
  ldouble u=pp[1];
  ldouble p= (GAMMA-1.)*(ldouble)u;
  ldouble Tgas=p*MU_GAS*M_PROTON/K_BOLTZ/rho;
  ldouble B = SIGMA_RAD*pow(Tgas,4.)/Pi;
  ldouble kappa=calc_kappa(rho,Tgas,-1.,-1.,-1.);
  
  //change in number of photons in fluid frame
  ldouble Cn=1./2.70118/K_BOLTZ;
  ldouble ndotff = -kappa*Cn*(Ehatrad/Thatrad - 4.*M_PI*B/Tgas);

  //relative gamma rad-fluid rest frames
  //ldouble relgamma = urfcon[0]*uffcov[0] + urfcon[1]*uffcov[1] +urfcon[2]*uffcov[2] +urfcon[3]*uffcov[3]; 
  //ldouble ndotrf = -ndotff/relgamma; 

  //the rate of change of number of photons invariant
  ldouble ndotrf = ndotff;

  return ndotrf; 
#else
  return 0;
#endif
}

//**********************************************************************
//* calculates the radiation temperature using
//* radiative energy density and the number of photons
//* assumes four-velocities known
//**********************************************************************
ldouble
calc_ncompt_Thatrad_4vel(ldouble *pp, void* ggg, ldouble Ehatrad, ldouble* urfcon, ldouble* uffcov)
{
  //number of photons in rad rest frame
  ldouble nphrf = pp[NF(0)];

  //relative gamma rad-fluid rest frames
  ldouble relgamma = urfcon[0]*uffcov[0] + urfcon[1]*uffcov[1] +urfcon[2]*uffcov[2] +urfcon[3]*uffcov[3]; 

  //ff quantities
  ldouble nphhat = - relgamma*nphrf;
  ldouble Thatrad = 1./2.70118/K_BOLTZ * Ehatrad / nphhat;

  return Thatrad;

}


//**********************************************************************
//* calculates the radiation temperature using
//* radiative energy density and the number of photons
//**********************************************************************
ldouble
calc_ncompt_Thatrad(ldouble *pp, void* ggg, ldouble Ehatrad)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;
  
  //the four-velocity of fluid in lab frame
  ldouble uffcon[4],utcon[4],uffcov[4],vpr[3];
  utcon[1]=pp[2];
  utcon[2]=pp[3];
  utcon[3]=pp[4];
  conv_vels_both(utcon,uffcon,uffcov,VELPRIM,VEL4,gg,GG);

  //the four-velocity of radiation in lab frame
  ldouble urfcon[4];
  urfcon[0]=0.;
  urfcon[1]=pp[FX0];
  urfcon[2]=pp[FY0];
  urfcon[3]=pp[FZ0];
  conv_vels(urfcon,urfcon,VELPRIMRAD,VEL4,geom->gg,geom->GG);

  return calc_ncompt_Thatrad_4vel(pp,ggg,Ehatrad,urfcon,uffcov);

}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//************** tests *************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************

int
test_solve_implicit_lab()
{
  FILE *in = fopen("imp.problem.0","r");

  int i1,i2,iv;
  ldouble uu0[NV],pp0[NV],pp[NV],dt;
  struct geometry geom;

  for (i1=0;i1<NV;i1++)
    iv=fscanf(in,"%lf",&uu0[i1]);
  for (i1=0;i1<NV;i1++)
    iv=fscanf(in,"%lf",&pp0[i1]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<5;i2++)
      iv=fscanf(in,"%lf ",&geom.gg[i1][i2]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<5;i2++)
      iv=fscanf(in,"%lf ",&geom.GG[i1][i2]);
  iv=fscanf(in,"%lf ",&dt);
  iv=fscanf(in,"%lf ",&geom.alpha);
  iv=fscanf(in,"%lf ",&geom.gdet);
  //for imp.problems > 21
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
      iv=fscanf(in,"%lf ",&geom.tup[i1][i2]);
  for (i1=0;i1<4;i1++)
    for (i2=0;i2<4;i2++)
      iv=fscanf(in,"%lf ",&geom.tlo[i1][i2]);

  fclose(in);

  geom.ix=geom.iy=geom.iz=0;

  ldouble Rttcov[4]={uu0[EE0],uu0[FX0],uu0[FY0],uu0[FZ0]};
  ldouble Rttcon[4];
  indices_12(Rttcov,Rttcon,geom.GG);
  //print_4vector(Rttcov);
  //print_4vector(Rttcon);
  ldouble vcon[4],ucon[4],ucov[4];

  //converting to 4-velocity
  vcon[1]=pp0[2];
  vcon[2]=pp0[3];
  vcon[3]=pp0[4];
  vcon[0]=0.;  
  conv_vels_both(vcon,ucon,ucov,VELPRIM,VEL4,geom.gg,geom.GG);
  //conv_velscov(vcon,ucov,VELPRIM,VEL4,geom.gg,geom.GG);
  //indices_21(ucon,ucov,geom.gg);
  print_4vector(ucon);
  print_4vector(ucov);

   //converting to 4-velocity
  /*
  vcon[1]=pp0[FX0];
  vcon[2]=pp0[FY0];
  vcon[3]=pp0[FZ0];
  vcon[0]=0.;  
  conv_vels(vcon,ucon,VELPRIM,VEL4,geom.gg,geom.GG);
  indices_21(ucon,ucov,geom.gg);
  print_4vector(ucon);
  print_4vector(ucov);
  */
  /*
  ldouble pp2[NV];
  int i;
  PLOOP(i) pp2[i]=pp0[i];
  u2p_solver(uu0,pp2,&geom,U2P_HOT,2);
  */

  //print_metric(geom.gg);
  //print_metric(geom.GG); 
  //print_Nvector(uu0,NV);

  //print_Nvector(pp0,NVMHD);
  //p2u(pp0,uu0,&geom);
  //getchar();
  //print_Nvector(uu0,NVMHD);

  /*
  print_metric(geom.gg);
  print_metric(geom.GG);
  printf("%e %e %e\n",dt,geom.alpha,geom.gdet);
  getchar();
  */
 
  /*
  pp0[FX0]/=10.;
  pp0[FY0]/=10.;
  pp0[FZ0]/=10.;
  */

  p2u(pp0,uu0,&geom);

  ldouble deltas[4];
  int verbose=2;
  int params[4];
  
  //return solve_explicit_lab_core(uu,pp,&geom,dt,deltas,verbose);
  
  //solve_implicit_lab_1dprim(uu,pp,&geom,dt,deltas,verbose,pp);
  //  return solve_implicit_lab_4dcon(uu0,pp0,&geom,dt,deltas,verbose,pp);

  

 
  params[0]=MHD;
  params[1]=RADIMPLICIT_ENERGYEQ;
  params[2]=RADIMPLICIT_LAB;
  params[3]=0; 
  solve_implicit_lab_4dprim(uu0,pp0,&geom,dt,deltas,verbose,params,pp);

  exit(1);
  
  
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
      int params[4];
      
      solve_explicit_lab_core(uu,pp,&geom,dt,deltas,verbose);
      params[1]=RADIMPLICIT_ENERGYEQ;
      params[2]=RADIMPLICIT_LAB;
      params[3]=1;
      ldouble ppret[NV];
      solve_implicit_lab_4dprim(uu,pp,&geom,dt,deltas,verbose,params,ppret);
      //solve_implicit_lab_4dcon(uu,pp,&geom,dt,deltas,verbose);

      getchar();
    }

  fclose(in);

  return 0;
}


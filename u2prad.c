//KORAL - u2p.c
//radiative routines related to u2p inversion

#include "ko.h"



// get's gamma assuming fixed E rather than using original R^t_t that we assume is flawed near floor regions.  We want to preserve R^t_i (i.e conserved momentum)
static int get_m1closure_gammarel2_cold(int verbose,
					void *ggg, 
					FTYPE *Avcon, 
					FTYPE *Avcov, 
					FTYPE *gammarel2return, 
					FTYPE *deltareturn, 
					FTYPE *numeratorreturn, 
					FTYPE *divisorreturn, 
					FTYPE *Erfreturn, 
					FTYPE *urfconrel)
{
   struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  FTYPE gamma2,gammarel2,delta;
  FTYPE Erf;
  FTYPE alpha=geom->alpha;
  int jj;

  static FTYPE gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=GG[0][0];
  gn12=GG[0][1];
  gn13=GG[0][2];
  gn14=GG[0][3];
  gn22=GG[1][1];
  gn23=GG[1][2];
  gn24=GG[1][3];
  gn33=GG[2][2];
  gn34=GG[2][3];
  gn44=GG[3][3]; 

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  // choose gamma
  if(gammarel2return==NULL){
    FTYPE gammamax=GAMMAMAXRAD;
    gammarel2=gammamax*gammamax;
  }
  else gammarel2=*gammarel2return; // feed in desired gammarel2

  FTYPE utsq=gammarel2/(alpha*alpha);

  FTYPE Avcovorig[NDIM],Avconorig[NDIM];
  DLOOPA(jj) Avcovorig[jj]=Avcov[jj];
  DLOOPA(jj) Avconorig[jj]=Avcon[jj];

  // get new Avcov[0]=R^t_t

  // NOTEMARK: Note that Sqrt() is only ever negative when gammarel2<0, so never has to be concern.
  Avcov[0]=(0.25*(-4.*(gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*utsq*(gn11 + utsq) + 
                        Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                              1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                       gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2))))/(gn11*utsq*(gn11 + utsq));

  Erf=(0.75*Sqrt((Power(gn12,2)*Power(Rdtx,2) + 2.*gn12*Rdtx*(gn13*Rdty + gn14*Rdtz) + Power(gn13*Rdty + gn14*Rdtz,2) - 
                           1.*gn11*(gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 2.*gn34*Rdty*Rdtz + 
                                    gn44*Power(Rdtz,2)))*utsq*(gn11 + utsq)*Power(gn11 + 4.*utsq,2)))/(utsq*(gn11 + utsq)*(gn11 + 4.*utsq));

  if(verbose) printf("NOR SOL: Avcov0new=%g Avcov0old=%g Erf=%g :: %g %g %g\n",Avcov[0],Avcovorig[0],Erf,Rtx,Rty,Rtz);

  //modify Avcon
  indices_12(Avcov,Avcon,GG);
  if(verbose) DLOOPA(jj) printf("jj=%d Avconorig=%g Avconnew=%g\n",jj,Avconorig[jj],Avcon[jj]);

 
  delta=0; // not yet

  *gammarel2return=gammarel2;
  *deltareturn=delta;

  // get good relative velocity
  FTYPE gammarel=sqrt(gammarel2);

  // get relative 4-velocity
  if(Erf>0.0) SLOOPA(jj) urfconrel[jj] = alpha * (Avcon[jj] + 1./3.*Erf*GG[0][jj]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
  else SLOOPA(jj) urfconrel[jj] = 0.0;

  *Erfreturn=Erf; // pass back new Erf to pointer

  return(0);
}




// get's gamma^2 for lab-frame gamma  using Rd and gcon
static int get_m1closure_gammarel2(int verbose,void *ggg, ldouble *Avcon, ldouble *Avcov, ldouble *gammarel2return,ldouble *deltareturn, ldouble *numeratorreturn, ldouble *divisorreturn)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble gamma2,gammarel2,delta,numerator,divisor;
  ldouble gamma2a,gamma2b;

  // mathematica solution that avoids catastrophic cancellation when Rtt very small (otherwise above gives gamma2=1/2 oddly when gamma2=1) -- otherwise same as above
  // well, then had problems for R~1E-14 for some reason when near BH.  Couldn't quickly figure out, so use no replacement of gv11.
  // see u2p_inversion.nb
  static ldouble gctt, gn11, gn12,  gn13,  gn14,  gn22,  gn23,  gn24,  gn33,  gn34,  gn44,  Rtt,  Rtx,  Rty,  Rtz,  Rdtt,  Rdtx,  Rdty,  Rdtz;
  gn11=GG[0][0];
  gn12=GG[0][1];
  gn13=GG[0][2];
  gn14=GG[0][3];
  gn22=GG[1][1];
  gn23=GG[1][2];
  gn24=GG[1][3];
  gn33=GG[2][2];
  gn34=GG[2][3];
  gn44=GG[3][3];

  Rtt=Avcon[0];
  Rtx=Avcon[1];
  Rty=Avcon[2];
  Rtz=Avcon[3];

  Rdtt=Avcov[0];
  Rdtx=Avcov[1];
  Rdty=Avcov[2];
  Rdtz=Avcov[3];

  gamma2a=(-0.25*(2.*Power(gn11,2)*Power(Rdtt,2) + (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*
        (gn12*Rdtx + gn13*Rdty + gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
            gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
               8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))) + 
       gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 2.*gn24*Rdtx*Rdtz + 
          2.*gn34*Rdty*Rdtz + gn44*Power(Rdtz,2) + Rdtt*
           (4.*gn13*Rdty + 4.*gn14*Rdtz + Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
               gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                  8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2)))))))/
   (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
    2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));


  if( gamma2a<GAMMASMALLLIMIT || isinf(gamma2a) ){
    gamma2b=(0.25*(-2.*Power(gn11,2)*Power(Rdtt,2) - 1.*gn11*(4.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 
                                                              Rdty*(4.*gn13*Rdtt + 2.*gn23*Rdtx + gn33*Rdty) + 2.*(2.*gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2)) + 
                   gn11*Rdtt*Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                  gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                        8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))) + 
                   (gn12*Rdtx + gn13*Rdty + gn14*Rdtz)*(-1.*gn12*Rdtx - 1.*gn13*Rdty - 1.*gn14*Rdtz + 
                                                        Sqrt(4.*Power(gn11,2)*Power(Rdtt,2) + Power(gn12*Rdtx + gn13*Rdty + gn14*Rdtz,2) + 
                                                             gn11*(8.*gn12*Rdtt*Rdtx + 3.*gn22*Power(Rdtx,2) + 8.*gn13*Rdtt*Rdty + 6.*gn23*Rdtx*Rdty + 3.*gn33*Power(Rdty,2) + 
                                                                   8.*gn14*Rdtt*Rdtz + 6.*gn24*Rdtx*Rdtz + 6.*gn34*Rdty*Rdtz + 3.*gn44*Power(Rdtz,2))))))/
      (gn11*Power(Rdtt,2) + 2.*gn12*Rdtt*Rdtx + gn22*Power(Rdtx,2) + 2.*gn13*Rdtt*Rdty + 2.*gn23*Rdtx*Rdty + gn33*Power(Rdty,2) + 
       2.*(gn14*Rdtt + gn24*Rdtx + gn34*Rdty)*Rdtz + gn44*Power(Rdtz,2));
    gamma2=gamma2b;
  }
  else{
    // choose
    gamma2=gamma2a;
  }

  //  if(isnan(gamma2a) && !isnan(gamma2b))
  //    gamma2=gamma2b;

  //  if(!isnan(gamma2a) && isnan(gamma2b))
  //    gamma2=gamma2a;

  ////////////////////////
  //
  //cap on u^t
  //
  ///////////////////////
  ldouble alpha=geom->alpha;


  // get relative 4-velocity, that is always >=1 even in GR
  gammarel2 = gamma2*alpha*alpha;

  //to consider it a failure - leads to instability
  //if(isnan(gammarel2)) return -1;
  
  // check for machine error away from 1.0 that happens sometimes
  if(isnan(gammarel2) || (gammarel2>GAMMASMALLLIMIT && gammarel2<1.0))
    //if((gammarel2>GAMMASMALLLIMIT && gammarel2<1.0))
  {
    if(verbose) printf("Hit machine error of gammarel2=%27.20g fixed to be 1.0\n",gammarel2);

    gammarel2=1.0;
  }

  /*
  if(isnan(gammarel2))
    {
      printf("nan in get_m1closure_gammarel2\n");
      printf("%e %e %e\n",gamma2,gamma2a,gamma2b);
      print_4vector(Avcon);
      print_4vector(Avcov);
      //      getchar();
    }
  */

  *gammarel2return=gammarel2; 
  *deltareturn=delta=0;
  *numeratorreturn=numerator=0;
  *divisorreturn=divisor=0;
  return(0);
}


// get Erf
static int get_m1closure_Erf(void *ggg, ldouble *Avcon, ldouble gammarel2, ldouble *Erfreturn)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble alpha=geom->alpha;

  ////////////
  //
  // get initial attempt for Erf
  // If delta<0, then gammarel2=nan and Erf<RADLIMIT check below will fail as good.
  //
  ////////////
  *Erfreturn = 3.*Avcon[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

  /*
  if(isnan(*Erfreturn)) 
    {
      printf("nan in get_m1closure_Erf\n %d %d %d %e %e\n",geom->ix,geom->iy,geom->iz,*Erfreturn,gammarel2);
      print_4vector(Avcon);
      getchar();
    }
  */

  return(0);
}

// get contravariant relative 4-velocity in lab frame
static int get_m1closure_urfconrel(int verbose, 
				   void *ggg, 
				   ldouble *pp, 
				   ldouble *Avcon, 
				   ldouble *Avcov, 
				   ldouble gammarel2, 
				   ldouble delta, 
				   ldouble numerator, 
				   ldouble divisor, 
				   ldouble *Erfreturn,
				   ldouble *urfconrel, 
				   int *corflag)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;

  ldouble alpha=geom->alpha;

  ldouble Erf=*Erfreturn; // get initial Erf
  ldouble Erfslow,Erffast;
  ldouble gammamax=GAMMAMAXRAD; 
  int ii,jj,kk;

  //////////////////////
  //
  // Fix-up inversion if problem with gamma (i.e. velocity) or energy density in radiation rest-frame (i.e. Erf)
  //
  //////////////////////

  // NOTE: gammarel2 just below 1.0 already fixed to be =1.0
  int nonfailure=gammarel2>=1.0 && Erf>ERADFLOOR && gammarel2<=gammamax*gammamax/GAMMASMALLLIMIT/GAMMASMALLLIMIT;
  // falilure1 : gammarel2 normal, but already Erf<ERADFLOOR (note for M1 that gammarel2>=1/4 for any reasonable chance for correct non-zero Erf
  int failure1=Avcon[0]<0.0 || (gammarel2>0.0 && gammarel2<=0.25 && delta>=0.0 && divisor!=0.0) || numerator==0.0 || gammarel2>=1.0 && delta>=0.0 && divisor!=0.0 && Erf<ERADFLOOR;
  // gamma probably around 1
  int failure2=gammarel2<1.0 && gammarel2>0.0 && delta>=0.0;
  // i.e. all else, so not really used below.
  int failure3=gammarel2>gammamax*gammamax && Erf>=ERADFLOOR || gammarel2<0.0 || delta<0.  || divisor==0.0 && numerator==0.0 || divisor==0.0 && numerator!=0.0;

  // any failure
  int failure=!nonfailure || isinf(gammarel2) || isinf(Erf);

  if(failure && (failure1==0 && failure2==0 && failure3==0)){
    printf("Undetected failure, now considered\n");
  }


  if(nonfailure){
    // get good relative velocity
    ldouble gammarel=sqrt(gammarel2);
 
    for(ii=1;ii<4;ii++)
      {	  
	urfconrel[ii] = alpha * (Avcon[ii] + 1./3.*Erf*GG[0][ii]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
      }

    *Erfreturn=Erf; // pass back new Erf to pointer
    *corflag=0;
    return(0);
  }
  else{
    if(verbose) 
      {
	printf("failure: %d %d %d %d %d\n",nonfailure,failure1,failure2,failure3,failure);
	printf("%e %e %e\n",Erf,gammarel2,Avcov[0]);
	print_4vector(Avcon);
	print_4vector(Avcov);
	if(isnan(Erf))	getchar();
      }
    // get \gammarel=1 case
    ldouble gammarel2slow=pow(1.0+10.0*NUMEPSILON,2.0);
    ldouble Avconslow[4],Avcovslow[4],urfconrelslow[4];
    for(jj=0;jj<4;jj++)
      {
	Avconslow[jj]=Avcon[jj];
	Avcovslow[jj]=Avcov[jj];
      }
    Erfslow=Erf;
    get_m1closure_gammarel2_cold(verbose,ggg,Avconslow,Avcovslow,&gammarel2slow,&delta,&numerator,&divisor,&Erfslow,urfconrelslow);

    // get \gammarel=gammamax case
    ldouble gammarel2fast=gammamax*gammamax;
    ldouble Avconfast[4],Avcovfast[4],urfconrelfast[4];
    for(jj=0;jj<4;jj++)
      {
	Avconfast[jj]=Avcon[jj];
	Avcovfast[jj]=Avcov[jj];
      }
    Erffast=Erf;
    get_m1closure_gammarel2_cold(verbose,ggg,Avconfast,Avcovfast,&gammarel2fast,&delta,&numerator,&divisor,&Erffast,urfconrelfast);

    int usingfast=1;
    // choose by which Avcov[0] is closest to original&&
    //    if( fabs(Avcovslow[0]-Avcov[0])>fabs(Avcovfast[0]-Avcov[0])){
    if( fabs(Avcovslow[0]-Avcov[0])>fabs(Avcovfast[0]-Avcov[0])){
      usingfast=1;
      for(jj=0;jj<4;jj++)
	{
	  Avcon[jj]=Avconfast[jj];
	  Avcov[jj]=Avcovfast[jj];
	  urfconrel[jj]=urfconrelfast[jj];
	}
      gammarel2=gammarel2fast;
      Erf=Erffast;
    }
    else{
      usingfast=0;
      for(jj=0;jj<4;jj++)
	{
	  Avcon[jj]=Avconslow[jj];
	  Avcov[jj]=Avcovslow[jj];
	  urfconrel[jj]=urfconrelslow[jj];
	}
      gammarel2=gammarel2slow;
      Erf=Erfslow;
    }    

    // report
    *corflag=1;
    if(verbose) printf("CASEGEN: gammarel>gammamax (cold, usingfast=%d): gammarel2=%g Erf=%g : i=%d j=%d k=%d\n",usingfast,gammarel2,Erf,geom->ix,geom->iy,geom->iz);
  }


  *Erfreturn=Erf; // pass back new Erf to pointer

  if(verbose) 
      {
	printf("end: %e %e %e %e\n",Erf,*Erfreturn,Erfslow,Erffast);
      }

  if(isinf(Erf) || isinf(gammarel2) || isinf(urfconrel[0])|| isinf(urfconrel[1])|| isinf(urfconrel[2])|| isinf(urfconrel[3]) )
    {
      if(verbose)      printf("JONNAN: ijk=%d %d %d :  %g %g : %g %g %g : %d %d %d %d : %g %g %g %g\n",geom->ix,geom->iy,geom->iz,Erf,gammarel2,urfconrel[1],urfconrel[2],urfconrel[3],failure1,failure2,failure3,failure,Avcon[0],Avcon[1],Avcon[2],Avcon[3]);
      return -1;
  }

  return(0);
}



//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//radiative primitives: (E,\tilde u^i)
//  E - radiative energy density in the rad.rest frame
//  u^i - relative velocity of the rad.rest frame
//takes conserved R^t_mu in uu
//**********************************************************************
//**********************************************************************
int
u2p_rad_urf(ldouble *uu, ldouble *pp,void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int irf,verbose=0;

  //multi-rad-fluids
  for(irf=0;irf<NRF;irf++)
    {
      //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
      *corrected=0;

      ldouble Rij[4][4];
      ldouble urfcon[4],urfcov[4],Erf;
      //conserved - R^t_mu
      ldouble Avcov[4]={uu[EE(irf)]/gdetu,uu[FX(irf)]/gdetu,uu[FY(irf)]/gdetu,uu[FZ(irf)]/gdetu};
      ldouble Avcon[4];
      //indices up - R^tmu
      indices_12(Avcov,Avcon,GG);

      ldouble gammarel2,delta,numerator,divisor;

      // get \gamma^2 for relative 4-velocity
      if(get_m1closure_gammarel2(verbose,ggg,Avcon,Avcov,&gammarel2,&delta,&numerator,&divisor)<0)
	{
	  printf("get_m1closure_gammarel2 failed\n");
	  
	  /*
	  // compute old \gammarel using pp[]
	  urfcon[0]=0.;
	  urfcon[1]=pp[7]; //single fluid only!
	  urfcon[2]=pp[8]; 
	  urfcon[3]=pp[9];
	  conv_vels(urfcon,urfcon,VELPRIMRAD,VELR,gg,GG);
	  ldouble qsq=0.;
	  int i,j;
	  for(i=1;i<4;i++)
	    for(j=1;j<4;j++)
	      qsq+=urfcon[i]*urfcon[j]*gg[i][j];
	  ldouble gammaprev=sqrt(1.+qsq);

	  printf("correcting g2: %e -> %e\n",gammarel2,gammaprev);

	  gammarel2 = gammaprev*gammaprev;
	  */
	  return -1;
	}

      // get E in radiation frame
      get_m1closure_Erf(ggg,Avcon,gammarel2,&Erf);

      //if(verbose) printf("erf init: %e\n",Erf);

      // get relative 4-velocity
      if(get_m1closure_urfconrel(verbose,ggg,pp,Avcon,Avcov,gammarel2,delta,numerator,divisor,&Erf,urfcon,corrected)<0)
	{
	  printf("get_m1closure_urfconrel failed\n");
	  return -1;
	}

      /*
      if(Erf/pp[EE(irf)]>1.e3 && Erf>1.e-8) 
	{
	  print_4vector(Avcon);
	  print_4vector(Avcov);
	  printf("%d %d %d\n",geom->ix,geom->iy,geom->iz);
	  printf("erf: %d %e->%e %e\n",*corrected,pp[6],Erf,urfcon[0]);
	  getchar();	  
	}
      */
      
      conv_vels(urfcon,urfcon,VELR,VELPRIMRAD,gg,GG);

      //new primitives
      pp[EE(irf)]=Erf;
      pp[FX(irf)]=urfcon[1];
      pp[FY(irf)]=urfcon[2];
      pp[FZ(irf)]=urfcon[3];
    }

  return 0;
}

//**********************************************************************
//**********************************************************************
//basic conserved to primitives solver for radiation
//uses M1 closure in arbitrary frame/metric
//radiative primitives: (E,\tilde u^i)
//  E - radiative energy density in the rad.rest frame
//  u^i - relative velocity of the rad.rest frame
//takes conserved R^t_mu in uu
//**********************************************************************
//**********************************************************************
int
u2p_rad_urf_old(ldouble *uu, ldouble *pp,void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4],gdet,gdetu;
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  int irf;

  //multi-rad-fluids
  for(irf=0;irf<NRF;irf++)
    {

      //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
      *corrected=0;

      int verbose=0,debug=0;
      int i,j;
      ldouble Rij[4][4];
      ldouble urfcon[4],urfcov[4],Erf;
      ldouble alpha = sqrt(-1./GG[0][0]);
      //conserved - R^t_mu
      ldouble Av[4]={uu[EE(irf)]/gdetu,uu[FX(irf)]/gdetu,uu[FY(irf)]/gdetu,uu[FZ(irf)]/gdetu};
      //indices up - R^tmu
      indices_12(Av,Av,GG);

#if(0)
      //ortonormal classic formulation
      ldouble Avlen2=Av[1]*Av[1]+Av[2]*Av[2]+Av[3]*Av[3];
      if(Avlen2>Av[0]*Av[0])
	{
	  ldouble AE=sqrt(Avlen2/Av[0]/Av[0]);
	  Av[1]/=1.00001*AE;
	  Av[2]/=1.00001*AE;
	  Av[3]/=1.00001*AE;
	  //printf("F/E: %e %e %e f:%e\n",Av[1]/Av[0],Av[2]/Av[0],Av[3]/Av[0],sqrt(Av[1]*Av[1]+Av[2]*Av[2]+Av[3]*Av[3])/Av[0]);getchar();
	}

      //g_munu R^tmu R^tnu
      ldouble gRR=gg[0][0]*Av[0]*Av[0]+gg[0][1]*Av[0]*Av[1]+gg[0][2]*Av[0]*Av[2]+gg[0][3]*Av[0]*Av[3]+
	gg[1][0]*Av[1]*Av[0]+gg[1][1]*Av[1]*Av[1]+gg[1][2]*Av[1]*Av[2]+gg[1][3]*Av[1]*Av[3]+
	gg[2][0]*Av[2]*Av[0]+gg[2][1]*Av[2]*Av[1]+gg[2][2]*Av[2]*Av[2]+gg[2][3]*Av[2]*Av[3]+
	gg[3][0]*Av[3]*Av[0]+gg[3][1]*Av[3]*Av[1]+gg[3][2]*Av[3]*Av[2]+gg[3][3]*Av[3]*Av[3];
 
      //the quadratic equation for u^t of the radiation rest frame (urf[0])
      //supposed to provide two roots for (u^t)^2 of opposite signs
      ldouble a,b,c,delta,gamma2;
      a=16.*gRR;
      b=8.*(gRR*GG[0][0]+Av[0]*Av[0]);
      c=gRR*GG[0][0]*GG[0][0]-Av[0]*Av[0]*GG[0][0];
      delta=b*b-4.*a*c;
      gamma2=  (-b-sqrt(delta))/2./a;
      //if unphysical try the other root
      if(gamma2<0.) gamma2=  (-b+sqrt(delta))/2./a; 

      //cap on u^t
      ldouble gammamax=GAMMAMAXRAD;

      //gamma in relative velocity definition
      ldouble gammarel2=gamma2*alpha*alpha;
      if((gammarel2-1.0)>-1E-13 && gammarel2<1.0) gammarel2=1.0;
      
      debug=0;
      Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

      if(Erf<EEFLOOR)
	{
	  Erf=EEFLOOR;      printf("imposing EEFLOOR 1\n");

	  urfcon[0]=0.;
	  urfcon[1]=0.;
	  urfcon[2]=0.;
	  urfcon[3]=0.;
	  if(verbose) {printf("nocapbad: gammarel2=%g\n",gammarel2);}
	}				
      else
	{ 
	  if(1) //VELR
	    {
	      ldouble gammarel=sqrt(gammarel2);
	      for(i=1;i<4;i++)
		{	  
		  urfcon[i] = alpha * (Av[i] + 1./3.*Erf*GG[0][i]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
		}
	    }
	  else
	    {
	      //relative velocity
	      ldouble gamma=urfcon[0]*alpha;
	      for(i=1;i<4;i++)
		{	  
		  urfcon[i]=(3.*Av[i]-Erf*GG[0][i])/(3.*Av[0]-Erf*GG[0][0])/alpha-GG[0][i]/GG[0][0]/alpha;
		  urfcon[i]*=gamma;
		}
	      urfcon[0]=0.;
	    }
	}


#else


      //g_munu R^tmu R^tnu
      ldouble gRR=gg[0][0]*Av[0]*Av[0]+gg[0][1]*Av[0]*Av[1]+gg[0][2]*Av[0]*Av[2]+gg[0][3]*Av[0]*Av[3]+
	gg[1][0]*Av[1]*Av[0]+gg[1][1]*Av[1]*Av[1]+gg[1][2]*Av[1]*Av[2]+gg[1][3]*Av[1]*Av[3]+
	gg[2][0]*Av[2]*Av[0]+gg[2][1]*Av[2]*Av[1]+gg[2][2]*Av[2]*Av[2]+gg[2][3]*Av[2]*Av[3]+
	gg[3][0]*Av[3]*Av[0]+gg[3][1]*Av[3]*Av[1]+gg[3][2]*Av[3]*Av[2]+gg[3][3]*Av[3]*Av[3];
 
      //the quadratic equation for u^t of the radiation rest frame (urf[0])
      //supposed to provide two roots for (u^t)^2 of opposite signs
      ldouble a,b,c,delta,gamma2;
      a=16.*gRR;
      b=8.*(gRR*GG[0][0]+Av[0]*Av[0]);
      c=gRR*GG[0][0]*GG[0][0]-Av[0]*Av[0]*GG[0][0];
      delta=b*b-4.*a*c;
      gamma2=  (-b-sqrt(delta))/2./a;
      //if unphysical try the other root
      if(gamma2<0.) gamma2=  (-b+sqrt(delta))/2./a; 

      //cap on u^t
      ldouble gammamax=GAMMAMAXRAD;

      //gamma in relative velocity definition
      ldouble gammarel2=gamma2*alpha*alpha;

      //to avoid the machine precision error
      if((gammarel2-1.0)>-1E-13 && gammarel2<1.0) gammarel2=1.0;

      if((gammarel2>1.0*gammamax*gammamax || gammarel2<.9 || delta<0.) && verbose) 
	{
	  debug=1;

	  printf("ixy: %d %d irf: %d\n",geom->ix,geom->iy,irf);
	  print_Nvector(uu,NV);

	  print_4vector(Av);
	  printf("F/E: %e %e %e f:%e\n",Av[1]*sqrt(gg[1][1])/Av[0],Av[2]*sqrt(gg[2][2])/Av[0],sqrt(gg[3][3])*Av[3]/Av[0],sqrt(Av[1]*Av[1]+Av[2]*Av[2]+Av[3]*Av[3])/Av[0]);
	  printf("gamma2: %.20e %e\n", (-b-sqrt(delta))/2./a, (-b+sqrt(delta))/2./a); 
	  printf("delta: %e gRR: %e Erf: %e\n",delta,gRR,3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0));
      
	  //	  if(geom->ix==0&&geom->iy==0&&irf==1)getchar();//getchar();
	  
	}

      if(gammarel2>1.0*gammamax*gammamax || gammarel2<0. || delta<0.) 
	{      
      
	  //top cap
	  *corrected=1;
	  if(verbose) printf("topcap %e %e %e \n",gammarel2,gammamax*gammamax,delta);
			 
	  //urfcon[0]=gammamax;
	  ldouble gammarel;
      
	  gammarel=gammamax;

	  gammarel2=gammarel*gammarel;

	  //proper direction for the radiation rest frame, will be normalized later      
	  //Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
	  Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

	  if(Erf<EEFLOOR)
	    {
	  
	      Erf=EEFLOOR;      printf("imposing EEFLOOR 2\n");

	      urfcon[0]=0.;
	      urfcon[1]=0.;
	      urfcon[2]=0.;
	      urfcon[3]=0.;
	      if(verbose) {printf("topcapErfneg: gammarel2=%g gamma2=%g\n",gammarel2,gamma2);}
	    }					
	  else
	    {
	      if(1){ //VELR
		// lab-frame radiation relative 4-velocity
		ldouble Aradrel[4];
		for(i=1;i<4;i++)
		  Aradrel[i] = alpha * (Av[i] + 1./3.*Erf*GG[0][i]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);

		// skipping as Aradrel[] constructed assuming gammamax
		// compute \gammarel using this
		ldouble qsq=0.;
		for(i=1;i<4;i++)
		  for(j=1;j<4;j++)
		    qsq+=Aradrel[i]*Aradrel[j]*gg[i][j];
		ldouble gammatemp=sqrt(1.+qsq);

		//MYFUN(gamma_calc_fromuconrel(Aradrel,ptrgeom,&gammatemp,&qsqtemp),"ucon_calc_rel4vel_fromuconrel: gamma_calc_fromuconrel failed\n","phys.tools.rad.c",1);

		// now rescale Aradrel[i] so will give desired \gammamax
		for(i=1;i<4;i++)
		  {
		    Aradrel[i] *= (gammamax/gammatemp);
		  }
	    

		for(i=1;i<4;i++)
		  {
		    urfcon[i]=Aradrel[i];
		  }
		if(verbose) {printf("topcapgamma Erf=%g gammaorg=%.7e gammatemp=%g gammanow=%g\n",Erf,gamma2/alpha/alpha,gammatemp,gammamax);}
	      }
	      else if(0){ //going through VEL4
		ldouble Arad[4];
		for(i=1;i<4;i++)
		  {
		    Arad[i]=(Av[i]-1./3.*Erf*GG[0][i])/(4./3.*Erf*gammamax);
		  }
       
		//is normalized now
		ldouble Afac;
		a=0.; c=0.; b=0.;
		for(i=1;i<4;i++)
		  {
		    a+=Arad[i]*Arad[i]*gg[i][i];
		    b+=2.*Arad[i]*gg[0][i]*gammamax;
		  }
		c=gg[0][0]*gammamax*gammamax+1.;
		delta=b*b-4.*a*c;
		Afac= (-b+sqrt(delta))/2./a;
      
		urfcon[0]=gammamax;
		urfcon[1]=Afac*Arad[1];
		urfcon[2]=Afac*Arad[2];
		urfcon[3]=Afac*Arad[3];

		//converting to relative four velocity
		conv_vels(urfcon,urfcon,VEL4,VELPRIMRAD,gg,GG);

		if(verbose) {printf("topcapolek: Erf=%g Afac=%g Arad123=%g %g %g\n",Erf,Afac,Arad[1],Arad[2],Arad[3]);}
	      }// end Olek method
	    }// end else if Erf>0
	}
      else if(gammarel2<1.)
	{
	  //low cap
     
	  //this usually happens when gamma=0.99999999 so enforcing 1. makes sense
	  *corrected=0;
	  if(verbose) {printf("midcapalt: gamma: %.20e\n",sqrt(gammarel2));}
	  if(0)
	    {
	      //zeros for relative velocity
	      urfcon[0]=urfcon[1]=urfcon[2]=urfcon[3]=0.;
	  
	      //calculating time component of lab 4-vel
	      ldouble gammarel=1.0;
	      ldouble urflab[4];
	  
	      urflab[0]=gammarel/alpha;
	  
	      //radiative energy density in the radiation rest frame
	      Erf=3.*Av[0]/(4.*urflab[0]*urflab[0]+GG[0][0]);
	    }
	  else //in VELR
	    {
	      gammarel2=1.0;
	      Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM
	      //zeros for relative velocity
	      urfcon[0]=urfcon[1]=urfcon[2]=urfcon[3]=0.;
	    }

	  if(verbose) {printf("midcapalt: Erf=%g\n",Erf);}
	 
	  if(Erf<EEFLOOR)
	    { 
	      // Can't have Erf<0.  Like floor on internal energy density.  If leave Erf<0, then will drive code crazy with free energy.
	      if(verbose) {printf("midcapaltnegErf: Erf=%g\n",Erf);}
	      Erf=EEFLOOR;      printf("imposing EEFLOOR 3\n");

	    }	
	}
      else
	{
	  //regular calculation
	  //urfcon[0]=sqrt(gamma2);
	  //radiative energy density in the radiation rest frame
	  //Erf=3.*Av[0]/(4.*urfcon[0]*urfcon[0]+GG[0][0]);
	  Erf=3.*Av[0]*alpha*alpha/(4.*gammarel2-1.0);  // JCM

	  if(Erf<EEFLOOR)
	    {
	      Erf=EEFLOOR;      printf("imposing EEFLOOR 4\n");

	      urfcon[0]=0.;
	      urfcon[1]=0.;
	      urfcon[2]=0.;
	      urfcon[3]=0.;

	      *corrected=1;
	      if(verbose) {printf("nocapbad: gammarel2=%g\n",gammarel2);}
	    }	
	  else
	    {
	      if(1) //VELR
		{
		  ldouble gammarel=sqrt(gammarel2);
		  for(i=1;i<4;i++)
		    {	  
		      urfcon[i] = alpha * (Av[i] + 1./3.*Erf*GG[0][i]*(4.0*gammarel2-1.0) )/(4./3.*Erf*gammarel);
		    }
		}
	      else
		{
		  //relative velocity
		  ldouble gamma=urfcon[0]*alpha;
		  for(i=1;i<4;i++)
		    {	  
		      urfcon[i]=(3.*Av[i]-Erf*GG[0][i])/(3.*Av[0]-Erf*GG[0][0])/alpha-GG[0][i]/GG[0][0]/alpha;
		      urfcon[i]*=gamma;
		    }
		  urfcon[0]=0.;
		}
	    }
	}
    
#endif

      conv_vels(urfcon,urfcon,VELR,VELPRIMRAD,gg,GG);
  
      /*
   if(debug==1 && verbose)
     {
#ifndef LABRADFLUXES
       ldouble Rij[4][4];
       calc_Rij(pp,ggg,Rij);
       trans22_cc2on(Rij,Rij,tup);
       ldouble Avv[4]={Rij[0][0],Rij[0][1], Rij[0][2],Rij[0][3]};
       print_4vector(Avv);
       printf("F/E: %e %e %e f:%e\n",Avv[1]/Avv[0],Avv[2]/Avv[0],Avv[3]/Avv[0],sqrt(Avv[1]*Avv[1]+Avv[2]*Avv[2]+Avv[3]*Avv[3])/Avv[0]);
#endif

       getchar();
       }
      */
      
      //new primitives
      pp[EE(irf)]=Erf;
      pp[FX(irf)]=urfcon[1];
      pp[FY(irf)]=urfcon[2];
      pp[FZ(irf)]=urfcon[3];
   
      if(debug==1 && verbose)
	{
	  ldouble Rij[4][4];
	  calc_Rij(pp,ggg,Rij);
	  trans22_cc2on(Rij,Rij,tup);
	  ldouble Avv[4]={Rij[0][0],Rij[0][1], Rij[0][2],Rij[0][3]};
	  print_4vector(Avv);
	  printf("F/E: %e %e %e f:%e\n",Avv[1]/Avv[0],Avv[2]/Avv[0],Avv[3]/Avv[0],sqrt(Avv[1]*Avv[1]+Avv[2]*Avv[2]+Avv[3]*Avv[3])/Avv[0]);
 
	  getchar();
	}
    }

  return 0;
}

int
u2p_rad(ldouble *uu, ldouble *pp, void *ggg, int *corrected)
{
  //whether primitives corrected for caps, floors etc. - if so, conserved will be updated
  *corrected=0;

#ifdef LABRADFLUXES
  //primitives = R^t_mu
  pp[EE0]=uu[EE0];
  pp[FX0]=uu[FX0];
  pp[FY0]=uu[FY0];
  pp[FZ0]=uu[FZ0];
  return 0;
#endif

#ifdef EDDINGTON_APR
  int ii;
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],gdet,gdetu,W;
  gg=geom->gg;
  GG=geom->GG;
  gdet=geom->gdet;gdetu=gdet;
#if (GDETIN==0) //gdet out of derivatives
  gdetu=1.;
#endif

  ldouble ucov[4],ucon[4]={0,pp[2],pp[3],pp[4]};
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  indices_21(ucon,ucov,gg);
  W=gdet*ucon[0];
  ldouble Rcon[4],Rcov[4]={uu[EE0]/gdetu*gdet, uu[FX0]/gdetu*gdet, uu[FY0]/gdetu*gdet, uu[FZ0]/gdetu*gdet};
  indices_12(Rcov,Rcon,GG);
  ldouble EE, Fcon[4],Fcov[4];

  //Fragile's formulae
  EE=3.*((Rcon[0]+2.*ucon[0]*dot(ucov,Rcon))/(gdet*GG[0][0]-2.*W*ucon[0]));
  Fcon[0]=-1./gdet*(W*EE+dot(ucov,Rcon));
  for(ii=1;ii<4;ii++)
    Fcon[ii]=Rcon[ii]/W - Fcon[0]*ucon[ii]/ucon[0] - 4./3.*EE*ucon[ii] - GG[0][ii]*EE/3./ucon[0];
  indices_21(Fcon,Fcov,gg);

  int iter=0;
  while(dot(Fcon,Fcov)>=EE*EE && iter<50)
    {     
      iter++;
      Fcon[1]/=1.1;
      Fcon[2]/=1.1;
      Fcon[3]/=1.1;
      Fcon[0]=-1./ucov[0]*(Fcon[1]*ucov[1]+Fcon[2]*ucov[2]+Fcon[3]*ucov[3]); //F^0 u_0 = - F^i u_i

      indices_21(Fcon,Fcov,gg);

      *corrected=1;
    }
  
  if(iter>=50) printf("flux damping didn't work in u2p_rad for Edd\n");

  pp[EE0]=EE;
  pp[FX0]=Fcon[1];
  pp[FY0]=Fcon[2];
  pp[FZ0]=Fcon[3];
  return 0;
#endif
  
  //M1
  return u2p_rad_urf(uu,pp,ggg,corrected);
  //  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//numerical conserved to primitives solver for radiation
//works in ortonormal fluid frame
//used e.g. for not-frame-invariant  Eddington apr. 
//solves in 4dimensions using frame boosts etc.
int f_u2prad_num(ldouble *uu,ldouble *pp, void* ggg,ldouble *f)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  ldouble Rij[4][4];

  calc_Rij_ff(pp,Rij);  
  trans22_on2cc(Rij,Rij,tlo);  
  boost22_ff2lab(Rij,Rij,pp,gg,GG); 
  indices_2221(Rij,Rij,gg);  

  f[0]=-Rij[0][0]+uu[EE0];
  f[1]=-Rij[0][1]+uu[FX0];
  f[2]=-Rij[0][2]+uu[FY0];
  f[3]=-Rij[0][3]+uu[FZ0];

  return 0;
} 

int
print_state_u2prad_num (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .3e % .3e % .3e % .3e "
	  "f(x) = % .3e % .3e % .3e % .3e\n",
	  iter,
	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
  return 0;
}

int
u2p_rad_onff(ldouble *uu, ldouble *pp, void* ggg, int *corrected)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble pp0[NV],pporg[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;

  ldouble EPS = 1.e-6;
  ldouble CONV = 1.e-6;

  int verbose=0;

  for(i=EE0;i<NV;i++)
    {
      pporg[i]=pp[i];
    }

  //converting radiative primitives to fluid frame ortonormal
  //ad_lab2ff(pp,pp,geom);

  if(verbose!=0)   print_Nvector(uu,NV);
  do
    {
      iter++;
      for(i=EE0;i<NV;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_u2prad_num(uu,pp,geom,f1);
 
      //calculating approximate Jacobian
      for(j=0;j<4;j++)
	{
	  pp[j+EE0]=pp[j+EE0]+EPS*pp[EE0];
	    
	  if(verbose>0)    print_Nvector(pp,NV);
 	  f_u2prad_num(uu,pp,geom,f2);
	  if(verbose>0)    print_state_u2prad_num (iter,x,f2); 
     
	  for(i=0;i<4;i++)
	    {
	      J[i][j]=(f2[i] - f1[i])/(EPS*pp[EE0]);
	    }

	  pp[j+EE0]=pp0[j+EE0];
	}
	

      //inversion
      inverse_44matrix(J,iJ);

      //updating x
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+EE0];
	}

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}
      if(verbose>0)    print_state_u2prad_num (iter,x,f1); 

      for(i=0;i<4;i++)
	{
	  pp[i+EE0]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+EE0]-pp0[i+EE0]);
	  f3[i]=fabs(f3[i]/pp0[EE0]);
	}

      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	break;

      if(iter>20)
	{
	  printf("iter exceeded in u2prad_num() %d %d %d\n",geom->ix,geom->iy,geom->iz);getchar();

	  pp[EE0]=pporg[EE0];
	  pp[FX0]=pp[FY0]=pp[FZ0]=0.;
	  
	  *corrected=1;
	  return -1;

	  break;
	}
     
    }
  while(1);
  
  if(pp[EE0]<EEFLOOR) 
    {
      printf("enegative u2prad() %d %d %d\n",geom->ix,geom->iy,geom->iz); getchar();
      pp[EE0]=EEFLOOR;
      pp[7]=pp[8]=pp[9]=0.;
      *corrected=1;
      return -1;
    }
  
  //converting to lab primitives
  //no - in EDD_APR I use fluid frame fluxes as primitives
  //prad_ff2lab(pp,pp,geom);
  
  if(verbose!=0)   {print_Nvector(pp,NV);}
  if(verbose>0)   {printf("----\n");}//getchar();}

  *corrected=0;
  return 0;

}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//checks if rad primitives make sense
int
check_floors_rad(ldouble *pp, int whichvel,void *ggg)
{
  //skip floors for some time
  //return 0;

  int verbose=0;
  int ret=0;

  struct geometry *geom
    = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5];
  gg=geom->gg;
  GG=geom->GG;


#ifdef RADIATION  
#ifndef MULTIRADFLUID  
  ldouble pp2[NV];
  int iv;
  for(iv=0;iv<NV;iv++)
    pp2[iv]=pp[iv];
  //  ldouble Rij[4][4];
  //  calc_Rij(pp2,ggg,Rij);

  //to fluid frame to compare with gas
  //prad_lab2ff(pp, pp2, ggg); does not work

  //absolute EE:
  if(pp2[EE0]<ERADFLOOR) 
    {
      if(verbose) printf("rhd_floors CASE R0 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp2[0],pp2[EE0]);
      pp2[EE0]=ERADFLOOR;
      ret=-1;
     }

  #ifndef SKIPRADSOURCE
   //Erf/rho ratios (inconsistent?)
  if(pp2[EE0]<EERHORATIOMIN*pp2[0]) 
    {
      if(verbose) printf("hd_floors CASE R2 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp2[0],pp2[EE0]);
      pp2[EE0]=EERHORATIOMIN*pp2[0];
      ret=-1;
    }

  if(pp2[EE0]>EERHORATIOMAX*pp2[0]) 
    {
      if(verbose) printf("hd_floors CASE R3 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp2[0],pp2[EE0]);
      pp2[0]=1./EERHORATIOMAX*pp2[EE0];
      ret=-1;
    }
  
  //EE/uint ratios  
  if(pp2[EE0]<EEUURATIOMIN*pp2[1]) 
    {
      if(verbose) printf("hd_floors CASE R4 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp2[1],pp2[EE0]);
      pp2[EE0]=EEUURATIOMIN*pp2[1];
      ret=-1;
    }

  if(pp2[EE0]>EEUURATIOMAX*pp2[1]) 
    {
      if(verbose) printf("hd_floors CASE R5 at (%d,%d,%d): %e %e\n",geom->ix,geom->iy,geom->iz,pp2[1],pp2[EE0]);
      pp2[1]=1./EEUURATIOMAX*pp2[EE0];
      ret=-1;
    }
 #endif
 
  for(iv=0;iv<NV;iv++)
    pp[iv]=pp2[iv];

  //prad_ff2lab(pp2, pp, ggg);  
#endif
#endif
 

  return ret;

  //skip checking velocities - primitives should work fine

  int i,j,k,correct;
  //velocities
  ldouble a,b,c,delta;
  ldouble u1[4],u2[4];
  if(whichvel==VEL4)
    {
      //assumes u^t unknown
      u1[0]=0.;
      u1[1]=pp[2];
      u1[2]=pp[3];
      u1[3]=pp[4];
      a=gg[0][0];
      b=0.;
      c=1.;
      for(i=1;i<4;i++)
	{
	  b+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      c+=u1[i]*u1[j]*gg[i][j];
	    }
	}
      delta=b*b-4.*a*c;
      correct=0;
      if(delta<0.) 
	correct=1;
      else
	{
	  u1[0]=(-b-sqrt(delta))/2./a;
	  if(u1[0]<1.) u1[0]=(-b+sqrt(delta))/2./a;
	  if(u1[0]<1.) 
	    correct=1;
	}
    }
  else if(whichvel==VEL3)
    {
      //assumes u^t unknown
      u1[0]=0.;
      u1[1]=pp[2];
      u1[2]=pp[3];
      u1[3]=pp[4];
      a=b=0.;
      for(i=1;i<4;i++)
	{
	  a+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      b+=u1[i]*u1[j]*gg[i][j];
	    }
	}

      correct=0;
      if(-1./(gg[0][0]+a+b)<0.)
	correct=1;
    }
  else if(whichvel==VELR)
    {
      

      u1[0]=0.;
      u1[1]=pp[2];
      u1[2]=pp[3];
      u1[3]=pp[4];

      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      ldouble alpha2=-1./GG[0][0];
      
      correct=0;
      if(gamma2/alpha2<0. || gamma2>GAMMAMAXHD*GAMMAMAXHD)
	{
	  if(verbose) 
	    {
	      printf("hd_floors CASE 4 %e %e %e\n",gamma2/alpha2,gamma2,GAMMAMAXHD*GAMMAMAXHD);
	      print_4vector(u1);//getchar();
	    }
	  correct=1;
	}
    }
  else
    {
      my_err("VEL not implemented in check_floors_hd()\n");
    }

  if(correct==0) return ret;

  //converting to whichvel
  conv_vels(u1,u1,whichvel,VEL4,gg,GG);

  //correcting and imposing gammamax keeping the direction given by spatial components
  ldouble Afac;
  ldouble gammamaxhd=GAMMAMAXHD;

  if(whichvel==VELR)
    {
      //convert spatial u1 to lab-frame
      for(i=1;i<4;i++)
	u1[i]=u1[i]+gammamaxhd*GG[0][i]/GG[0][0];
    }

  //normalizing to u^t=gammamaxhd
  c=0.; b=0.;
  for(i=1;i<4;i++)
    {
      a+=u1[i]*u1[i]*gg[i][i];
      b+=2.*u1[i]*gg[0][i]*gammamaxhd;
    }
  c=gg[0][0]*gammamaxhd*gammamaxhd+1.;
  delta=b*b-4.*a*c;
      
  Afac= (-b+sqrt(delta))/2./a;

  u2[0]=gammamaxhd;
  u2[1]=Afac*u1[1];
  u2[2]=Afac*u1[2];
  u2[3]=Afac*u1[3];

   //converting to whichvel
  conv_vels(u2,u2,VEL4,whichvel,gg,GG);

  /*
  print_4vector(u1);
  print_4vector(u2);
  printf("which: %d det: %e\n",whichvel,dot(u2,u1));
  //  getchar();
  */

  //back to primitives
  pp[2]=u2[1];
  pp[3]=u2[2];
  pp[4]=u2[3];

  //to let the others know to correct conserved
  ret=-1;
  return ret;
}

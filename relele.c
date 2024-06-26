//KORAL - relele.c
//some relativistic routines

#include "ko.h" 

int //calculates only ucon, assumes ut unknown where applicable
conv_vels(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
#ifdef NONRELMHD //only four-velocity used;
  u2[1]=u1[1];u2[2]=u1[2];u2[3]=u1[3];
  fill_utinucon(u2,gg,GG); //calculates the proper u^t for four-velocity
  return 0;
#endif

  ldouble ucov[4];
  ldouble alpgam;
  if(which1==VELR)
    {
      int i,j;
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];

      ldouble gamma2=(1.+qsq);
      ldouble alpha2=(-1./GG[0][0]);
      alpgam=sqrt(alpha2*gamma2);
    }

  conv_vels_core(u1,u2,ucov,which1,which2,gg,GG,alpgam,0);
  return 0;
}

int //calculates only ucon, assumes ut known where applicable
conv_vels_ut(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
#ifdef NONRELMHD //only four-velocity used;
  u2[1]=u1[1];u2[2]=u1[2];u2[3]=u1[3];
  fill_utinucon(u2,gg,GG); //calculates the proper u^t for four-velocity
  return 0;
#endif

  ldouble ucov[4];

  ldouble alpgam;
  if(which1==VELR)
    {
      int i,j;
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];

      ldouble gamma2=(1.+qsq);
      ldouble alpha2=(-1./GG[0][0]);
      alpgam=sqrt(alpha2*gamma2);
    }

  conv_vels_core(u1,u2,ucov,which1,which2,gg,GG,alpgam,1);
  return 0;
}

int //calculates both ucon and ucov, assumes ut unknown where applicable
conv_vels_both(ldouble *u1,ldouble *u2con,ldouble *u2cov,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
 #ifdef NONRELMHD //only four-velocity used;
  u2con[1]=u1[1];u2con[2]=u1[2];u2con[3]=u1[3];
  fill_utinucon(u2con,gg,GG); //calculates the proper u^t for four-velocity
  indices_21(u2con,u2cov,gg);
  return 0;
  #endif

  ldouble ucov[4];

  ldouble alpgam;
  if(which1==VELR)
    {
      int i,j;
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];

      ldouble gamma2=(1.+qsq);
      ldouble alpha2=(-1./GG[0][0]);
      alpgam=sqrt(alpha2*gamma2);
    }

  conv_vels_core(u1,u2con,u2cov,which1,which2,gg,GG,alpgam,0);
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts contravariant velocities u1
//to contravariant u2con
//and covariant u2cov (if which2==VEL4)
int
conv_vels_core(ldouble *u1,ldouble *u2conout,ldouble *u2covout,int which1,int which2,ldouble gg[][5],ldouble GG[][5],ldouble alpgam,int utknown)
{
  int i,j;
  ldouble u2con[4],u2cov[4];
  int verbose=0;
  if(verbose)
    {
      printf("conv_vels: %d -> %d\n",which1,which2);
      print_4vector(u1);      
    }

  /*************** VEL3 -> VEL3 ***************/
  if(which1==VEL3 && which2==VEL3)
    {
      for(i=0;i<4;i++) u2con[i]=u1[i];      
    }
  /*************** VEL4 -> VEL4 ***************/
  else if(which1==VEL4 && which2==VEL4)
    {
      if(utknown==0)
	{
	  //assumes u^t unknown
	  ldouble a,b,c;
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
	  ldouble delta=b*b-4.*a*c;
	  if(delta<0.) my_err("delta.lt.0 in VEL4->VEL4\n");
	  u2con[0]=(-b-sqrt(delta))/2./a;
	  //TODO: more strict criterion
	  if(u2con[0]<0.) u2con[0]=(-b+sqrt(delta))/2./a;
	}
      else
	u2con[0]=u1[0];
     
      for(i=1;i<4;i++) u2con[i]=u1[i]; 

      indices_21(u2con,u2cov,gg);
    }
  /*************** VELR -> VELR ***************/
  else if(which1==VELR && which2==VELR)
    {
      for(i=0;i<4;i++) u2con[i]=u1[i];      
    }
  /*************** VEL4 -> VEL3 ***************/
  else if(which1==VEL4 && which2==VEL3)
    {
      if(utknown==0)
	{
	  //assumes u^t unknown
	  ldouble a,b,c;
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
	  ldouble delta=b*b-4.*a*c;
	  if(delta<0.) {printf("delta.lt.0 in VEL4->VEL3\n");return -1;}
	  u2con[0]=(-b-sqrt(delta))/2./a;
	  if(u2con[0]<1.) u2con[0]=(-b+sqrt(delta))/2./a;
	}
      else
	u2con[0]=u1[0];

      for(i=1;i<4;i++) u2con[i]=u1[i]/u2con[0];      
    }
  /*************** VEL3 -> VEL4 ***************/
  else if(which1==VEL3 && which2==VEL4)
    {
      ldouble a,b,c;
      a=b=0.;
      for(i=1;i<4;i++)
	{
	  a+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      b+=u1[i]*u1[j]*gg[i][j];
	    }
	}
						
      u2con[0]=sqrt(-1./(gg[0][0]+a+b));
      if(u2con[0]<1. || isnan(u2con[0]))
	{
	  printf("ut.nan in conv_vels(%d,%d) VEL3->VEL4 - returning error\n",which1,which2); //getchar();
	  return -1;//u2con[0]=1.0;
	}
      u2con[1]=u1[1]*u2con[0];
      u2con[2]=u1[2]*u2con[0];
      u2con[3]=u1[3]*u2con[0];  

      indices_21(u2con,u2cov,gg);
 
    }
  /*************** VEL3 -> VELR ***************/
  else if(which1==VEL3 && which2==VELR)
    {
      ldouble a,b,c;
      a=b=0.;
      for(i=1;i<4;i++)
	{
	  a+=2.*u1[i]*gg[0][i];
	  for(j=1;j<4;j++)
	    {
	      b+=u1[i]*u1[j]*gg[i][j];
	    }
	}
						
      u2con[0]=sqrt(-1./(gg[0][0]+a+b));
      if(u2con[0]<1. || isnan(u2con[0]))
	{
	  //	  printf("ut.nan in conv_vels(%d,%d) VEL3->VELR - returning error\n",which1,which2);
	  //	  print_4vector(u1);//getchar();
	  return -1;//u2con[0]=1.0;
	}
      //to 4-velocity
      u2con[1]=u1[1]*u2con[0];
      u2con[2]=u1[2]*u2con[0];
      u2con[3]=u1[3]*u2con[0];   
      //to relative velocity
      for(i=1;i<4;i++)
	u2con[i]=u2con[i]-u2con[0]*GG[0][i]/GG[0][0];
    }
  /*************** VEL4 -> VELR ***************/
  else if (which1==VEL4 && which2==VELR)
    {
      if(utknown==0)
	{
	  //assumes u^t unknown
	  ldouble a,b,c;
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
	  ldouble delta=b*b-4.*a*c;
	  if(delta<0.) {printf("delta.lt.0 in VEL4->VELR\n");return -1;}
	  ldouble ut1,ut2;
	  ut1=(-b-sqrt(delta))/2./a;
	  ut2=(-b+sqrt(delta))/2./a;

	  if(ut1>0. && ut2>0.) //ambiguous solution in the ergosphere
	    {
	      printf("ambigous ut in VEL4->VELR: %e %e\n",ut1,ut2);
	      u2con[0]=ut1;
	    }
	  else
	    u2con[0]=my_max(ut1,ut2);
	}
      else
	u2con[0]=u1[0];
       
      for(i=1;i<4;i++)
	u2con[i]=u1[i]-u2con[0]*GG[0][i]/GG[0][0];
    }

  /*************** VELR -> VEL4 ***************/
  else if (which1==VELR && which2==VEL4)
    {

      //test
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];

      ldouble gamma2=(1.+qsq);
      ldouble alpha2=(-1./GG[0][0]);
      
      u2con[0]=sqrt(gamma2/alpha2);
      u2con[1]=u1[1]-alpgam*GG[0][1];
      u2con[2]=u1[2]-alpgam*GG[0][2];
      u2con[3]=u1[3]-alpgam*GG[0][3];
      
      indices_21(u2con,u2cov,gg);

      /*
      ldouble u1cov[4];
      u1[0]=0.;
      int i,k;
      for(i=0;i<4;i++)
	{
	  u1cov[i]=0.;
	  for(k=0;k<4;k++)
	    {
	      u1cov[i]+=u1[k]*gg[i][k];
	    }	  
	}

      for(i=0;i<4;i++)
	u2cov[i]=u1cov[i]-alpgam*delta(0,i);

      //indices_12(u2cov,u2con,GG);
      for(i=0;i<4;i++)
	{
	  u2con[i]=0.;
	  for(k=0;k<4;k++)
	    {
	      u2con[i]+=u2cov[k]*GG[i][k];
	    }	  
	}
      */

    }
  /*************** VELR -> VEL3 ***************/
  else if (which1==VELR && which2==VEL3)
    {
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=1.+qsq;
      ldouble alpha2=-1./GG[0][0];
      u2con[0]=sqrt(gamma2/alpha2);
      for(i=1;i<4;i++)
	u2con[i]=u1[i]+u2con[0]*GG[0][i]/GG[0][0];
      u2con[1]/=u2con[0];
      u2con[2]/=u2con[0];
      u2con[3]/=u2con[0];
    }

  /*************** not supported  ***************/
  else
    {
      my_err("velocity conversion not supported.\n");
      return -1;
    }

  

  //write to index 0 only when needed

  u2conout[0]=u2con[0];
  u2conout[1]=u2con[1];
  u2conout[2]=u2con[2];
  u2conout[3]=u2con[3];

  u2covout[0]=u2cov[0];
  u2covout[1]=u2cov[1];
  u2covout[2]=u2cov[2];
  u2covout[3]=u2cov[3];
  
  if(verbose)
    {
      print_4vector(u2con);      
      printf("conv_vels done %d %d\n",which1,which2);
    }
  
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates u^t from spatial components of four-velocity u^mu
int
fill_utinucon(ldouble *u1,double gg[][5],ldouble GG[][5])
{
  //assumes u^t unknown
  ldouble a,b,c;
  int i,j;
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
  ldouble delta=b*b-4.*a*c;
  if(delta<0.) my_err("delta.lt.0 in VEL4->VEL4\n");
  u1[0]=(-b-sqrt(delta))/2./a;
  //TODO: more strict criterion
  if(u1[0]<0.) u1[0]=(-b+sqrt(delta))/2./a;
  
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates u_t from spatial components of four-velocity u_mu
int
fill_utinucov(ldouble *u1,double gg[][5],ldouble GG[][5])
{
  //TODO : this may be incorrect - verify!
  //assumes u_t unknown
  ldouble a,b,c;
  int i,j;
  a=GG[0][0];
  b=0.;
  c=1.;
  for(i=1;i<4;i++)
    {
      b+=2.*u1[i]*GG[0][i];
      for(j=1;j<4;j++)
	{
	  c+=u1[i]*u1[j]*GG[i][j];
	}
    }
  ldouble delta=b*b-4.*a*c;
  if(delta<0.) my_err("delta.lt.0 in VEL4->VEL4 in fill_utinucov\n");
  u1[0]=(-b-sqrt(delta))/2./a;
  //TODO: more strict criterion
  if(u1[0]>0.) u1[0]=(-b+sqrt(delta))/2./a;
  
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts contravariant u1 -> covariant u2 (takes advantage of VELR)
int
conv_velscov(ldouble *u1,ldouble *u2,int which1,int which2,ldouble gg[][5],ldouble GG[][5])
{
  #ifdef NONRELMHD
  printf("conv_velscov() not ready for NONRELMHD\n");exit(1);
  #endif

  int i,j;
  ldouble ut[4];
  int verbose=0;
  if(verbose)
    {
      printf("conv_velscov: %d -> %d\n",which1,which2);
      print_4vector(u1);      
    }

  if(which2!=VEL4)
    {
      my_err("conv_velscov() outputs only VEL4 u_mu\n");      
    }

 

  if(which1!=VELR)
    {
      conv_vels(u1,u2,which1,VEL4,gg,GG);
      indices_21(u2,ut,gg);
    }
  else
    {
      ldouble u1cov[4],utcov[4];
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=(1.+qsq);
      ldouble alpha2=(-1./GG[0][0]);

      u1[0]=0.;
      indices_21(u1,u1cov,gg); //lowering indices in utilda
      
      for(i=0;i<4;i++)
	ut[i]=u1cov[i]-sqrt(alpha2*gamma2)*delta(0,i);
    }

      /*
	if(which1!=VELR || GG[0][0]>0. || 1) //when out of domain do simple
    {
      conv_vels(u1,u2,which1,VEL4,gg,GG);
      indices_21(u2,ut,gg);
    }
  else
    {
      ldouble u1cov[4],utcon[4];
      ldouble qsq=0.;
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  qsq+=u1[i]*u1[j]*gg[i][j];
      ldouble gamma2=(1.+qsq);
      ldouble alpha2=(-1./GG[0][0]);

      indices_21(u1,u1cov,gg); //lowering indices in utilda

      utcon[0]=sqrt(gamma2/alpha2);
      for(i=1;i<4;i++)
	utcon[i]=u1[i]+utcon[0]*GG[0][i]/GG[0][0];
 
      for(i=0;i<4;i++)
	ut[i]=u1cov[i]-sqrt(alpha2*gamma2)*delta(0,i);

      //normalizing to ucon ucov = -1
      ldouble norm=dot(ut,utcon);

      if(fabs(norm+1)>0.1 || !isfinite(norm)) //large error
	{
	  conv_vels(u1,u2,which1,VEL4,gg,GG);
	  indices_21(u2,ut,gg);
	}
      else //correcting norm
	{
	  for(i=0;i<4;i++)
	    ut[i]/=norm;
	}
    }
  */
 
  u2[0]=ut[0];
  u2[1]=ut[1];
  u2[2]=ut[2];
  u2[3]=ut[3];

  if(verbose)
    {
      print_4vector(u2);      
      printf("conv_velscov done %d %d\n",which1,which2);
    }
  
  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts hydro velocities as above but takes full vector of primitives
//as an input and outputs to primitives the spatial components;
int
conv_velsinprims(ldouble *pp,int which1, int which2,ldouble gg[][5],ldouble GG[][5])
{
  int ret=0;
  ldouble v1[4],v2[4];
  v1[0]=0.; //not considered
  v1[1]=pp[2];
  v1[2]=pp[3];
  v1[3]=pp[4];
  ret=conv_vels(v1,v2,which1,which2,gg,GG);
  if(ret==0)
    {
      pp[2]=v2[1];
      pp[3]=v2[2];
      pp[4]=v2[3];
      return 0;
    }
  else
    {
      return -1;
    }
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns contravariant four-velocity of a normal observer
//n_mu=(-alp,0,0,0)
int
calc_normalobs_4vel(ldouble GG[][5], ldouble *ncon)
{
  ldouble alp=1./sqrt(-GG[0][0]);
  ldouble ncov[4]={-alp,0.,0.,0.};
  indices_12(ncov,ncon,GG);
  return 0.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns contravariant relative-velocity of a normal observer
int
calc_normalobs_relvel(ldouble GG[][5], ldouble *ncon)
{
  ldouble ucon[4];
  calc_normalobs_4vel(GG,ucon);

  int i;
  for(i=1;i<4;i++)
   ncon[i]=ucon[i]-ucon[0]*GG[0][i]/GG[0][0];

  return 0.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns hydro primitives for an atmosphere
//velocities already in VELPRIM
int
set_hdatmosphere(ldouble *pp,ldouble *xx,ldouble gg[][5],ldouble GG[][5],int atmtype)
{
  if(atmtype==0) //normal observer plur rho \propto r^-1.5
    {
      //normal observer
      ldouble ucon[4];
      ldouble xx2[4];
      calc_normalobs_relvel(GG,ucon);
      conv_vels(ucon,ucon,VELR,VELPRIM,gg,GG);
      pp[2]=ucon[1];
      pp[3]=ucon[2];
      pp[4]=ucon[3];
      // Bondi-like atmosphere
      coco_N(xx,xx2,MYCOORDS,BLCOORDS);
      ldouble r=xx2[1];

      ldouble rout=2.; //RHOATMMIN etc. given at rout=2

      pp[0] = RHOATMMIN*pow(r/rout,-1.5);
      pp[1] = UINTATMMIN*pow(r/rout,-2.5);

      #ifdef MAGNFIELD
      pp[B1]=pp[B2]=pp[B3]=0.;
      #endif

      return 0;
    }
  else if(atmtype==1) //normal observer plur rho \propto r^-2.0
    {
      //normal observer
      ldouble ucon[4],r;
      ldouble xx2[4];
      calc_normalobs_relvel(GG,ucon);
      conv_vels(ucon,ucon,VELR,VELPRIM,gg,GG);
      pp[2]=ucon[1];
      pp[3]=ucon[2];
      pp[4]=ucon[3];

      // Bondi-like atmosphere
      coco_N(xx,xx2,MYCOORDS,BLCOORDS);
      r=xx[1];

      ldouble rout=2.; //RHOATMMIN etc. given at rout=2

      pp[0] = RHOATMMIN*pow(r/rout,-2.0);
      pp[1] = UINTATMMIN*pow(r/rout,-2.5);

      #ifdef MAGNFIELD
      pp[B1]=pp[B2]=pp[B3]=0.;
      #endif
  
      return 0;
    }
  if(atmtype==2) //normal observer constant density/pressure
    {
      //normal observer
      ldouble ucon[4];
      ldouble xx2[4];
      calc_normalobs_relvel(GG,ucon);
      conv_vels(ucon,ucon,VELR,VELPRIM,gg,GG);
      pp[2]=ucon[1];
      pp[3]=ucon[2];
      pp[4]=ucon[3];
      pp[0] = RHOATMMIN;
      pp[1] = UINTATMMIN;

      #ifdef MAGNFIELD
      pp[B1]=pp[B2]=pp[B3]=0.;
      #endif

      return 0;
    }
  if(atmtype==3)
    {
#ifndef PAR_D
#define PAR_D 1.e0
#endif
#ifndef PAR_E
#define PAR_E 1.e-4
#endif
      ldouble r=xx[1];
      ldouble D=PAR_D/(r*r*sqrtl(2./r*(1.-2./r)));
      ldouble E=PAR_E/(pow(r*r*sqrt(2./r),GAMMA)*pow(1.-2./r,(GAMMA+1.)/4.));
      ldouble V=sqrtl(2./r)*(1.-2./r);
      ldouble W=1./sqrtl(1.-V*V*gg[1][1]);
      ldouble rho=D/W;
      ldouble uint=E;
      //corrected rho:
      rho=PAR_D/(r*r*sqrtl(2./r));    

      pp[0]=rho; pp[1]=uint; pp[2]=-V; pp[3]=pp[4]=0.; 
      conv_velsinprims(pp,VEL3,VELPRIM,gg,GG);
  
      #ifdef MAGNFIELD
      pp[B1]=pp[B2]=pp[B3]=0.;
      #endif

      return 0;
    }
  if(atmtype==4) //zero BL velocity
    {
      //normal observer
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

      // to VELPRIM
      conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);

      pp[2]=ucon[1];
      pp[3]=ucon[2];
      pp[4]=ucon[3];

      //density etc.
      ldouble r=xx2[1];

      ldouble rout=2.; //RHOATMMIN etc. given at rout=2

      pp[0] = RHOATMMIN*pow(r/rout,-1.5);
      pp[1] = UINTATMMIN*pow(r/rout,-2.5);
  
      #ifdef MAGNFIELD
      pp[B1]=pp[B2]=pp[B3]=0.;
      #endif

      return 0;
    }
 
    my_err("atmtype value not handled in set_atmosphere()\n");
  return 0.;
}





//**********************************************************************
//**********************************************************************
//**********************************************************************
//picks metric like tensor from cell face arr at ix,iy,iz
int
pick_Tb(ldouble* arr,int ix,int iy,int iz,int idim,ldouble T[][4])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=get_Tb(arr,i,j,ix,iy,iz,idim);
  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks metric like tensor from arr at ix,iy,iz
int
pick_T(ldouble* arr,int ix,int iy,int iz,ldouble T[][4])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      T[i][j]=get_T(arr,i,j,ix,iy,iz);
  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks a metric at ix,iy,iz
int
pick_g(int ix,int iy,int iz,ldouble gg[][5])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<5;j++)
      gg[i][j]=get_g(g,i,j,ix,iy,iz);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//picks gdet at ix,iy,iz
ldouble
pick_gdet(int ix,int iy,int iz)
{
  int i,j;
  return get_g(g,3,4,ix,iy,iz);
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks an inversed metric at ix,iy,iz
int
pick_G(int ix,int iy,int iz,ldouble GG[][5])
{
  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<5;j++)
      GG[i][j]=get_g(G,i,j,ix,iy,iz);

  return 0;
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks a metric at cell faces
int
pick_gb(int ix,int iy,int iz,int idim,ldouble gg[][5])
{
  ldouble g00,g03,g11,g22,g33,gdet,dlgdet0,dlgdet1,dlgdet2;

  //ix,iy,iz correspond to indices in cell-faces arrays
  if(idim==0)
    {
      gg[0][0]=get_gb(gbx,0,0,ix,iy,iz,0);
      gg[0][1]=get_gb(gbx,0,1,ix,iy,iz,0);
      gg[0][2]=get_gb(gbx,0,2,ix,iy,iz,0);
      gg[0][3]=get_gb(gbx,0,3,ix,iy,iz,0);
      gg[1][0]=get_gb(gbx,1,0,ix,iy,iz,0);
      gg[1][1]=get_gb(gbx,1,1,ix,iy,iz,0);
      gg[1][2]=get_gb(gbx,1,2,ix,iy,iz,0);
      gg[1][3]=get_gb(gbx,1,3,ix,iy,iz,0);
      gg[2][0]=get_gb(gbx,2,0,ix,iy,iz,0);  
      gg[2][1]=get_gb(gbx,2,1,ix,iy,iz,0);  
      gg[2][2]=get_gb(gbx,2,2,ix,iy,iz,0);  
      gg[2][3]=get_gb(gbx,2,3,ix,iy,iz,0);  
      gg[3][0]=get_gb(gbx,3,0,ix,iy,iz,0);
      gg[3][1]=get_gb(gbx,3,1,ix,iy,iz,0);
      gg[3][2]=get_gb(gbx,3,2,ix,iy,iz,0);
      gg[3][3]=get_gb(gbx,3,3,ix,iy,iz,0);
      gg[0][4]=get_gb(gbx,0,4,ix,iy,iz,0);
      gg[1][4]=get_gb(gbx,1,4,ix,iy,iz,0);
      gg[2][4]=get_gb(gbx,2,4,ix,iy,iz,0);
      gg[3][4]=get_gb(gbx,3,4,ix,iy,iz,0);      
    }
  if(idim==1)
    {
      gg[0][0]=get_gb(gby,0,0,ix,iy,iz,1);
      gg[0][1]=get_gb(gby,0,1,ix,iy,iz,1);
      gg[0][2]=get_gb(gby,0,2,ix,iy,iz,1);
      gg[0][3]=get_gb(gby,0,3,ix,iy,iz,1);
      gg[1][0]=get_gb(gby,1,0,ix,iy,iz,1);
      gg[1][1]=get_gb(gby,1,1,ix,iy,iz,1);
      gg[1][2]=get_gb(gby,1,2,ix,iy,iz,1);
      gg[1][3]=get_gb(gby,1,3,ix,iy,iz,1);
      gg[2][0]=get_gb(gby,2,0,ix,iy,iz,1);  
      gg[2][1]=get_gb(gby,2,1,ix,iy,iz,1);  
      gg[2][2]=get_gb(gby,2,2,ix,iy,iz,1);  
      gg[2][3]=get_gb(gby,2,3,ix,iy,iz,1);  
      gg[3][0]=get_gb(gby,3,0,ix,iy,iz,1);
      gg[3][1]=get_gb(gby,3,1,ix,iy,iz,1);
      gg[3][2]=get_gb(gby,3,2,ix,iy,iz,1);
      gg[3][3]=get_gb(gby,3,3,ix,iy,iz,1);
      gg[0][4]=get_gb(gby,0,4,ix,iy,iz,1);
      gg[1][4]=get_gb(gby,1,4,ix,iy,iz,1);
      gg[2][4]=get_gb(gby,2,4,ix,iy,iz,1);
      gg[3][4]=get_gb(gby,3,4,ix,iy,iz,1);      
    }
  if(idim==2)
    {
      gg[0][0]=get_gb(gbz,0,0,ix,iy,iz,2);
      gg[0][1]=get_gb(gbz,0,1,ix,iy,iz,2);
      gg[0][2]=get_gb(gbz,0,2,ix,iy,iz,2);
      gg[0][3]=get_gb(gbz,0,3,ix,iy,iz,2);
      gg[1][0]=get_gb(gbz,1,0,ix,iy,iz,2);
      gg[1][1]=get_gb(gbz,1,1,ix,iy,iz,2);
      gg[1][2]=get_gb(gbz,1,2,ix,iy,iz,2);
      gg[1][3]=get_gb(gbz,1,3,ix,iy,iz,2);
      gg[2][0]=get_gb(gbz,2,0,ix,iy,iz,2);  
      gg[2][1]=get_gb(gbz,2,1,ix,iy,iz,2);  
      gg[2][2]=get_gb(gbz,2,2,ix,iy,iz,2);  
      gg[2][3]=get_gb(gbz,2,3,ix,iy,iz,2);  
      gg[3][0]=get_gb(gbz,3,0,ix,iy,iz,2);
      gg[3][1]=get_gb(gbz,3,1,ix,iy,iz,2);
      gg[3][2]=get_gb(gbz,3,2,ix,iy,iz,2);
      gg[3][3]=get_gb(gbz,3,3,ix,iy,iz,2);
      gg[0][4]=get_gb(gbz,0,4,ix,iy,iz,2);
      gg[1][4]=get_gb(gbz,1,4,ix,iy,iz,2);
      gg[2][4]=get_gb(gbz,2,4,ix,iy,iz,2);
      gg[3][4]=get_gb(gbz,3,4,ix,iy,iz,2);      
  }  

  return 0;

}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//picks an inversed metric at cell faces
int
pick_Gb(int ix,int iy,int iz,int idim,ldouble gg[][5])
{
  ldouble g00,g03,g11,g22,g33,gdet,dlgdet0,dlgdet1,dlgdet2;

  //ix,iy,iz correspond to indices in cell-faces arrays
  if(idim==0)
    {
      gg[0][0]=get_gb(Gbx,0,0,ix,iy,iz,0);
      gg[0][1]=get_gb(Gbx,0,1,ix,iy,iz,0);
      gg[0][2]=get_gb(Gbx,0,2,ix,iy,iz,0);
      gg[0][3]=get_gb(Gbx,0,3,ix,iy,iz,0);
      gg[1][0]=get_gb(Gbx,1,0,ix,iy,iz,0);
      gg[1][1]=get_gb(Gbx,1,1,ix,iy,iz,0);
      gg[1][2]=get_gb(Gbx,1,2,ix,iy,iz,0);
      gg[1][3]=get_gb(Gbx,1,3,ix,iy,iz,0);
      gg[2][0]=get_gb(Gbx,2,0,ix,iy,iz,0);  
      gg[2][1]=get_gb(Gbx,2,1,ix,iy,iz,0);  
      gg[2][2]=get_gb(Gbx,2,2,ix,iy,iz,0);  
      gg[2][3]=get_gb(Gbx,2,3,ix,iy,iz,0);  
      gg[3][0]=get_gb(Gbx,3,0,ix,iy,iz,0);
      gg[3][1]=get_gb(Gbx,3,1,ix,iy,iz,0);
      gg[3][2]=get_gb(Gbx,3,2,ix,iy,iz,0);
      gg[3][3]=get_gb(Gbx,3,3,ix,iy,iz,0);
      gg[0][4]=get_gb(Gbx,0,4,ix,iy,iz,0);
      gg[1][4]=get_gb(Gbx,1,4,ix,iy,iz,0);
      gg[2][4]=get_gb(Gbx,2,4,ix,iy,iz,0);
      gg[3][4]=get_gb(Gbx,3,4,ix,iy,iz,0);      
    }
  if(idim==1)
    {
      gg[0][0]=get_gb(Gby,0,0,ix,iy,iz,1);
      gg[0][1]=get_gb(Gby,0,1,ix,iy,iz,1);
      gg[0][2]=get_gb(Gby,0,2,ix,iy,iz,1);
      gg[0][3]=get_gb(Gby,0,3,ix,iy,iz,1);
      gg[1][0]=get_gb(Gby,1,0,ix,iy,iz,1);
      gg[1][1]=get_gb(Gby,1,1,ix,iy,iz,1);
      gg[1][2]=get_gb(Gby,1,2,ix,iy,iz,1);
      gg[1][3]=get_gb(Gby,1,3,ix,iy,iz,1);
      gg[2][0]=get_gb(Gby,2,0,ix,iy,iz,1);  
      gg[2][1]=get_gb(Gby,2,1,ix,iy,iz,1);  
      gg[2][2]=get_gb(Gby,2,2,ix,iy,iz,1);  
      gg[2][3]=get_gb(Gby,2,3,ix,iy,iz,1);  
      gg[3][0]=get_gb(Gby,3,0,ix,iy,iz,1);
      gg[3][1]=get_gb(Gby,3,1,ix,iy,iz,1);
      gg[3][2]=get_gb(Gby,3,2,ix,iy,iz,1);
      gg[3][3]=get_gb(Gby,3,3,ix,iy,iz,1);
      gg[0][4]=get_gb(Gby,0,4,ix,iy,iz,1);
      gg[1][4]=get_gb(Gby,1,4,ix,iy,iz,1);
      gg[2][4]=get_gb(Gby,2,4,ix,iy,iz,1);
      gg[3][4]=get_gb(Gby,3,4,ix,iy,iz,1);      
    }
  if(idim==2)
    {
      gg[0][0]=get_gb(Gbz,0,0,ix,iy,iz,2);
      gg[0][1]=get_gb(Gbz,0,1,ix,iy,iz,2);
      gg[0][2]=get_gb(Gbz,0,2,ix,iy,iz,2);
      gg[0][3]=get_gb(Gbz,0,3,ix,iy,iz,2);
      gg[1][0]=get_gb(Gbz,1,0,ix,iy,iz,2);
      gg[1][1]=get_gb(Gbz,1,1,ix,iy,iz,2);
      gg[1][2]=get_gb(Gbz,1,2,ix,iy,iz,2);
      gg[1][3]=get_gb(Gbz,1,3,ix,iy,iz,2);
      gg[2][0]=get_gb(Gbz,2,0,ix,iy,iz,2);  
      gg[2][1]=get_gb(Gbz,2,1,ix,iy,iz,2);  
      gg[2][2]=get_gb(Gbz,2,2,ix,iy,iz,2);  
      gg[2][3]=get_gb(Gbz,2,3,ix,iy,iz,2);  
      gg[3][0]=get_gb(Gbz,3,0,ix,iy,iz,2);
      gg[3][1]=get_gb(Gbz,3,1,ix,iy,iz,2);
      gg[3][2]=get_gb(Gbz,3,2,ix,iy,iz,2);
      gg[3][3]=get_gb(Gbz,3,3,ix,iy,iz,2);
      gg[0][4]=get_gb(Gbz,0,4,ix,iy,iz,2);
      gg[1][4]=get_gb(Gbz,1,4,ix,iy,iz,2);
      gg[2][4]=get_gb(Gbz,2,4,ix,iy,iz,2);
      gg[3][4]=get_gb(Gbz,3,4,ix,iy,iz,2);      
  }  

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//prints primitives
int
print_p(ldouble *p)
{
  printf("rho:   %10e\nuu:    %10e\nvr:    %10e\nvth:   %10e\nvph:   %10e\nS:     %10e\n",p[0],p[1],p[2],p[3],p[4],p[5]);
#ifdef RADIATION
  printf("E:     %10e\nFx:    %10e\nFy:    %10e\nFz:    %10e\n\n",p[EE0],p[FX0],p[FY0],p[FZ0]);
#endif
  return 0;
}

//prints conserved
int
print_u(ldouble *u)
{
  printf("rhout: %10e\nTtt:   %10e\nTtr:   %10e\nTtth:  %10e\nTtph:  %10e\nSut:   %10e\n",u[0],u[1]-u[0],u[2],u[3],u[4],u[5]);
#ifdef RADIATION
  printf("Rtt:   %10e\nRt1:   %10e\nRt2:   %10e\nRt3:   %10e\n\n",u[EE0],u[FX0],u[FY0],u[FZ0]);
#endif
  return 0;
}



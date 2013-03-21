//KORAL - frames.c
//boosting, moving indices etc.

#include "ko.h"

/*****************************************************************/
/********** all primitives between coordinates  ***************/
/*****************************************************************/
int 
trans_pall_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, ldouble gg1[][5], ldouble GG1[][5], ldouble gg2[][5], ldouble GG2[][5])
{
  trans_phd_coco(pp1, pp2, CO1,CO2,xxvec,gg1,GG1,gg2,GG2);
#ifdef RADIATION
  trans_prad_coco(pp1, pp2, CO1,CO2,xxvec,gg1,GG1,gg2,GG2);
#endif
  return 0;
}

/*****************************************************************/
/********** hydro primitives (E,F^i) between coordinates  *******/
/********** does not touch radiative primitives ***********************/
/*****************************************************************/
int 
trans_phd_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, ldouble gg1[][5], ldouble GG1[][5], ldouble gg2[][5], ldouble GG2[][5])
{
  int i;
  for(i=NVHD;i<NV;i++)
    pp2[i]=pp1[i];
     
  if(CO1==CO2)
    {
      pp2[0]=pp1[0];
      pp2[1]=pp1[1];
      pp2[2]=pp1[2];
      pp2[3]=pp1[3];
      pp2[4]=pp1[4];
    }
  else
    {
      pp2[0]=pp1[0];
      pp2[1]=pp1[1];
      //velocity in CO1
      ldouble ucon[4];
      ucon[0]=0;
      ucon[1]=pp1[2];
      ucon[2]=pp1[3];
      ucon[3]=pp1[4];

      conv_vels(ucon,ucon,VELPRIM,VEL4,gg1,GG1);
      //converting to CO2
      trans2_coco(xxvec,ucon,ucon,CO1,CO2);
      //to VELPRIM
      conv_vels(ucon,ucon,VEL4,VELPRIM,gg2,GG2);

      pp2[2]=ucon[1]; 
      pp2[3]=ucon[2];
      pp2[4]=ucon[3];
    }
  
  return 0;
}
 
/*****************************************************************/
/********** radiative primitives (E,F^i) between coordinates  *******/
/********** does not touch hydro primitives ***********************/
/*****************************************************************/
int 
trans_prad_coco(ldouble *pp1, ldouble *pp2, int CO1,int CO2, ldouble *xxvec, ldouble gg1[][5], ldouble GG1[][5], ldouble gg2[][5], ldouble GG2[][5])
{
  int i;
  for(i=0;i<NVHD;i++)
    pp2[i]=pp1[i];

  if(CO1==CO2)
    {
      pp2[6]=pp1[6];
      pp2[7]=pp1[7];
      pp2[8]=pp1[8];
      pp2[9]=pp1[9];
    }
  else
    {
      //Erf unchanged
      pp2[6]=pp1[6];

       //velocity in CO1
      ldouble ucon[4];
      ucon[0]=0;
      ucon[1]=pp1[7];
      ucon[2]=pp1[8];
      ucon[3]=pp1[9];

      conv_vels(ucon,ucon,VELPRIMRAD,VEL4,gg1,GG1);
      //converting to CO2
      trans2_coco(xxvec,ucon,ucon,CO1,CO2);
      //to VELPRIM
      conv_vels(ucon,ucon,VEL4,VELPRIMRAD,gg2,GG2);

      pp2[7]=ucon[1]; 
      pp2[8]=ucon[2];
      pp2[9]=ucon[3];

   }
  
  return 0;
}


/*****************************************************************/
/****** radiative ff primitives (E,F^i) -> primitives in lab frame  *******/
/*****************************************************************/
int prad_ff2lab(ldouble *pp1, ldouble *pp2, void* ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  /*
print_metric(gg);print_metric(GG);
print_tensor(tlo);print_tensor(tup);getchar();
  */

  ldouble Rij[4][4];
  int i,j;

  int verbose=0;
 
  calc_Rij_ff(pp1,Rij);  
  trans22_on2cc(Rij,Rij,tlo);  
  boost22_ff2lab(Rij,Rij,pp1,gg,GG); 
  indices_2221(Rij,Rij,gg);  

  for(i=0;i<NVHD;i++)
    pp2[i]=pp1[i];

  //temporarily store conserved in pp2[]
  pp2[6]=Rij[0][0];
  pp2[7]=Rij[0][1];
  pp2[8]=Rij[0][2];
  pp2[9]=Rij[0][3];

  //  print_tensor(Rij);
      
  //convert to real primitives
  int corrected;

  u2p_rad(pp2,pp2,geom,&corrected);

  return 0;
} 

//*****************************************************************/
/********** radiative primitives lab -> (E,F^i) in fluid frame *********/
/*****************************************************************/
int prad_lab2ff(ldouble *pp1, ldouble *pp2, void *ggg)
{
  struct geometry *geom
   = (struct geometry *) ggg;

  ldouble (*gg)[5],(*GG)[5],(*tlo)[4],(*tup)[4];
  gg=geom->gg;
  GG=geom->GG;
  tlo=geom->tlo;
  tup=geom->tup;

  ldouble Rij[4][4];
  int i,j;  

  calc_Rij(pp1,ggg,Rij);
  boost22_lab2ff(Rij,Rij,pp1,gg,GG);
  trans22_cc2on(Rij,Rij,tup);

  for(i=0;i<NVHD;i++)
    pp2[i]=pp1[i];

  //E,F^i
  pp2[6]=Rij[0][0];
  pp2[7]=Rij[0][1];
  pp2[8]=Rij[0][2];
  pp2[9]=Rij[0][3];

  return 0;
} 

/***********************************************************************/
/******radiative primitives (E,F^i) fluid frame -> (E,F^i) in ZAMO *********/
/***********************************************************************/
int prad_ff2zamo(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4])
{
  ldouble Rij[4][4];
  int i,j;

  calc_Rij_ff(pp1,Rij);
  boost22_ff2zamo(Rij,Rij,pp1,gg,GG,eup);

  for(i=0;i<NVHD;i++)
    pp2[i]=pp1[i];

  //(E,F^i)_ZAMO
  pp2[6]=Rij[0][0];
  pp2[7]=Rij[0][1];
  pp2[8]=Rij[0][2];
  pp2[9]=Rij[0][3];

  return 0;
} 

/***********************************************************************/
/******radiative primitives (E,F^i) fluid frame -> (E,F^i) in ZAMO *********/
/***********************************************************************/
int prad_zamo2ff(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble eup[][4])
{
  ldouble Rij[4][4];
  int i,j;

  //infact, closure in ZAMO flat space
  calc_Rij_ff(pp1,Rij);
  boost22_zamo2ff(Rij,Rij,pp1,gg,GG,eup);

  for(i=0;i<NVHD;i++)
    pp2[i]=pp1[i];

  //(E,F^i)_ff
  pp2[6]=Rij[0][0];
  pp2[7]=Rij[0][1];
  pp2[8]=Rij[0][2];
  pp2[9]=Rij[0][3];

  return 0;
} 

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//calculates Lorenz matrix for lab -> ff
int
calc_Lorentz_lab2ff(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble L[][4])
{
  int i,j,k;
  int verbose=0;

  //calculating the four-velocity of fluid in lab frame
  ldouble ucon[4],ucov[4],vpr[3];
  ucon[1]=pp[2];
  ucon[2]=pp[3];
  ucon[3]=pp[4];
  conv_vels(ucon,ucon,VELPRIM,VEL4,gg,GG);
  //covariant four-velocity
  indices_21(ucon,ucov,gg);  

  if(verbose>0) print_4vector(ucon);

  //four velocity of the lab frame
  ldouble alpha=sqrt(-1./GG[0][0]);
  ldouble wcon[4],wcov[4]={-alpha,0.,0.,0.};
  indices_12(wcov,wcon,GG);

  //wcon[0]=1./sqrt(-gg[0][0]);
  //wcon[1]=wcon[2]=wcon[3]=0.;
  //indices_21(wcon,wcov,gg);

  if(verbose>0) print_4vector(wcon);

  //temporary Om matrix
  ldouble Om[4][4];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Om[i][j]=ucon[i]*wcov[j]-wcon[i]*ucov[j];
  
  //Lorentz factor = -w^mu u_mu
  ldouble gam=-dot(wcon,ucov);

  ldouble Omsum;
  //Lorentz matrix components
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Omsum=0.;
	for(k=0;k<4;k++)
	  Omsum+=Om[i][k]*Om[k][j];
	
	L[i][j]=kron(i,j)+1./(1.+gam)*Omsum+Om[i][j];
      }
  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//calculates Lorenz matrix for ff -> lab
int
calc_Lorentz_ff2lab(ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble L[][4])
{
  int i,j,k;
  int verbose=0;

  //calculating the four-velocity of fluid in lab frame
  ldouble wcon[4],wcov[4];
  wcon[1]=pp[2];
  wcon[2]=pp[3];
  wcon[3]=pp[4];
  conv_vels(wcon,wcon,VELPRIM,VEL4,gg,GG);
  //covariant four-velocity
  indices_21(wcon,wcov,gg);  
 
  if(verbose>0) print_4vector(wcon);

  //four velocity of the lab frame
  ldouble alpha=sqrt(-1./GG[0][0]);
  ldouble ucon[4],ucov[4]={-alpha,0.,0.,0.};
  indices_12(ucov,ucon,GG);

  //ucon[0]=1./sqrt(-gg[0][0]);
  //ucon[1]=ucon[2]=ucon[3]=0.;
  //indices_21(ucon,ucov,gg);

  if(verbose>0) print_4vector(ucon);

  //temporary Om matrix
  ldouble Om[4][4];

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      Om[i][j]=ucon[i]*wcov[j]-wcon[i]*ucov[j];
  
  //Lorentz factor = -w^mu u_mu
  ldouble gam=-dot(wcon,ucov);

  ldouble Omsum;
  //Lorentz matrix components
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	Omsum=0.;
	for(k=0;k<4;k++)
	  Omsum+=Om[i][k]*Om[k][j];
	
	L[i][j]=kron(i,j)+1./(1.+gam)*Omsum+Om[i][j];
      }
  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost from lab to fluid frame
int
boost22_lab2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;

  if(verbose>0) print_tensor(T1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff(pp,gg,GG,L);

  //copying and multiplying by lapse to express T1 in ZAMO
  ldouble alpha=sqrt(-1./GG[0][0]);  
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j]*alpha;
	}
    }
  
  //multiplying by lapse to express T1 in ZAMO
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost from lab to fluid frame
int
boost22_ff2lab(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;

  if(verbose>0) print_tensor(T1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_ff2lab(pp,gg,GG,L);

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
 
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    } 

  //dividing by lapse to express T2 in no-frame
  ldouble alpha=sqrt(-1./GG[0][0]);  
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=T2[i][j]/alpha;
	}
    }

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost from lab to fluid frame
int
boost2_lab2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_4vector(A1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_lab2ff(pp,gg,GG,L);

  //copying and multiplying by lapse to express A1 in ZAMO
  ldouble alpha=sqrt(-1./GG[0][0]);  
  for(i=0;i<4;i++)
    {
      At[i]=A1[i]*alpha;
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }

  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0; 
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost from fluid to lab frame
int
boost2_ff2lab(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_4vector(A1);

  //Lorentz transformation matrix
  ldouble L[4][4];
  calc_Lorentz_ff2lab(pp,gg,GG,L);

  //copying
  for(i=0;i<4;i++)
    {
      At[i]=A1[i];
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }

  //dividing by lapse to express A2 in no-frame
  ldouble alpha=sqrt(-1./GG[0][0]);  
  for(i=0;i<4;i++)
    {      
      A2[i]=A2[i]/alpha;	
    }


  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0; 
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost ZAMO -> ortonormal fluid frame
//eup currently defined only for Kerr-like metric
int
boost2_zamo2ff(ldouble* A1,ldouble* A2,ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_tensor(eup);
  if(verbose>0) print_4vector(A1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ulab[0]=0.;
  ulab[1]=pp[2];
  ulab[2]=pp[3];
  ulab[3]=pp[4];

  conv_vels(ulab,ulab,VELPRIM,VEL4,gg,GG);
 
  if(verbose>0) print_4vector(ulab);
  
  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  if(verbose>0) print_4vector(uzamo);

  //proper velocity for ZAMO
  vpr[0]=-uzamo[1]/uzamo[0];
  vpr[1]=-uzamo[2]/uzamo[0];
  vpr[2]=-uzamo[3]/uzamo[0];

  //Lorentz transformation matrix
  ldouble L[4][4];
  
  //Lorentz factor
  ldouble vpr2=dot3(vpr,vpr); 
  ldouble gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2;
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }

  //copying
  for(i=0;i<4;i++)
    {
      At[i]=A1[i];
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }

  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost fluid frame -> ZAMO
int
boost2_ff2zamo(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_tensor(eup);
  if(verbose>0) print_4vector(A1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ulab[0]=0.;
  ulab[1]=pp[2];
  ulab[2]=pp[3];
  ulab[3]=pp[4];

  conv_vels(ulab,ulab,VELPRIM,VEL4,gg,GG);
 
  if(verbose>0) print_4vector(ulab);

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  if(verbose>0) print_4vector(uzamo);

  //proper velocity for ZAMO
  vpr[0]=uzamo[1]/uzamo[0];
  vpr[1]=uzamo[2]/uzamo[0];
  vpr[2]=uzamo[3]/uzamo[0];

  //Lorentz transformation matrix
  ldouble L[4][4];
  
  //Lorentz factor
  ldouble vpr2=dot3(vpr,vpr); 
  ldouble gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2;
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }

  //copying
  for(i=0;i<4;i++)
    {
      At[i]=A1[i];
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      A2[i]=0.;
      for(k=0;k<4;k++)
	{
	  A2[i]+=L[i][k]*At[k];
	}
    }

  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost ZAMO -> fluid frame
int
boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;
  //debug
  //  if(pp[2]<0.) verbose=1;

  if(verbose>0) print_tensor(eup);

  if(verbose) printf("sqrtg: %e %e %e %e\n",sqrt(-gg[0][0]),sqrt(gg[1][1]),sqrt(gg[2][2]),sqrt(gg[3][3]));

  if(verbose>0) print_tensor(T1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ulab[0]=0.;
  ulab[1]=pp[2];
  ulab[2]=pp[3];
  ulab[3]=pp[4];

  conv_vels(ulab,ulab,VELPRIM,VEL4,gg,GG);
 
  if(verbose>0) print_4vector(ulab);

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  if(verbose>0) print_4vector(uzamo);

  //proper velocity for ZAMO
  vpr[0]=-uzamo[1]/uzamo[0];
  vpr[1]=-uzamo[2]/uzamo[0];
  vpr[2]=-uzamo[3]/uzamo[0];

  if(verbose) print_4vector(vpr);

  //Lorentz transformation matrix
  ldouble L[4][4];
  
  //Lorentz factor
  ldouble vpr2=dot3(vpr,vpr); 
  ldouble gam=uzamo[0];


  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2;
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij Lorentz boost fluid frame -> ZAMO
int
boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble GG[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;
  
  if(verbose>0) print_tensor(eup);
  if(verbose>0) print_tensor(T1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ulab[0]=0.;
  ulab[1]=pp[2];
  ulab[2]=pp[3];
  ulab[3]=pp[4];

  conv_vels(ulab,ulab,VELPRIM,VEL4,gg,GG);
   
  if(verbose>0) print_4vector(ulab);

  //transforming 4-vector lab->zamo
  trans2_lab2zamo(ulab,uzamo,eup);

  if(verbose>0) print_4vector(uzamo);

  //proper velocity for ZAMO
  vpr[0]=uzamo[1]/uzamo[0];
  vpr[1]=uzamo[2]/uzamo[0];
  vpr[2]=uzamo[3]/uzamo[0];

  if(verbose>0) print_Nvector(vpr,3);

  //Lorentz transformation matrix
  ldouble L[4][4];
  
  //Lorentz factor
  ldouble vpr2=dot3(vpr,vpr); 
  ldouble gam=uzamo[0];

  //unchanged sign of vpr  
  L[0][0]=gam;
  L[1][0]=L[0][1]=gam*vpr[0];
  L[2][0]=L[0][2]=gam*vpr[1];
  L[3][0]=L[0][3]=gam*vpr[2];

  //Lorentz matrix components
  if(vpr2>SMALL)
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j)+vpr[i-1]*vpr[j-1]*(gam-1.)/vpr2;
    }
  else
    {
      for(i=1;i<4;i++)
	for(j=1;j<4;j++)
	  L[i][j]=kron(i,j);
    }

  //copying
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }
  
  if(verbose>0) print_tensor(L);

  //boosting
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=L[i][k]*L[j][l]*Tt[k][l];
		}
	    }
	}
    }

  if(verbose>0) print_tensor(T2);

  if(verbose>0) getchar();

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//multiplies 22 tensor T1 by 21 tensor A
//T2^ij = A^i_k A^j_l T1^kl
int
multiply22(ldouble T1[][4],ldouble T2[][4],ldouble A[][4])
{
  int i,j,k,l;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=T1[i][j];
	}
    }

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      for(l=0;l<4;l++)
		{
		  T2[i][j]+=A[i][k]*A[j][l]*Tt[k][l];
		}
	    }
	}
    }
  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//multiplies 2 vector u1 by 21 tensor A
//u2^i = A^i_j u1^j
int
multiply2(ldouble *u1,ldouble *u2,ldouble A[][4])
{
  int i,j;
  ldouble ut[4];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=A[i][j]*ut[j];
	}
    }

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation ZAMO -> lab
int
trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble elo[][4])
{
  multiply22(T1,T2,elo);
 
  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation lab -> ZAMO
int
trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble eup[][4])
{
  multiply22(T1,T2,eup);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation lab -> ZAMO
int
trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble eup[][4])
{
  multiply2(u1,u2,eup);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation ZAMO -> lab
int
trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble elo[][4])
{
  multiply2(u1,u2,elo);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation between coordinates
int
trans2_coco(ldouble *xx,ldouble *u1,ldouble *u2,int CO1, int CO2)
{
  ldouble dxdx[4][4];
  ldouble xx2[4];
  if(CO1==CO2)
    {
      u2[0]=u1[0];
      u2[1]=u1[1];
      u2[2]=u1[2];
      u2[3]=u1[3];
    }
  else if(((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==SPHCOORDS) ||
	  ((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==SPHCOORDS))
    {
      u2[0]=u1[0];
      u2[1]=u1[1];
      u2[2]=u1[2];
      u2[3]=u1[3];
    }
   else if(CO1==KSCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_KS2BL(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==KSCOORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MKS1COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS1COORDS)
    {
      dxdx_KS2MKS1(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==MCYL1COORDS && CO2==CYLCOORDS)
    {
      dxdx_MCYL12CYL(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if(CO1==CYLCOORDS && CO2==MCYL1COORDS)
    {
      dxdx_CYL2MCYL1(xx,dxdx);
      multiply2(u1,u2,dxdx);
    }
  else if (CO1==MKS1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply2(u2,u2,dxdx);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MKS1COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply2(u1,u2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS1(xx2,dxdx);
      multiply2(u2,u2,dxdx);  
    }
  else
    my_err("transformation not implemented in trans2_coco()\n");

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation between coordinates
int
trans22_coco(ldouble *xx,ldouble T1[][4],ldouble T2[][4],int CO1, int CO2)
{
  ldouble dxdx[4][4];
  ldouble xx2[4];
  if(CO1==CO2)
    {
      int i,j;
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  T2[i][j]=T1[i][j];
    }
  else if(((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==SPHCOORDS) ||
	  ((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==SPHCOORDS))
    {
      int i,j;
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  T2[i][j]=T1[i][j];
    }
  else if(CO1==KSCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_KS2BL(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==KSCOORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MKS1COORDS && CO2==KSCOORDS)
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==KSCOORDS && CO2==MKS1COORDS)
    {
      dxdx_KS2MKS1(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==MCYL1COORDS && CO2==CYLCOORDS)
    {
      dxdx_MCYL12CYL(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if(CO1==CYLCOORDS && CO2==MCYL1COORDS)
    {
      dxdx_CYL2MCYL1(xx,dxdx);
      multiply22(T1,T2,dxdx);
    }
  else if (CO1==MKS1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    {
      dxdx_MKS12KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2BL(xx2,dxdx);
      multiply22(T2,T2,dxdx);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MKS1COORDS)
    {
      dxdx_BL2KS(xx,dxdx);
      multiply22(T1,T2,dxdx);
      coco_N(xx,xx2,CO1,KSCOORDS);
      dxdx_KS2MKS1(xx2,dxdx);
      multiply22(T2,T2,dxdx);  
    }
  else
    my_err("transformation not implemented in trans22_coco()\n");

  return 0;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation ortonormal to code coordinates
int
trans22_on2cc(ldouble T1[][4],ldouble T2[][4],ldouble tlo[][4])
{
  multiply22(T1,T2,tlo);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation code coords -> ortonormal
int
trans22_cc2on(ldouble T1[][4],ldouble T2[][4],ldouble tup[][4])
{
  multiply22(T1,T2,tup);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation cc -> on
int
trans2_cc2on(ldouble *u1,ldouble *u2,ldouble tup[][4])
{
  multiply2(u1,u1,tup);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation on -> cc
int
trans2_on2cc(ldouble *u1,ldouble *u2,ldouble tlo[][4])
{
  multiply2(u1,u2,tlo);

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// T^ij -> T^i_j
int
indices_2221(ldouble T1[][4],ldouble T2[][4],ldouble gg[][5])
{
  int i,j,k;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      Tt[i][j]+=T1[i][k]*gg[k][j];
	    }	  
	}
    }

   for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}

//*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// T^i_j -> T^ij
int
indices_2122(ldouble T1[][4],ldouble T2[][4],ldouble GG[][5])
{
  int i,j,k;
  ldouble Tt[4][4];

  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  Tt[i][j]=0.;
	  for(k=0;k<4;k++)
	    {
	      Tt[i][j]+=T1[i][k]*GG[k][j];
	    }	  
	}
    }

   for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  T2[i][j]=Tt[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A^i -> A^_j
int
indices_21(ldouble A1[4],ldouble A2[4],ldouble gg[][5])
{
  int i,j,k;
  ldouble At[4];

  for(i=0;i<4;i++)
    {
      At[i]=0.;
      for(k=0;k<4;k++)
	{
	  At[i]+=A1[k]*gg[i][k];
	}	  
    }

   for(i=0;i<4;i++)
    {
	  A2[i]=At[i];
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// A_i -> A^_j
int
indices_12(ldouble A1[4],ldouble A2[4],ldouble GG[][5])
{
  int i,j,k;
  ldouble At[4];

  for(i=0;i<4;i++)
    {
      At[i]=0.;
      for(k=0;k<4;k++)
	{
	  At[i]+=A1[k]*GG[i][k];
	}	  
    }

   for(i=0;i<4;i++)
    {
	  A2[i]=At[i];
    }

  return 0;
}


/*****************************************************************/
/* prints tensor to screen *****************************************/
/*****************************************************************/
int
print_tensor(ldouble T[][4])
{
  int i;
  printf("============\n");
  for(i=0;i<4;i++)
    printf("%10.3e %10.3e %10.3e %10.3e\n",T[0][i],T[1][i],T[2][i],T[3][i]);
  printf("============\n");
  return 0;  
}

/*****************************************************************/
/* prints metric to screen *****************************************/
/*****************************************************************/
int
print_metric(ldouble T[][5])
{
  int i;
  printf("============\n");
  for(i=0;i<4;i++)
    printf("%10.3e %10.3e %10.3e %10.3e\n",T[i][0],T[i][1],T[i][2],T[i][3]);
  printf("============\n");
  return 0;  
}

/*****************************************************************/
/* prints 4vector to screen ****************************************/
/*****************************************************************/
int
print_4vector(ldouble v[4])
{
  int i;
  printf("============\n");
  printf("%10.3e %10.3e %10.3e %10.3e\n",v[0],v[1],v[2],v[3]);
  printf("============\n");
  return 0;  
}

/*****************************************************************/
/* prints Nvector to screen ****************************************/
/*****************************************************************/
int
print_Nvector(ldouble v[4],int N)
{
  int i;
  printf("============\n");
  for(i=0;i<N;i++)
  printf("%10.7e ",v[i]);
  printf("\n============\n");
  return 0;  
}

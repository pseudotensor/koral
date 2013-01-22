//KORAL - frames.c
//boosting, moving indices etc.

#include "ko.h"

/*****************************************************************/
/****** radiative ff primitives (E,F^i) -> primitives in lab frame  *******/
/*****************************************************************/
int prad_ff2lab(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble tlo[][4])
{
  ldouble Rij[4][4];
  int i,j;

  int verbose=0;
  if(fabsl(pp1[2])>0.) verbose=1;
 if(verbose) printf("poczatek prad_ff2lab\n");
 
 
  if(verbose) print_Nvector(pp1,NV);
  calc_Rij_ff(pp1,Rij);
  if(verbose) print_tensor(Rij);
  trans22_on2cc(Rij,Rij,tlo);
  if(verbose) print_tensor(tlo);
  if(verbose) print_tensor(Rij);
  boost22_ff2lab(Rij,Rij,pp1,gg,GG);

  if(verbose) print_tensor(Rij);
  indices_2221(Rij,Rij,gg);
  if(verbose) print_tensor(Rij);
  if(verbose) printf("koniec prad_ff2lab\n");
  if(verbose) getchar();

  for(i=0;i<NVHD;i++)
    pp2[i]=pp1[i];

  //temporarily store conserved in pp2[]
  pp2[6]=Rij[0][0];
  pp2[7]=Rij[0][1];
  pp2[8]=Rij[0][2];
  pp2[9]=Rij[0][3];

  //convert to real primitives
  int corrected;
  //u2p_rad(pp2,pp2,gg,GG,&corrected);

  return 0;
} 

//*****************************************************************/
/********** radiative primitives lab -> (E,F^i) in fluid frame *********/
/*****************************************************************/
int prad_lab2ff(ldouble *pp1, ldouble *pp2, ldouble gg[][5], ldouble GG[][5], ldouble tup[][4])
{
  ldouble Rij[4][4];
  int i,j;  

  calc_Rij(pp1,gg,GG,Rij);
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
  boost22_ff2zamo(Rij,Rij,pp1,gg,eup);

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
  boost22_zamo2ff(Rij,Rij,pp1,gg,eup);

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
/********** (E,F^i) ZAMO -> (E,F^i) fluid frame ********************/
/*****************************************************************/
int f_prad_zamo2ff(ldouble *ppff, ldouble *ppzamo, ldouble gg[][5], ldouble eup[][4],ldouble *f)
{
  ldouble Rij[4][4];
  calc_Rij_ff(ppff,Rij);
  boost22_ff2zamo(Rij,Rij,ppff,gg,eup);

  f[0]=-Rij[0][0]+ppzamo[6];
  f[1]=-Rij[0][1]+ppzamo[7];
  f[2]=-Rij[0][2]+ppzamo[8];
  f[3]=-Rij[0][3]+ppzamo[9];

  return 0;
} 

int
print_state_prad_zamo2ff (int iter, ldouble *x, ldouble *f)
{
  printf ("iter = %3d x = % .3Le % .3Le % .3Le % .3Le "
	  "f(x) = % .3Le % .3Le % .3Le % .3Le\n",
	  iter,
	  x[0],x[1]/x[0],x[2]/x[0],x[3]/x[0],f[0],f[1],f[2],f[3]);
  return 0;
}

int prad_zamo2ff_num(ldouble *ppzamo, ldouble *ppff, ldouble gg[][5],ldouble GG[][5], ldouble eup[][4])
{

  ldouble pp0[NV],pp[NV];
  ldouble J[4][4],iJ[4][4];
  ldouble x[4],f1[4],f2[4],f3[4];
  int i,j,k,iter=0;

  ldouble EPS = 1.e-6;
  ldouble CONV = 1.e-8;

  int verbose=0;
  
  //initial guess
  for(i=0;i<NV;i++)
    {
      pp[i]=ppzamo[i];
    }

  if(verbose!=0)   print_Nvector(ppff,NV);

  //debug
  f_prad_zamo2ff(ppzamo,pp,gg,eup,f1);
  for(i=0;i<4;i++)
    {
      x[i]=pp[i+6];
    }  
  if(verbose>0) print_state_prad_zamo2ff (iter,x,f1); 
 
  do
    {
      iter++;
      for(i=6;i<NV;i++)
	{
	  pp0[i]=pp[i];
	}

      //valueas at zero state
      f_prad_zamo2ff(pp,ppzamo,gg,eup,f1);
 
      //calculating approximate Jacobian
      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      pp[j+6]=pp[j+6]+EPS*pp[6];
	    
	      f_prad_zamo2ff(pp,ppzamo,gg,eup,f2);
     
	      J[i][j]=(f2[i] - f1[i])/(EPS*pp[6]);

	      pp[j+6]=pp0[j+6];
	    }
	}
  
      //inversion
      inverse_44matrix(J,iJ);

      //updating unknowns
      for(i=0;i<4;i++)
	{
	  x[i]=pp0[i+6];
	}      

      for(i=0;i<4;i++)
	{
	  for(j=0;j<4;j++)
	    {
	      x[i]-=iJ[i][j]*f1[j];
	    }
	}
      if(verbose>0)    print_state_prad_zamo2ff (iter,x,f1); 

      for(i=0;i<4;i++)
	{
	  pp[i+6]=x[i];
	}
  
      //test convergence
      for(i=0;i<4;i++)
	{
	  f3[i]=(pp[i+6]-pp0[i+6]);
	  f3[i]=fabs(f3[i]/pp0[6]);
	}

      if(f3[0]<CONV && f3[1]<CONV && f3[2]<CONV && f3[3]<CONV)
	break;

      if(iter>50)
	{
	  printf("iter exceeded in prad_zamo2ff()\n");
	  break;
	}
    }
  while(1);

  if(verbose!=0)   {print_Nvector(pp,NV);getchar(); }

  if(verbose>0)   printf("----\n");

  //returning prad
  for(i=0;i<NV;i++)
    {
      ppzamo[i]=pp[i];
    }

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
  ldouble wcon[4],wcov[4];
  wcon[0]=1./sqrtl(-gg[0][0]);
  wcon[1]=wcon[2]=wcon[3]=0.;
  indices_21(wcon,wcov,gg);

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
  ldouble ucon[4],ucov[4];
  ucon[0]=1./sqrtl(-gg[0][0]);
  ucon[1]=ucon[2]=ucon[3]=0.;
  indices_21(ucon,ucov,gg);

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

  if(verbose>0) print_4vector(A2);

  if(verbose>0) getchar();

  return 0; 
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//A^i Lorentz boost ZAMO -> fluid frame
int
boost2_zamo2ff(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_tensor(eup);
  if(verbose>0) print_4vector(A1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];


  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];

  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vth*vth*gthth + vph*vph*gphph );
  if(ut2<0.) 
    {
      ut2=0.;
     }
  ldouble ut=sqrtl(ut2);

  ulab[0]=ut;
  ulab[1]=vr*ut;
  ulab[2]=vth*ut;
  ulab[3]=vph*ut;
  
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
boost2_ff2zamo(ldouble A1[4],ldouble A2[4],ldouble *pp,ldouble gg[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble At[4]   ;

  int verbose=0;

  if(verbose>0) print_tensor(eup);
  if(verbose>0) print_4vector(A1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];


  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];

  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vth*vth*gthth + vph*vph*gphph );
  if(ut2<0.) 
    {
      ut2=0.;
     }
  ldouble ut=sqrtl(ut2);

  ulab[0]=ut;
  ulab[1]=vr*ut;
  ulab[2]=vth*ut;
  ulab[3]=vph*ut;
  
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
boost22_zamo2ff(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;
  //debug
  //  if(pp[2]<0.) verbose=1;

  if(verbose>0) print_tensor(eup);

  if(verbose) printf("sqrtg: %Le %Le %Le %Le\n",sqrtl(-gg[0][0]),sqrtl(gg[1][1]),sqrtl(gg[2][2]),sqrtl(gg[3][3]));

  if(verbose>0) print_tensor(T1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];


  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];

  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vth*vth*gthth + vph*vph*gphph );
  if(ut2<0.) 
    {
      ut2=0.;
    }

  ldouble ut=sqrtl(ut2);

  ulab[0]=ut;
  ulab[1]=vr*ut;
  ulab[2]=vth*ut;
  ulab[3]=vph*ut;
  
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
boost22_ff2zamo(ldouble T1[][4],ldouble T2[][4],ldouble *pp,ldouble gg[][5],ldouble eup[][4])
{ 
  int i,j,k,l;
  ldouble Tt[4][4];

  int verbose=0;
  
  if(verbose>0) print_tensor(eup);
  if(verbose>0) print_tensor(T1);

  //calculating the proper velocity of fluid as measured from ZAMO
  ldouble ulab[4],uzamo[4],vpr[3];

  ldouble vr=pp[2];
  ldouble vth=pp[3];
  ldouble vph=pp[4];

  ldouble gtt=gg[0][0];
  ldouble gtph=gg[0][3];
  ldouble grr=gg[1][1];
  ldouble gthth=gg[2][2];
  ldouble gphph=gg[3][3];

  ldouble ut2=-1./(gtt + 2.*vph*gtph + vr*vr*grr + vth*vth*gthth + vph*vph*gphph );
  if(ut2<0.) 
    {
      ut2=0.;
    }
  ldouble ut=sqrtl(ut2);

  ulab[0]=ut;
  ulab[1]=vr*ut;
  ulab[2]=vth*ut;
  ulab[3]=vph*ut;
  
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
//T^ij transfromation ZAMO -> lab
int
trans22_zamo2lab(ldouble T1[][4],ldouble T2[][4],ldouble elo[][4])
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
		  T2[i][j]+=elo[i][k]*elo[j][l]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation lab -> ZAMO
int
trans22_lab2zamo(ldouble T1[][4],ldouble T2[][4],ldouble eup[][4])
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
		  T2[i][j]+=eup[i][k]*eup[j][l]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation lab -> ZAMO
int
trans2_lab2zamo(ldouble *u1,ldouble *u2,ldouble eup[][4])
{
  int i,j,k;
  ldouble ut[4];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*eup[j][i];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation ZAMO -> lab
int
trans2_zamo2lab(ldouble *u1,ldouble *u2,ldouble elo[][4])
{
  int i,j,k;
  ldouble ut[4];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*elo[i][j];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation ortonormal to code coordinates
int
trans22_on2cc(ldouble T1[][4],ldouble T2[][4],ldouble tlo[][4])
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
		  T2[i][j]+=tlo[i][k]*tlo[j][l]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//T^ij transfromation code coords -> ortonormal
int
trans22_cc2on(ldouble T1[][4],ldouble T2[][4],ldouble tup[][4])
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
		  T2[i][j]+=tup[i][k]*tup[j][l]*Tt[k][l];
		}
	    }
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation cc -> on
int
trans2_cc2on(ldouble *u1,ldouble *u2,ldouble tup[][4])
{
  int i,j,k;
  ldouble ut[4];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*tup[j][i];
	}
    }

  return 0;
}

/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
//u^i transfromation on -> cc
int
trans2_on2cc(ldouble *u1,ldouble *u2,ldouble tlo[][4])
{
  int i,j,k;
  ldouble ut[4];

  for(i=0;i<4;i++)
    ut[i]=u1[i];

  for(i=0;i<4;i++)
    {
      u2[i]=0.;
      for(j=0;j<4;j++)
	{
	  u2[i]+=ut[j]*tlo[i][j];
	}
    }

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
    printf("%10.3Le %10.3Le %10.3Le %10.3Le\n",T[i][0],T[i][1],T[i][2],T[i][3]);
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
    printf("%10.3Le %10.3Le %10.3Le %10.3Le\n",T[i][0],T[i][1],T[i][2],T[i][3]);
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
  printf("%10.3Le %10.3Le %10.3Le %10.3Le\n",v[0],v[1],v[2],v[3]);
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
  printf("%10.7Le ",v[i]);
  printf("\n============\n");
  return 0;  
}

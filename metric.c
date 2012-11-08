//KORAL - metric.c
//metric-related routines
#include "ko.h"

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns metric determinant sqrt(-g)
ldouble
calc_gdet(ldouble *xx)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
 
#ifdef SPHERICAL
  return sqrtl(Power(x1,4)*Power(Sin(x2),2));
#endif

#ifdef CYLINDRICAL
  return x1;
#endif

#ifdef MINKOWSKI
  return 1.;
#endif

#ifdef KERR
 ldouble a=BHSPIN;
  return Sqrt(Power(Power(a,2) + 2*Power(x1,2) + 
       Power(a,2)*Cos(2*x2),2)*Power(Sin(x2),2))
    /2.;
#endif

#ifdef SCHWARZSCHILD
  return sqrtl(Power(x1,4)*Power(Sin(x2),2));
#endif  
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns D[gdet,x^idim]/gdet
ldouble
calc_dlgdet(ldouble *xx, int idim)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
 
#ifdef SPHERICAL
;if(idim==0) return  2/x1
;if(idim==1) return  Cot(x2)
;if(idim==2) return  0
;
#endif

#ifdef CYLINDRICAL
  if(idim==0) return 1./x1;
  if(idim==1) return 0.;
  if(idim==2) return 0.;
#endif

#ifdef MINKOWSKI
  return 0.;
#endif

#ifdef KERR
 ldouble a=BHSPIN;
;if(idim==0) return  (4*x1)/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(x1,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))
;if(idim==2) return  0
;
#endif

#ifdef SCHWARZSCHILD
;if(idim==0) return  2/x1
;if(idim==1) return  Cot(x2)
;if(idim==2) return  0
;
#endif  
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates base vectors and 1-forms of LNRF to transform lab <--> LNRF
int
calc_LNRFes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4])
{
  ldouble e2nu,e2psi,e2mu1,e2mu2,omega;
  ldouble gtt,gtph,gphph,grr,gthth;
  int i,j;

  gtt=g[0][0];
  gtph=g[0][3];
  gphph=g[3][3];
  grr=g[1][1];
  gthth=g[2][2];

  //Bardeen's 72 coefficients:
  e2nu=-gtt+gtph*gtph/gphph;
  e2psi=gphph;
  e2mu1=grr;
  e2mu2=gthth;
  omega=-gtph/gphph;

  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      {
	emuup[i][j]=0.;
	emulo[i][j]=0.;
      }

  emuup[0][0]=sqrtl(e2nu);
  emuup[1][1]=sqrtl(e2mu1);
  emuup[2][2]=sqrtl(e2mu2);
  emuup[0][3]=-omega*sqrtl(e2psi);
  emuup[3][3]=sqrtl(e2psi);

  emulo[3][0]=omega*1./sqrtl(e2nu);
  emulo[0][0]=1./sqrtl(e2nu);
  emulo[1][1]=1./sqrtl(e2mu1);
  emulo[2][2]=1./sqrtl(e2mu2);
  emulo[3][3]=1./sqrtl(e2psi);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//g_ij
int
calc_g(ldouble *xx, ldouble g[][5])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  g[3][4]=calc_gdet(xx);

#ifdef SCHWARZSCHILD
;g[0][0]= -1 + 2/x1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1/(1 - 2/x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2)
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(x1,2)*Power(Sin(x2),2)
;
#endif

#ifdef KERR
 ldouble a=BHSPIN;
;g[0][0]= -1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= (2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][0]= 0
;g[1][1]= (Power(x1,2) + Power(a,2)*Power(Cos(x2),2))/(Power(a,2) + (-2 + x1)*x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(Sin(x2),2)*(Power(a,2) + Power(x1,2) + (2*Power(a,2)*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;
#endif

#ifdef SPHERICAL
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2)
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
   ;g[3][3]= Power(x1,2)*Power(Sin(x2),2)
;
#endif


#ifdef CYLINDRICAL
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= 1
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= x1*x1 
;
#endif

#ifdef MINKOWSKI
;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= 1
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= 1
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= 1
;
#endif

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//g^ij
int
calc_G(ldouble *xx, ldouble G[][5])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

#ifdef SCHWARZSCHILD
;G[0][0]= x1/(2 - x1)
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= (-2 + x1)/x1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= Power(x1,-2)
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/Power(x1,2)
;

#endif

#ifdef KERR
 ldouble a=BHSPIN;
;G[0][0]= -((Power(a,4) + 2*Power(x1,4) + Power(a,2)*x1*(2 + 3*x1) + Power(a,2)*(Power(a,2) + (-2 + x1)*x1)*Cos(2*x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))))
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= (4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;G[1][0]= 0
;G[1][1]= (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= (4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= (2*((-2 + x1)*x1 + Power(a,2)*Power(Cos(x2),2))*Power(Csc(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;



#endif

#ifdef SPHERICAL
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= 1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= Power(x1,-2)
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/Power(x1,2)
;
#endif

#ifdef CYLINDRICAL
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= 1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(x1,-2)
;
#endif

#ifdef MINKOWSKI
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= 1
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= 1
;
#endif

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//Christofells: \Gamma^i_jk
int
calc_Krzysie(ldouble *xx, ldouble Krzys[][4][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

#ifdef SCHWARZSCHILD
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 1/((-2 + x1)*x1)
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 1/((-2 + x1)*x1)
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= (-2 + x1)/Power(x1,3)
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 1/(2*x1 - Power(x1,2))
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 2 - x1
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -((-2 + x1)*Power(Sin(x2),2))
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 1/x1
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 1/x1
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -(Cos(x2)*Sin(x2))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 1/x1
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= Cot(x2)
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 1/x1
;Krzys[3][3][2]= Cot(x2)
;Krzys[3][3][3]= 0
;


#endif

#ifdef KERR
 ldouble a=BHSPIN;
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= (-2*(Power(a,2) + Power(x1,2))*(Power(a,2) - 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][0][2]= (-4*Power(a,2)*x1*Sin(2*x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= (-2*(Power(a,2) + Power(x1,2))*(Power(a,2) - 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= (-2*a*(Power(a,4) - 3*Power(a,2)*Power(x1,2) - 6*Power(x1,4) + Power(a,2)*(Power(a,2) - Power(x1,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][2][0]= (-4*Power(a,2)*x1*Sin(2*x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= (-8*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= (-2*a*(Power(a,4) - 3*Power(a,2)*Power(x1,2) - 6*Power(x1,4) + Power(a,2)*(Power(a,2) - Power(x1,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][3][2]= (-8*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= -(((Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= -((a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= ((Power(a,2) - x1)*x1 - Power(a,2)*(-1 + x1)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][2]= -((x1*(Power(a,2) + (-2 + x1)*x1))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= -((a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + (-2 + x1)*x1)*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][0][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= -((a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= (Power(a,2)*Cos(x2)*Sin(x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][1][2]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= -((a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -((Cos(x2)*Sin(x2)*(2*Power(a,2)*Power(x1,2)*(Power(a,2) + Power(x1,2))*Power(Cos(x2),2) + Power(a,4)*(Power(a,2) + Power(x1,2))*Power(Cos(x2),4) + x1*(Power(a,2)*Power(x1,3) + Power(x1,5) + 4*Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + 2*Power(a,4)*Power(Sin(x2),4) + Power(a,4)*Power(Sin(2*x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= (-4*a*Power(x1,2) + 4*Power(a,3)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][0][2]= (8*a*x1*Cot(x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= (-4*a*Power(x1,2) + 4*Power(a,3)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= (4*((-2 + x1)*Power(x1,4) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(a,2)*Power(Cos(x2),2)*(2*(-1 + x1)*Power(x1,2) + Power(a,2)*Power(Sin(x2),2))))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][2][0]= (8*a*x1*Cot(x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= ((3*Power(a,4) + 8*Power(a,2)*x1 + 8*Power(a,2)*Power(x1,2) + 8*Power(x1,4) + 4*Power(a,2)*(Power(a,2) + 2*(-1 + x1)*x1)*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= (4*((-2 + x1)*Power(x1,4) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(a,2)*Power(Cos(x2),2)*(2*(-1 + x1)*Power(x1,2) + Power(a,2)*Power(Sin(x2),2))))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][2]= ((3*Power(a,4) + 8*Power(a,2)*x1 + 8*Power(a,2)*Power(x1,2) + 8*Power(x1,4) + 4*Power(a,2)*(Power(a,2) + 2*(-1 + x1)*x1)*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][3]= 0
;
#endif


#ifdef SPHERICAL
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 0
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= -x1
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(x1*Power(Sin(x2),2))
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 1/x1
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 1/x1
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -(Cos(x2)*Sin(x2))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 1/x1
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= Cot(x2)
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 1/x1
;Krzys[3][3][2]= Cot(x2)
;Krzys[3][3][3]= 0
;
#endif

#ifdef CYLINDRICAL
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 0
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 0
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -x1
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 0
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 0
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= 0
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 1/x1
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= 0
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 1/x1
;Krzys[3][3][2]= 0
;Krzys[3][3][3]= 0
;
#endif

#ifdef MINKOWSKI
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= 0
;Krzys[0][0][2]= 0
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= 0
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= 0
;Krzys[0][2][0]= 0
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= 0
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= 0
;Krzys[0][3][2]= 0
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= 0
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= 0
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= 0
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 0
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= 0
;Krzys[2][0][0]= 0
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= 0
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= 0
;Krzys[2][1][2]= 0
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= 0
;Krzys[2][2][2]= 0
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= 0
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= 0
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= 0
;Krzys[3][0][2]= 0
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= 0
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= 0
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= 0
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= 0
;Krzys[3][3][2]= 0
;Krzys[3][3][3]= 0
;
#endif

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************

int
print_Krzysie(ldouble g[][4][4])
{
  int i,j,k;
  for(i=0;i<4;i++)
    {
      for(j=0;j<4;j++)
	{
	  for(k=0;k<4;k++)
	    {
	      printf("%10Lf ",g[i][j][k]);
	    }
	  printf("\n");
	}
      printf("\n");
    }
  printf("\n");
  return -1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************

int
print_g(ldouble g[][5])
{
  int i,j,k;
  for(j=0;j<4;j++)
    {
      for(k=0;k<4;k++)
	{
	  printf("%10Lf ",g[j][k]);
	}
      printf("\n");
    }
  printf("\n");
  return -1;
}

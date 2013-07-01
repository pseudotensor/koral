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
  return calc_gdet_arb(xx,MYCOORDS);
}
 
ldouble
calc_gdet_arb(ldouble *xx,int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
 
if(coords==SPHCOORDS) {
  return sqrt(Power(x1,4)*Power(Sin(x2),2));
 }

if(coords==CYLCOORDS) {
  return x1;
 }

if(coords==MINKCOORDS) {
  return 1.;
 }

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
  return Sqrt(Power(Power(a,2) + 2*Power(x1,2) + 
       Power(a,2)*Cos(2*x2),2)*Power(Sin(x2),2))
    /2.;
 }

if(coords==KSCOORDS) {
  ldouble a=BHSPIN;
return Sqrt(Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)*
	    Power(Sin(x2),2));
 }

/*
if(coords==SCHWCOORDS) {
  return sqrt(Power(x1,4)*Power(Sin(x2),2));
 } 
*/

if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif
  return Sqrt(Power(exp(1.0),2*x1)*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)*Power(Sin(x2),2))/2.;
 } 

if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif
  return Sqrt(Power(exp(1.0),2*x1)*Power(Power(exp(1.0),x1) + R0,2));
 } 

 return 0.;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns D[gdet,x^idim]/gdet
ldouble
calc_dlgdet(ldouble *xx, int idim)
{
  return calc_dlgdet_arb(xx,idim,MYCOORDS);
}

ldouble
calc_dlgdet_arb(ldouble *xx, int idim,int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
 
if(coords==SPHCOORDS) {
;if(idim==0) return  2/x1
;if(idim==1) return  Cot(x2)
;if(idim==2) return  0
;
}

if(coords==CYLCOORDS) {
  if(idim==0) return 1./x1;
  if(idim==1) return 0.;
  if(idim==2) return 0.;
}

if(coords==MINKCOORDS) {
  return 0.;
}

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;if(idim==0) return  (4*x1)/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(x1,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))
;if(idim==2) return  0
;
}

if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;if(idim==0) return  (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(x1,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(2.*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;if(idim==2) return  0
;
}
/*
if(coords==SCHWCOORDS) {
;if(idim==0) return  2/x1
;if(idim==1) return  Cot(x2)
;if(idim==2) return  0
;
}  */

if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif
;if(idim==0) return  (Power(a,2) + 6*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),x1)*R0 + 2*Power(R0,2) + Power(a,2)*Cos(2*x2))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))
;if(idim==1) return  ((-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + 3*Power(a,2)*Cos(2*x2))*Cot(x2))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))
;if(idim==2) return  0
;
 }

if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif
;if(idim==0) return  (2*Power(exp(1.0),x1) + R0)/(Power(exp(1.0),x1) + R0)
;if(idim==1) return  0
;if(idim==2) return  0
;
 }

 return 0;

}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//g_ij
int
calc_g(ldouble *xx, ldouble g[][5])
{
  calc_g_arb(xx,g,MYCOORDS);
  return 0;
}

int
calc_g_arb(ldouble *xx, ldouble g[][5],int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  g[3][4]=calc_gdet(xx);

if(coords==MCYL1COORDS) {
#if(MYCOORDS==MCYL1COORDS)
  ldouble R0;
  R0=MKS1R0;

;g[0][0]= -1
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= 0
;g[1][0]= 0
;g[1][1]= Power(exp(1.0),2*x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= 1
;g[2][3]= 0
;g[3][0]= 0
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(Power(exp(1.0),x1) + R0,2)
;
#endif
}

if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif
;g[0][0]= -((Power(a,2) + 2*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0) + Power(a,2)*Cos(2*x2))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2)))
;g[0][1]= (4*Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[0][2]= 0
;g[0][3]= (-4*a*(Power(exp(1.0),x1) + R0)*Power(Sin(x2),2))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[1][0]= (4*Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[1][1]= (4*Power(exp(1.0),2*x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;g[1][2]= 0
;g[1][3]= (-4*a*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-4*a*(Power(exp(1.0),x1) + R0)*Power(Sin(x2),2))/(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))
;g[3][1]= (-4*a*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;g[3][2]= 0
;g[3][3]= (4*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*((Power(exp(1.0),x1) + R0)*(Power(Power(exp(1.0),x1) + R0,3) + Power(a,2)*(2 + Power(exp(1.0),x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;
}


if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;g[0][0]= -1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][1]= 0
;g[0][2]= 0
;g[0][3]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][0]= 0
;g[1][1]= (Power(x1,2) + Power(a,2)*Power(Cos(x2),2))/(Power(a,2) + (-2 + x1)*x1)
;g[1][2]= 0
;g[1][3]= 0
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[3][1]= 0
;g[3][2]= 0
;g[3][3]= Power(Sin(x2),2)*(Power(a,2) + Power(x1,2) + (2*Power(a,2)*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;
}

/*
if(coords==SCHWCOORDS) {
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
}
*/

if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;g[0][0]= -1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][1]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[0][2]= 0
;g[0][3]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][0]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][1]= 1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[1][2]= 0
;g[1][3]= -(a*(1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))*Power(Sin(x2),2))
;g[2][0]= 0
;g[2][1]= 0
;g[2][2]= Power(x1,2) + Power(a,2)*Power(Cos(x2),2)
;g[2][3]= 0
;g[3][0]= (-2*a*x1*Power(Sin(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;g[3][1]= -(a*(1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))*Power(Sin(x2),2))
;g[3][2]= 0
;g[3][3]= Power(Sin(x2),2)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2) + Power(a,2)*(1 + (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))*Power(Sin(x2),2))
;
}

if(coords==SPHCOORDS) {
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
}


if(coords==CYLCOORDS) {
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
}

if(coords==MINKCOORDS) {
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
}

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//g^ij
int
calc_G(ldouble *xx, ldouble G[][5])
{
  calc_G_arb(xx,G,MYCOORDS);
  return 0;
}

int
calc_G_arb(ldouble *xx, ldouble G[][5],int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];

  if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif
;G[0][0]= -1
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= 0
;G[1][1]= Power(exp(1.0),-2*x1)
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= Power(Power(exp(1.0),x1) + R0,-2)
;

  }

  if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif
;G[0][0]= -(((Power(exp(1.0),x1) + R0)*(2 + Power(exp(1.0),x1) + R0) + Power(a,2)*Power(Cos(x2),2))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[0][1]= (2*(Power(exp(1.0),x1) + R0))/(Power(exp(1.0),x1)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= (2*(Power(exp(1.0),x1) + R0))/(Power(exp(1.0),x1)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[1][1]= (Power(a,2) + (-2 + Power(exp(1.0),x1) + R0)*(Power(exp(1.0),x1) + R0))/(Power(exp(1.0),2*x1)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[1][2]= 0
;G[1][3]= a/(Power(exp(1.0),x1)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= a/(Power(exp(1.0),x1)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;


  }

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;G[0][0]= -((Power(a,4) + 2*Power(x1,4) + Power(a,2)*x1*(2 + 3*x1) + Power(a,2)*(Power(a,2) + (-2 + x1)*x1)*Cos(2*x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2))))
;G[0][1]= 0
;G[0][2]= 0
;G[0][3]= (-4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;G[1][0]= 0
;G[1][1]= (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][2]= 0
;G[1][3]= 0
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= (-4*a*x1)/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;G[3][1]= 0
;G[3][2]= 0
;G[3][3]= (2*((-2 + x1)*x1 + Power(a,2)*Power(Cos(x2),2))*Power(Csc(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;
}

/*
if(coords==SCHWCOORDS) {
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

}
*/

if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;G[0][0]= -((x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;G[0][1]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[0][2]= 0
;G[0][3]= 0
;G[1][0]= (2*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][1]= (Power(a,2) + (-2 + x1)*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[1][2]= 0
;G[1][3]= a/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][0]= 0
;G[2][1]= 0
;G[2][2]= 1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[2][3]= 0
;G[3][0]= 0
;G[3][1]= a/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;G[3][2]= 0
;G[3][3]= Power(Csc(x2),2)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;
}


if(coords==SPHCOORDS) {
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
}

if(coords==CYLCOORDS) {
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
}

if(coords==MINKCOORDS) {
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
}

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//Christofells: \Gamma^i_jk
int
calc_Krzysie(ldouble *xx, ldouble Krzys[][4][4])
{
  calc_Krzysie_arb(xx,Krzys,MYCOORDS);
  return 0;
}

int
calc_Krzysie_at_center(int ix,int iy,int iz, ldouble Krzys[][4][4])
{
#if(MODYFIKUJKRZYSIE==0)
  ldouble xx[4];

  xx[0]=0.;
  xx[1]=get_x(ix,0);
  xx[2]=get_x(iy,1);
  xx[3]=get_x(iz,2);

  calc_Krzysie_arb(xx,Krzys,MYCOORDS);
#else
  ldouble xx[4];
  ldouble Krzys_org[4][4][4];

  //analytical at center

  xx[0]=0.;
  xx[1]=get_x(ix,0);
  xx[2]=get_x(iy,1);
  xx[3]=get_x(iz,2);

  calc_Krzysie_arb(xx,Krzys_org,MYCOORDS);
  calc_Krzysie_arb(xx,Krzys,MYCOORDS);

  //modifying \Gamma ^mu_mu_k

  int kappa,mu;
  ldouble Sk,Wk[4],dS[4],gdet[3],D[4],dxk,Ck;
  for(kappa=1;kappa<=3;kappa++)
    {
      Ck=Krzys_org[0][kappa][0]
	+ Krzys_org[1][kappa][1]
	+ Krzys_org[2][kappa][2]
	+ Krzys_org[3][kappa][3];
 
      Sk=1.e-300 + fabs(Krzys_org[0][kappa][0])
	+ fabs(Krzys_org[1][kappa][1])
	+ fabs(Krzys_org[2][kappa][2])
	+ fabs(Krzys_org[3][kappa][3]);

      for(mu=0;mu<4;mu++)
	Wk[mu]=fabs(Krzys_org[mu][kappa][mu])/Sk;

     
      //center
      xx[0]=0.;
      xx[1]=get_x(ix,0);
      xx[2]=get_x(iy,1);
      xx[3]=get_x(iz,2);
      gdet[0]=calc_gdet_arb(xx,MYCOORDS);

      //upper face
      xx[0]=0.;
      xx[1]=get_x(ix,0);
      xx[2]=get_x(iy,1);
      xx[3]=get_x(iz,2);
      if(kappa==1) xx[1]=get_xb(ix+1,0);
      if(kappa==2) xx[2]=get_xb(iy+1,1);
      if(kappa==3) xx[3]=get_xb(iz+1,2);
      gdet[1]=calc_gdet_arb(xx,MYCOORDS);

      //lower face
      xx[0]=0.;
      xx[1]=get_x(ix,0);
      xx[2]=get_x(iy,1);
      xx[3]=get_x(iz,2);
      if(kappa==1) xx[1]=get_xb(ix,0);
      if(kappa==2) xx[2]=get_xb(iy,1);
      if(kappa==3) xx[3]=get_xb(iz,2);
      gdet[2]=calc_gdet_arb(xx,MYCOORDS);

      //numerical differencing
      if(kappa==1) dxk=get_size_x(ix,0);
      if(kappa==2) dxk=get_size_x(iy,1);
      if(kappa==3) dxk=get_size_x(iz,2);
      D[kappa]=(gdet[1]-gdet[2])/(dxk*gdet[0]);

      //correcting Krzysie
      for(mu=0;mu<4;mu++)
	{
	  Krzys[mu][kappa][mu]+=(D[kappa]-Ck)*Wk[mu];
	  Krzys[mu][mu][kappa]=Krzys[mu][kappa][mu];

	  /*
	  printf("%f %f > %d %d >  %e %e %e | %e -> %e\n",get_x(ix,0),get_x(iy,1),
		 kappa,mu,
		 D[kappa],Ck,Wk[mu],Krzys_org[mu][kappa][mu],Krzys[mu][kappa][mu]); getchar();
	  */
	  
	}
    }

#endif
  return 0;
}

int
calc_Krzysie_arb(ldouble *xx, ldouble Krzys[][4][4],int coords)
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];


  if(coords==MCYL1COORDS) {
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif
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
;Krzys[1][1][1]= 1
;Krzys[1][1][2]= 0
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= 0
;Krzys[1][2][2]= 0
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= 0
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -1 - R0/Power(exp(1.0),x1)
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
;Krzys[3][1][3]= Power(exp(1.0),x1)/(Power(exp(1.0),x1) + R0)
;Krzys[3][2][0]= 0
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= 0
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= Power(exp(1.0),x1)/(Power(exp(1.0),x1) + R0)
;Krzys[3][3][2]= 0
;Krzys[3][3][3]= 0
;
  }

  if(coords==MKS1COORDS) {
  ldouble a=BHSPIN;
  ldouble R0;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif
;Krzys[0][0][0]= (8*(Power(exp(1.0),x1) + R0)*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][0][1]= (4*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][0][2]= (-4*Power(a,2)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][0][3]= (8*a*(Power(exp(1.0),x1) + R0)*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][1][0]= (4*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][1][1]= (2*Power(exp(1.0),2*x1)*(-3*Power(a,4) - 4*Power(a,2)*Power(exp(1.0),x1) + 8*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),4*x1) - 4*Power(a,2)*R0 + 24*Power(exp(1.0),2*x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 24*Power(exp(1.0),x1)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 8*Power(R0,3) + 32*Power(exp(1.0),x1)*Power(R0,3) + 8*Power(R0,4) - 4*Power(a,2)*(Power(a,2) + Power(exp(1.0),x1) + R0)*Cos(2*x2) - Power(a,4)*Cos(4*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][1][2]= (-4*Power(a,2)*Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][1][3]= (-4*a*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][2][0]= (-4*Power(a,2)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][1]= (-4*Power(a,2)*Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][2]= (-2*Power(Power(exp(1.0),x1) + R0,2))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[0][2][3]= (8*Power(a,3)*(Power(exp(1.0),x1) + R0)*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][0]= (8*a*(Power(exp(1.0),x1) + R0)*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][3][1]= (-4*a*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[0][3][2]= (8*Power(a,3)*(Power(exp(1.0),x1) + R0)*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][3]= (-2*(Power(exp(1.0),x1) + R0)*(Power(a,4) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][0][0]= (4*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/(Power(exp(1.0),x1)*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][0][1]= (-4*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2))*(2*(Power(exp(1.0),x1) + R0) - Power(a,2)*Power(Sin(x2),2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (4*a*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/(Power(exp(1.0),x1)*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][1][0]= (-4*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2))*(2*(Power(exp(1.0),x1) + R0) - Power(a,2)*Power(Sin(x2),2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][1][1]= (10*Power(a,6) + 8*Power(a,4)*Power(exp(1.0),x1) + 32*Power(a,2)*Power(exp(1.0),2*x1) + 36*Power(a,4)*Power(exp(1.0),2*x1) + 16*Power(a,2)*Power(exp(1.0),3*x1) - 64*Power(exp(1.0),4*x1) + 48*Power(a,2)*Power(exp(1.0),4*x1) - 32*Power(exp(1.0),5*x1) + 32*Power(exp(1.0),6*x1) + 32*Power(a,2)*Power(exp(1.0),x1)*R0 + 72*Power(a,4)*Power(exp(1.0),x1)*R0 + 32*Power(a,2)*Power(exp(1.0),2*x1)*R0 - 192*Power(exp(1.0),3*x1)*R0 + 192*Power(a,2)*Power(exp(1.0),3*x1)*R0 - 128*Power(exp(1.0),4*x1)*R0 + 192*Power(exp(1.0),5*x1)*R0 + 36*Power(a,4)*Power(R0,2) + 16*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) - 192*Power(exp(1.0),2*x1)*Power(R0,2) + 288*Power(a,2)*Power(exp(1.0),2*x1)*Power(R0,2) - 192*Power(exp(1.0),3*x1)*Power(R0,2) + 480*Power(exp(1.0),4*x1)*Power(R0,2) - 64*Power(exp(1.0),x1)*Power(R0,3) + 192*Power(a,2)*Power(exp(1.0),x1)*Power(R0,3) - 128*Power(exp(1.0),2*x1)*Power(R0,3) + 640*Power(exp(1.0),3*x1)*Power(R0,3) + 48*Power(a,2)*Power(R0,4) - 32*Power(exp(1.0),x1)*Power(R0,4) + 480*Power(exp(1.0),2*x1)*Power(R0,4) + 192*Power(exp(1.0),x1)*Power(R0,5) + 32*Power(R0,6) + Power(a,2)*(15*Power(a,4) + 16*Power(a,2)*(3*Power(exp(1.0),2*x1) + 3*Power(R0,2) + Power(exp(1.0),x1)*(1 + 6*R0)) + 16*(3*Power(exp(1.0),4*x1) + 3*Power(R0,4) + Power(exp(1.0),3*x1)*(-1 + 12*R0) + 2*Power(exp(1.0),2*x1)*(1 - R0 + 9*Power(R0,2)) + Power(exp(1.0),x1)*R0*(2 - R0 + 12*Power(R0,2))))*Cos(2*x2) + 2*Power(a,4)*(3*Power(a,2) + 6*Power(exp(1.0),2*x1) + 6*Power(R0,2) + 4*Power(exp(1.0),x1)*(1 + 3*R0))*Cos(4*x2) + Power(a,6)*Cos(6*x2))/(4.*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][1][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][3]= (a*(Power(a,4) - 8*Power(a,2)*Power(exp(1.0),x1) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) - 8*Power(a,2)*R0 + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),x1)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 16*Power(R0,3) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(-2 + Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][2]= -(((Power(exp(1.0),x1) + R0)*(Power(a,2) + (-2 + Power(exp(1.0),x1) + R0)*(Power(exp(1.0),x1) + R0)))/(Power(exp(1.0),x1)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (4*a*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/(Power(exp(1.0),x1)*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[1][3][1]= (a*(Power(a,4) - 8*Power(a,2)*Power(exp(1.0),x1) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) - 8*Power(a,2)*R0 + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),x1)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 16*Power(R0,3) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(-2 + Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*(Power(a,4) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/(Power(exp(1.0),x1)*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)))
;Krzys[2][0][0]= (-8*Power(a,2)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][0][1]= (-8*Power(a,2)*Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (8*a*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(Power(exp(1.0),x1) + R0,2))*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][0]= (-8*Power(a,2)*Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][1]= (-8*Power(a,2)*Power(exp(1.0),2*x1)*(Power(exp(1.0),x1) + R0)*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][1][2]= (Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= (a*Power(exp(1.0),x1)*(3*Power(a,4) + 16*Power(a,2)*Power(exp(1.0),x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),4*x1) + 16*Power(a,2)*R0 + 16*Power(a,2)*Power(exp(1.0),x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),x1)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 16*Power(R0,3) + 32*Power(exp(1.0),x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Sin(2*x2))/(2.*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= (Power(exp(1.0),x1)*(Power(exp(1.0),x1) + R0))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (8*a*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(Power(exp(1.0),x1) + R0,2))*Sin(2*x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[2][3][1]= (a*Power(exp(1.0),x1)*(3*Power(a,4) + 16*Power(a,2)*Power(exp(1.0),x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 16*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),4*x1) + 16*Power(a,2)*R0 + 16*Power(a,2)*Power(exp(1.0),x1)*R0 + 48*Power(exp(1.0),2*x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),x1)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 16*Power(R0,3) + 32*Power(exp(1.0),x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Sin(2*x2))/(2.*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= (4*Sin(x2)*(-(Cos(x2)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*((Power(exp(1.0),x1) + R0)*(Power(Power(exp(1.0),x1) + R0,3) + Power(a,2)*(2 + Power(exp(1.0),x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*Cos(x2)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2) + Power(a,2)*Cos(x2)*((Power(exp(1.0),x1) + R0)*(Power(Power(exp(1.0),x1) + R0,3) + Power(a,2)*(2 + Power(exp(1.0),x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2) - 2*Power(a,2)*(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*((Power(exp(1.0),x1) + R0)*(Power(Power(exp(1.0),x1) + R0,3) + Power(a,2)*(2 + Power(exp(1.0),x1) + R0)) + Power(a,2)*(Power(a,2) + Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(-1 + R0) + (-2 + R0)*R0)*Power(Cos(x2),2))*Sin(x2)*Sin(2*x2)))/((Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2))*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;Krzys[3][0][0]= (-4*a*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][0][1]= (4*a*Power(exp(1.0),x1)*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][0][2]= (-8*a*(Power(exp(1.0),x1) + R0)*Cot(x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][0][3]= (4*Power(a,2)*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][1][0]= (4*a*Power(exp(1.0),x1)*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][1][1]= (4*a*Power(exp(1.0),2*x1)*(-Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) - Power(a,2)*Cos(2*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][1][2]= (-4*a*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][1][3]= (Power(exp(1.0),x1)*(Power(a,4) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][2][0]= (-8*a*(Power(exp(1.0),x1) + R0)*Cot(x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][1]= (-4*a*Power(exp(1.0),x1)*(Power(exp(1.0),2*x1) + 2*Power(exp(1.0),x1)*(1 + R0) + R0*(2 + R0) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][2]= -((a*(Power(exp(1.0),x1) + R0))/(Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[3][2][3]= ((3*Power(a,4) + 8*Power(a,2)*Power(exp(1.0),x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),4*x1) + 8*Power(a,2)*R0 + 16*Power(a,2)*Power(exp(1.0),x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 32*Power(exp(1.0),x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*(Power(exp(1.0),2*x1) + (-1 + R0)*R0 + Power(exp(1.0),x1)*(-1 + 2*R0)))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][0]= (4*Power(a,2)*(Power(a,2) - 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][3][1]= (Power(exp(1.0),x1)*(Power(a,4) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2)))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3)
;Krzys[3][3][2]= ((3*Power(a,4) + 8*Power(a,2)*Power(exp(1.0),x1) + 8*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(exp(1.0),4*x1) + 8*Power(a,2)*R0 + 16*Power(a,2)*Power(exp(1.0),x1)*R0 + 32*Power(exp(1.0),3*x1)*R0 + 8*Power(a,2)*Power(R0,2) + 48*Power(exp(1.0),2*x1)*Power(R0,2) + 32*Power(exp(1.0),x1)*Power(R0,3) + 8*Power(R0,4) + 4*Power(a,2)*(Power(a,2) + 2*(Power(exp(1.0),2*x1) + (-1 + R0)*R0 + Power(exp(1.0),x1)*(-1 + 2*R0)))*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][3]= -((a*(Power(a,4) + 3*Power(a,4)*Power(exp(1.0),x1) - 4*Power(a,2)*Power(exp(1.0),2*x1) + 8*Power(a,2)*Power(exp(1.0),3*x1) + 8*Power(exp(1.0),5*x1) + 3*Power(a,4)*R0 - 8*Power(a,2)*Power(exp(1.0),x1)*R0 + 24*Power(a,2)*Power(exp(1.0),2*x1)*R0 + 40*Power(exp(1.0),4*x1)*R0 - 4*Power(a,2)*Power(R0,2) + 24*Power(a,2)*Power(exp(1.0),x1)*Power(R0,2) + 80*Power(exp(1.0),3*x1)*Power(R0,2) + 8*Power(a,2)*Power(R0,3) + 80*Power(exp(1.0),2*x1)*Power(R0,3) + 40*Power(exp(1.0),x1)*Power(R0,4) + 8*Power(R0,5) + 4*Power(a,2)*(Power(exp(1.0),x1) + R0)*(Power(a,2) + Power(exp(1.0),x1) + 2*Power(exp(1.0),2*x1) + R0 + 4*Power(exp(1.0),x1)*R0 + 2*Power(R0,2))*Cos(2*x2) + Power(a,4)*(-1 + Power(exp(1.0),x1) + R0)*Cos(4*x2))*Power(Sin(x2),2))/Power(Power(a,2) + 2*Power(Power(exp(1.0),x1) + R0,2) + Power(a,2)*Cos(2*x2),3))
;
}
  /*
if(coords==SCHWCOORDS) {
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


}
  */
if(coords==KSCOORDS) {
 ldouble a=BHSPIN;
;Krzys[0][0][0]= (2*x1*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][0][1]= ((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][0][2]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][0][3]= (2*a*x1*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][1][0]= ((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][1][1]= (2*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1 + Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][1][2]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][1][3]= (a*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][2][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][2][1]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][2][2]= (-2*Power(x1,2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[0][2][3]= (2*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][3][0]= (2*a*x1*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][3][1]= (a*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[0][3][2]= (2*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[0][3][3]= (-2*x1*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][0][0]= -(((Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][1]= -(((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(2*x1 - Power(a,2)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][1][0]= -(((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(2*x1 - Power(a,2)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][1][1]= -(((Power(x1,2) - Power(a,2)*Power(Cos(x2),2))*(x1*(2 + x1) + Power(a,2)*Cos(2*x2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][1][2]= -((Power(a,2)*Sin(2*x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;Krzys[1][1][3]= (a*Power(Sin(x2),2)*(Power(a,4)*x1*Power(Cos(x2),4) + Power(x1,2)*(2*x1 + Power(x1,3) - Power(a,2)*Power(Sin(x2),2)) + Power(a,2)*Power(Cos(x2),2)*(2*x1*(-1 + Power(x1,2)) + Power(a,2)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Sin(2*x2))/(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))
;Krzys[1][2][2]= -((x1*(Power(a,2) + (-2 + x1)*x1))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][3][1]= (a*Power(Sin(x2),2)*(Power(a,4)*x1*Power(Cos(x2),4) + Power(x1,2)*(2*x1 + Power(x1,3) - Power(a,2)*Power(Sin(x2),2)) + Power(a,2)*Power(Cos(x2),2)*(2*x1*(-1 + Power(x1,2)) + Power(a,2)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + (-2 + x1)*x1)*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][0][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][1]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][1]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][2]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= (a*Cos(x2)*Sin(x2)*(Power(x1,3)*(2 + x1) + 2*Power(a,2)*x1*(1 + x1)*Power(Cos(x2),2) + Power(a,4)*Power(Cos(x2),4) + 2*Power(a,2)*x1*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][3][1]= (a*Cos(x2)*Sin(x2)*(Power(x1,3)*(2 + x1) + 2*Power(a,2)*x1*(1 + x1)*Power(Cos(x2),2) + Power(a,4)*Power(Cos(x2),4) + 2*Power(a,2)*x1*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -((Cos(x2)*Sin(x2)*(Power(a,6)*Power(Cos(x2),6) + Power(Cos(x2),4)*(3*Power(a,4)*Power(x1,2) + Power(a,6)*Power(Sin(x2),2)) + Power(Cos(x2),2)*(3*Power(a,2)*Power(x1,4) + 2*Power(a,4)*Power(x1,2)*Power(Sin(x2),2)) + x1*(Power(x1,5) + Power(a,2)*Power(x1,2)*(4 + x1)*Power(Sin(x2),2) + 2*Power(a,4)*Power(Sin(x2),4) + Power(a,4)*Power(Sin(2*x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[3][0][0]= (a*Power(x1,2) - Power(a,3)*Power(Cos(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][0][1]= (a*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][0][2]= (-2*a*x1*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][0][3]= (Power(a,2)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][1][0]= (a*(Power(x1,2) - Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][1][1]= (a*Power(x1,2) - Power(a,3)*Power(Cos(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][1][2]= -((a*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2))
;Krzys[3][1][3]= (Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][2][0]= (-2*a*x1*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][2][1]= -((a*(x1*(2 + x1) + Power(a,2)*Power(Cos(x2),2))*Cot(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2))
;Krzys[3][2][2]= -((a*x1)/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[3][2][3]= ((Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)*Cot(x2))/4. + Power(a,2)*x1*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][3][0]= (Power(a,2)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][3][1]= (Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[3][3][2]= ((Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)*Cot(x2))/4. + Power(a,2)*x1*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),2)
;Krzys[3][3][3]= -((a*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;
}

if(coords==KERRCOORDS) {
 ldouble a=BHSPIN;
;Krzys[0][0][0]= 0
;Krzys[0][0][1]= (-2*(Power(a,2) + Power(x1,2))*(Power(a,2) - 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][0][2]= (-4*Power(a,2)*x1*Sin(2*x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][0][3]= 0
;Krzys[0][1][0]= (-2*(Power(a,2) + Power(x1,2))*(Power(a,2) - 2*Power(x1,2) + Power(a,2)*Cos(2*x2)))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][1][1]= 0
;Krzys[0][1][2]= 0
;Krzys[0][1][3]= (2*a*(Power(a,4) - 3*Power(a,2)*Power(x1,2) - 6*Power(x1,4) + Power(a,2)*(Power(a,2) - Power(x1,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][2][0]= (-4*Power(a,2)*x1*Sin(2*x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][2][1]= 0
;Krzys[0][2][2]= 0
;Krzys[0][2][3]= (8*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][0]= 0
;Krzys[0][3][1]= (2*a*(Power(a,4) - 3*Power(a,2)*Power(x1,2) - 6*Power(x1,4) + Power(a,2)*(Power(a,2) - Power(x1,2))*Cos(2*x2))*Power(Sin(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[0][3][2]= (8*Power(a,3)*x1*Cos(x2)*Power(Sin(x2),3))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[0][3][3]= 0
;Krzys[1][0][0]= -(((Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[1][0][1]= 0
;Krzys[1][0][2]= 0
;Krzys[1][0][3]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][1][0]= 0
;Krzys[1][1][1]= ((Power(a,2) - x1)*x1 - Power(a,2)*(-1 + x1)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][1][3]= 0
;Krzys[1][2][0]= 0
;Krzys[1][2][1]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][2]= -((x1*(Power(a,2) + (-2 + x1)*x1))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[1][2][3]= 0
;Krzys[1][3][0]= (a*(Power(a,2) + (-2 + x1)*x1)*(-Power(x1,2) + Power(a,2)*Power(Cos(x2),2))*Power(Sin(x2),2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[1][3][1]= 0
;Krzys[1][3][2]= 0
;Krzys[1][3][3]= -(((Power(a,2) + (-2 + x1)*x1)*Power(Sin(x2),2)*(Power(x1,5) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(Cos(x2),2)*(2*Power(a,2)*Power(x1,3) + Power(a,4)*Power(Sin(x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[2][0][0]= (-2*Power(a,2)*x1*Cos(x2)*Sin(x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][0][1]= 0
;Krzys[2][0][2]= 0
;Krzys[2][0][3]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][1][0]= 0
;Krzys[2][1][1]= (Power(a,2)*Cos(x2)*Sin(x2))/((Power(a,2) + (-2 + x1)*x1)*(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][1][2]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][1][3]= 0
;Krzys[2][2][0]= 0
;Krzys[2][2][1]= x1/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2))
;Krzys[2][2][2]= -((Power(a,2)*Cos(x2)*Sin(x2))/(Power(x1,2) + Power(a,2)*Power(Cos(x2),2)))
;Krzys[2][2][3]= 0
;Krzys[2][3][0]= (a*x1*(Power(a,2) + Power(x1,2))*Sin(2*x2))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3)
;Krzys[2][3][1]= 0
;Krzys[2][3][2]= 0
;Krzys[2][3][3]= -((Cos(x2)*Sin(x2)*(2*Power(a,2)*Power(x1,2)*(Power(a,2) + Power(x1,2))*Power(Cos(x2),2) + Power(a,4)*(Power(a,2) + Power(x1,2))*Power(Cos(x2),4) + x1*(Power(a,2)*Power(x1,3) + Power(x1,5) + 4*Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + 2*Power(a,4)*Power(Sin(x2),4) + Power(a,4)*Power(Sin(2*x2),2))))/Power(Power(x1,2) + Power(a,2)*Power(Cos(x2),2),3))
;Krzys[3][0][0]= 0
;Krzys[3][0][1]= (4*a*Power(x1,2) - 4*Power(a,3)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][0][2]= (-8*a*x1*Cot(x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][0][3]= 0
;Krzys[3][1][0]= (4*a*Power(x1,2) - 4*Power(a,3)*Power(Cos(x2),2))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][1][1]= 0
;Krzys[3][1][2]= 0
;Krzys[3][1][3]= (4*((-2 + x1)*Power(x1,4) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(a,2)*Power(Cos(x2),2)*(2*(-1 + x1)*Power(x1,2) + Power(a,2)*Power(Sin(x2),2))))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][2][0]= (-8*a*x1*Cot(x2))/Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2)
;Krzys[3][2][1]= 0
;Krzys[3][2][2]= 0
;Krzys[3][2][3]= ((3*Power(a,4) + 8*Power(a,2)*x1 + 8*Power(a,2)*Power(x1,2) + 8*Power(x1,4) + 4*Power(a,2)*(Power(a,2) + 2*(-1 + x1)*x1)*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][0]= 0
;Krzys[3][3][1]= (4*((-2 + x1)*Power(x1,4) + Power(a,4)*x1*Power(Cos(x2),4) - Power(a,2)*Power(x1,2)*Power(Sin(x2),2) + Power(a,2)*Power(Cos(x2),2)*(2*(-1 + x1)*Power(x1,2) + Power(a,2)*Power(Sin(x2),2))))/((Power(a,2) + (-2 + x1)*x1)*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][2]= ((3*Power(a,4) + 8*Power(a,2)*x1 + 8*Power(a,2)*Power(x1,2) + 8*Power(x1,4) + 4*Power(a,2)*(Power(a,2) + 2*(-1 + x1)*x1)*Cos(2*x2) + Power(a,4)*Cos(4*x2))*Cot(x2))/(2.*Power(Power(a,2) + 2*Power(x1,2) + Power(a,2)*Cos(2*x2),2))
;Krzys[3][3][3]= 0
;
}


if(coords==SPHCOORDS) {
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
}

if(coords==CYLCOORDS) {
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
}

if(coords==MINKCOORDS) {
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
}

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for cell ix,iy,iz
int 
fill_geometry(int ix,int iy,int iz,void *geom)
{
 struct geometry *ggg 
   = (struct geometry *) geom;

  pick_g(ix,iy,iz,ggg->gg);
  pick_G(ix,iy,iz,ggg->GG);
  pick_T(tmuup,ix,iy,iz,ggg->tup);
  pick_T(tmulo,ix,iy,iz,ggg->tlo);
  pick_T(emuup,ix,iy,iz,ggg->eup);
  pick_T(emulo,ix,iy,iz,ggg->elo);
  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;
  ggg->xxvec[0]=0.;
  ggg->xxvec[1]=get_x(ix,0);
  ggg->xxvec[2]=get_x(iy,1);
  ggg->xxvec[3]=get_x(iz,2);
  ggg->xx=ggg->xxvec[1];
  ggg->yy=ggg->xxvec[2];
  ggg->zz=ggg->xxvec[3];
  ggg->gdet=ggg->gg[3][4];

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for cell face ix,iy,iz in idim
int 
fill_geometry_face(int ix,int iy,int iz,int idim, void *geom)
{
  struct geometry *ggg 
    = (struct geometry *) geom;

  pick_gb(ix,iy,iz,idim,ggg->gg);
  pick_Gb(ix,iy,iz,idim,ggg->GG);
  if(idim==0)
    {
      /*
      //if these tetrads required at some point then they should be calculated on the go here
	
      pick_Tb(tmuupbx,ix,iy,iz,idim,ggg->tup);
      pick_Tb(tmulobx,ix,iy,iz,idim,ggg->tlo);
      pick_Tb(emuupbx,ix,iy,iz,idim,ggg->eup);
      pick_Tb(emulobx,ix,iy,iz,idim,ggg->elo);
      */
      ggg->xxvec[1]=get_xb(ix,0);
      ggg->xxvec[2]=get_x(iy,1);
      ggg->xxvec[3]=get_x(iz,2);
    }
  if(idim==1)
    {
      /*
      pick_Tb(tmuupby,ix,iy,iz,idim,ggg->tup);
      pick_Tb(tmuloby,ix,iy,iz,idim,ggg->tlo);
      pick_Tb(emuupby,ix,iy,iz,idim,ggg->eup);
      pick_Tb(emuloby,ix,iy,iz,idim,ggg->elo);
      */
      ggg->xxvec[1]=get_x(ix,0);
      ggg->xxvec[2]=get_xb(iy,1);
      ggg->xxvec[3]=get_x(iz,2);
    }
  if(idim==2)
    {
      /*
      pick_Tb(tmuupbz,ix,iy,iz,idim,ggg->tup);
      pick_Tb(tmulobz,ix,iy,iz,idim,ggg->tlo);
      pick_Tb(emuupbz,ix,iy,iz,idim,ggg->eup);
      pick_Tb(emulobz,ix,iy,iz,idim,ggg->elo);
      */
      ggg->xxvec[1]=get_x(ix,0);
      ggg->xxvec[2]=get_x(iy,1);
      ggg->xxvec[3]=get_xb(iz,2);
    }
  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;

  ggg->xxvec[0]=0.;
  ggg->xx=ggg->xxvec[1];
  ggg->yy=ggg->xxvec[2];
  ggg->zz=ggg->xxvec[3];

  ggg->gdet=ggg->gg[3][4];

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//fills geometry structure for cell ix,iy,iz in arbitrary metric
int 
fill_geometry_arb(int ix,int iy,int iz,void *geom,int COORDS)
{
  struct geometry *ggg 
    = (struct geometry *) geom;

  ldouble xxvec[4],xxvecBL[4];

  get_xx(ix,iy,iz,xxvec);
  coco_N(xxvec,xxvecBL,MYCOORDS,COORDS);

  calc_g_arb(xxvecBL,ggg->gg,COORDS);
  calc_G_arb(xxvecBL,ggg->GG,COORDS);

  calc_tetrades(ggg->gg,ggg->tup,ggg->tlo,COORDS);
  calc_ZAMOes(ggg->gg,ggg->eup,ggg->elo,COORDS);

  ggg->alpha=sqrt(-1./ggg->GG[0][0]);
  ggg->ix=ix;  ggg->iy=iy;  ggg->iz=iz;

  ggg->xxvec[0]=0.;
  ggg->xxvec[1]=xxvecBL[1];
  ggg->xxvec[2]=xxvecBL[2];
  ggg->xxvec[3]=xxvecBL[3];  

  ggg->xx=xxvecBL[1];
  ggg->yy=xxvecBL[2];
  ggg->zz=xxvecBL[3];
  
  ggg->gdet=ggg->gg[3][4];

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates orthonormal tetrad
//so far limited to grt.neq.0
int
calc_tetrades(ldouble g[][5], ldouble tmuup[][4], ldouble tmulo[][4],int coords)
{
  ldouble blob1,blob1sq;
  
  int method=1;
  int verbose=0;

  //numerical search for eigenvectors
  if(method==1)
    {

      double metric[]={g[0][0],g[0][1],g[0][2],g[0][3],
		       g[1][0],g[1][1],g[1][2],g[1][3],
		       g[2][0],g[2][1],g[2][2],g[2][3],
		       g[3][0],g[3][1],g[3][2],g[3][3]};		       

      gsl_matrix_view m = gsl_matrix_view_array (metric, 4, 4);     
      gsl_vector *eval = gsl_vector_alloc (4);
      gsl_matrix *evec = gsl_matrix_alloc (4, 4);     
      gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (4);       
      gsl_eigen_symmv (&m.matrix, eval, evec, w);     
      gsl_eigen_symmv_free (w);
       
      int i,j;
     
      for (i = 0; i < 4; i++)
	{
	  double eval_i 
	    = gsl_vector_get (eval, i);
	  gsl_vector_view evec_i 
	    = gsl_matrix_column (evec, i);
     
	  if(verbose) 
	    {
	      printf ("eigenvalue = %g\n", eval_i);
	      printf ("eigenvector = \n");
	      gsl_vector_fprintf (stdout, &evec_i.vector, "%g");
	    }

	  for(j=0;j<4;j++) tmulo[j][i]=gsl_matrix_get(evec,i,j)/sqrt(fabs(eval_i));

	}          


      //sorting to get maximal values on the diagonal
      int neworder[4];
      ldouble max;
      for(i=0;i<4;i++)
	{
	  max=-1.;
	  for(j=0;j<4;j++)
	    if(fabs(tmulo[i][j])>max)
	      {
		neworder[i]=j;
		max=fabs(tmulo[i][j]);
	      }
	}
      
      if(verbose)
	{
	  print_tensor(tmulo);
	  printf("no: %d %d %d %d\n",neworder[0],neworder[1],neworder[2],neworder[3]);
	}
      
      ldouble temp[4][4];
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  temp[neworder[i]][j]=tmulo[i][j];
      for(i=0;i<4;i++)
	for(j=0;j<4;j++)
	  tmulo[i][j]=temp[i][j];
      //changing orientation to positive
      for(i=0;i<4;i++) 
	if(tmulo[i][i]<0.) 
	  for(j=0;j<4;j++) tmulo[i][j]*=-1.;

      if(verbose) print_tensor(tmulo);

      //ortonormal -> lab
      for(i=0;i<4;i++)
	{     
	  for(j=0;j<4;j++)
	    tmuup[i][j]=tmulo[i][j];
	}

      //make them covariant
      for(i=0;i<4;i++)
	{	  
	  indices_21(tmuup[i],tmuup[i],g);
	  if(tmuup[i][i]<0.) 
	  for(j=0;j<4;j++) tmuup[i][j]*=-1.;
	}

      if(verbose) print_tensor(tmuup);

      gsl_vector_free (eval);
      gsl_matrix_free (evec);
    }

  //algorithm from HARM's tetrad.c from Mathematica assuming only grt non-zero
  if(method==2)
    {
      ldouble gtt,grr,grt,ghh,gpp;
      int i,j;

      gtt=g[0][0];
      grr=g[1][1];
      grt=g[0][1];
      ghh=g[2][2]; //gthth
      gpp=g[3][3]; //gphph

      //lab -> ortonormal
      blob1=grr - sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt;
      blob1sq=blob1*blob1;
      tmulo[0][0]= (sqrt(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt)/
		    pow((4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt)))*blob1sq,0.25));
      tmulo[0][1] =  (2.0*grt)/(sqrt(4.0*((grt)*(grt)) + (grr - gtt)*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt))*
				sqrt(fabs(grr - sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt)));
      tmulo[0][2] = 0.0;
      tmulo[0][3] = 0.0;
      tmulo[1][0] = (2.0*grt)/sqrt((4.0*((grt)*(grt)) + (grr - gtt)*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt))*
				   (grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt));
      tmulo[1][1] = sqrt((grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) - gtt)/
			 (sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt)))*(grr + sqrt(4.0*((grt)*(grt)) + ((grr - gtt)*(grr - gtt))) + gtt)));
      tmulo[1][2] = 0.0;
      tmulo[1][3] = 0.0;
      tmulo[2][0] = 0.0;
      tmulo[2][1] = 0.0;
      tmulo[2][2] = 1.0/fabs(sqrt(ghh));
      tmulo[2][3] = 0.0;
      tmulo[3][0] = 0.0;
      tmulo[3][1] = 0.0;
      tmulo[3][2] = 0.0;
      tmulo[3][3] = 1.0/fabs(sqrt(gpp));

      //ortonormal -> lab
      for(i=0;i<4;i++)
	{     
	  for(j=0;j<4;j++)
	    tmuup[i][j]=tmulo[i][j];
	}

      //make them covariant
      for(i=0;i<4;i++)
	{	  
	  indices_21(tmuup[i],tmuup[i],g);
	  if(tmuup[i][i]<0.) 
	  for(j=0;j<4;j++) tmuup[i][j]*=-1.;
	}
    }

  return 0;
}




//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates base vectors and 1-forms of LNRF to transform lab <--> LNRF
int
calc_ZAMOes(ldouble g[][5], ldouble emuup[][4], ldouble emulo[][4], int coords)
{
  ldouble e2nu,e2psi,e2mu1,e2mu2,omega;
  ldouble gtt,gtph,gphph,grr,gthth;
  int i,j;

  //TODO: to it in a general way
  //the following do not work but are not supposed to be used
  //so the error message suppressed
  //if(coords==KSCOORDS || coords==MKSCOORDS)
  //my_err("ZAMO tetrad calculation not implemented yet.\n");

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

  emuup[0][0]=sqrt(fabs(e2nu));
  emuup[1][1]=sqrt(fabs(e2mu1));
  emuup[2][2]=sqrt(fabs(e2mu2));
  emuup[0][3]=-omega*sqrt(fabs(e2psi));
  emuup[3][3]=sqrt(fabs(e2psi));

  emulo[3][0]=omega*1./sqrt(fabs(e2nu));
  emulo[0][0]=1./sqrt(fabs(e2nu));
  emulo[1][1]=1./sqrt(fabs(e2mu1));
  emulo[2][2]=1./sqrt(fabs(e2mu2));
  emulo[3][3]=1./sqrt(fabs(e2psi));

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for BL -> KS
int
coco_BL2KS(ldouble *xBL, ldouble *xKS)
{
  ldouble r=xBL[1];
  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;
  ldouble sqrta=sqrt(1.-a*a);
  //t
  xKS[0]=xBL[0]+2./sqrta*atanh(sqrta/(1.-r))+log(delta);
  //r
  xKS[1]=xBL[1];
  //theta
  xKS[2]=xBL[2];
  //phi
  xKS[3]=xBL[3]+a/sqrta*atanh(sqrta/(1.-r));

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> BL
int
coco_KS2BL(ldouble *xKS, ldouble *xBL)
{
  ldouble r=xKS[1];
  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;
  ldouble sqrta=sqrt(1.-a*a);
  //t
  xBL[0]=xKS[0]-2./sqrta*atanh(sqrta/(1.-r))-log(delta);
  //r
  xBL[1]=xKS[1];
  //theta
  xBL[2]=xKS[2];
  //phi
  xBL[3]=xKS[3]-a/sqrta*atanh(sqrta/(1.-r));

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for KS -> MKS1
int
coco_KS2MKS1(ldouble *xKS, ldouble *xMKS1)
{
  ldouble KSx0=xKS[0];
  ldouble KSx1=xKS[1];
  ldouble KSx2=xKS[2];
  ldouble KSx3=xKS[3];
  ldouble R0;

#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif

  xMKS1[0]
    = KSx0
    ;
  xMKS1[1]
    = log(KSx1-R0)
    ;
  xMKS1[2]
    = KSx2
    ;
  xMKS1[3]
    = KSx3
    ;

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MKS1 -> KS
int
coco_MKS12KS(ldouble *xMKS1, ldouble *xKS)
{
  ldouble x0=xMKS1[0];
  ldouble x1=xMKS1[1];
  ldouble x2=xMKS1[2];
  ldouble x3=xMKS1[3];
  ldouble R0;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif

  xKS[0]
    = x0
    ;
  xKS[1]
    = exp(x1) + R0
    ;
  xKS[2]
    = x2
    ;
  xKS[3]
    = x3
    ;

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MCYL1 -> CYL
int
coco_MCYL12CYL(ldouble *xMCYL1, ldouble *xCYL)
{
  ldouble x0=xMCYL1[0];
  ldouble x1=xMCYL1[1];
  ldouble x2=xMCYL1[2];
  ldouble x3=xMCYL1[3];
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif

  xCYL[0]
    = x0
    ;
  xCYL[1]
    = exp(x1) + R0
    ;
  xCYL[2]
    = x2
    ;
  xCYL[3]
    = x3
    ;

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for CYL -> MCYL1
int
coco_CYL2MCYL1(ldouble *xCYL, ldouble *xMCYL1)
{
  ldouble CYLx0=xCYL[0];
  ldouble CYLx1=xCYL[1];
  ldouble CYLx2=xCYL[2];
  ldouble CYLx3=xCYL[3];
  ldouble R0=0.;

#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif

  xMCYL1[0]
    = CYLx0
    ;
  xMCYL1[1]
    = log(CYLx1-R0)
    ;
  xMCYL1[2]
    = CYLx2
    ;
  xMCYL1[3]
    = CYLx3
    ;

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for SPH -> CARTESIAN
int
coco_SPH2MINK(ldouble *xSPH, ldouble *xMINK)
{
  ldouble r=xSPH[1];
  ldouble th=xSPH[2];
  ldouble ph=xSPH[3];
  
  xMINK[0]=xSPH[0];
  xMINK[1]=r*sin(th)*cos(ph);
  xMINK[2]=r*sin(th)*sin(ph);
  xMINK[3]=r*cos(th);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//converts coordinates
//for MINK -> SPH
int
coco_MINK2SPH(ldouble *xMINK, ldouble *xSPH)
{
  ldouble x=xMINK[1];
  ldouble y=xMINK[2];
  ldouble z=xMINK[3];
  
  xSPH[0]=xMINK[0];
  xSPH[1]=sqrt(x*x+y*y+z*z);
  xSPH[2]=acos(z/xSPH[1]);
  xSPH[3]=my_atan2(y,x);

  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//wrapper to convert coordinates
int
coco_N(ldouble *x1, ldouble *x2,int CO1, int CO2)
{
  if(CO1==CO2)
    {
      x2[0]=x1[0];
      x2[1]=x1[1];
      x2[2]=x1[2];
      x2[3]=x1[3];
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==CYLCOORDS)
    {
      ldouble r,th,ph;
      r=x1[1];
      th=x1[2];
      ph=x1[3];
            
      x2[0]=x1[0];
      x2[3]=ph;
      x2[2]=r*cos(th);
      x2[1]=r*sin(th);
    }
  else if((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==MCYL1COORDS)
    {
      //to CYL
      coco_MCYL12CYL(x1,x1); 

      //to BL/SPH
      ldouble R,z,ph;
      R=x1[1];
      z=x1[2];
      ph=x1[3];
            
      x2[0]=x1[0];
      x2[3]=ph;
      x2[1]=sqrt(R*R+z*z);
      x2[2]=asin(R/x2[1]);
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MCYL1COORDS)
    {
      ldouble r,th,ph;

      r=x1[1];
      th=x1[2];
      ph=x1[3];

      //to CYL
      x2[0]=x1[0];
      x2[3]=ph;
      x2[2]=r*cos(th);
      x2[1]=r*sin(th);
      
      //to MCYL1
      coco_CYL2MCYL1(x2,x2);  
    }
  else if((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==CYLCOORDS)
    {
      ldouble R,z,ph;
      R=x1[1];
      z=x1[2];
      ph=x1[3];
            
      x2[0]=x1[0];
      x2[3]=ph;
      x2[1]=sqrt(R*R+z*z);
      x2[2]=asin(R/x2[1]);
    }
  
  else if(((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==SPHCOORDS) ||
	  ((CO2==SCHWCOORDS || CO2==KERRCOORDS) && CO1==SPHCOORDS))
    {
      x2[0]=x1[0];
      x2[1]=x1[1];
      x2[2]=x1[2];
      x2[3]=x1[3];
    }
  else if((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==KSCOORDS)
    coco_BL2KS(x1,x2);
  else if (CO1==KSCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS))
    coco_KS2BL(x1,x2);
  else if (CO1==KSCOORDS && CO2==MKS1COORDS)
    coco_KS2MKS1(x1,x2);
  else if (CO1==MKS1COORDS && CO2==KSCOORDS)
    coco_MKS12KS(x1,x2);
  else if (CO1==MCYL1COORDS && CO2==CYLCOORDS)
    coco_MCYL12CYL(x1,x2);
  else if (CO1==CYLCOORDS && CO2==MCYL1COORDS)
    coco_CYL2MCYL1(x1,x2);  
  else if (CO1==MKS1COORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_MKS12KS(x1,x2);
      coco_KS2BL(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS) && CO2==MKS1COORDS)
    {
      coco_BL2KS(x1,x2);
      coco_KS2MKS1(x2,x2);
    }
  else if ((CO1==SCHWCOORDS || CO1==KERRCOORDS || CO1==SPHCOORDS) && CO2==MINKCOORDS)
    {
      coco_SPH2MINK(x1,x2);
    }
  else if (CO1==MINKCOORDS && (CO2==SCHWCOORDS || CO2==KERRCOORDS || CO2==SPHCOORDS))
    {
      coco_MINK2SPH(x1,x2);
    }
  else if (CO1==MINKCOORDS && CO2==MKS1COORDS)
    {
      coco_MINK2SPH(x1,x2);
      coco_BL2KS(x2,x2);
      coco_KS2MKS1(x2,x2);
    }
 else if (CO1==KSCOORDS & CO2==MINKCOORDS)
    {      
      coco_KS2BL(x1,x2);
      coco_SPH2MINK(x2,x2);
    }
  else if (CO1==MKS1COORDS & CO2==MINKCOORDS)
    {
      coco_MKS12KS(x1,x2);
      coco_KS2BL(x2,x2);
      coco_SPH2MINK(x2,x2);
    }
  else
    {
      printf("coco: %d -> %d\n",CO1,CO2);
      my_err("coco coordinate conversion not implemented\n");
    }
  return 0;
}


//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for BL -> KS
int
dxdx_BL2KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble t=xx[0];
  ldouble r=xx[1];
  ldouble th=xx[2];
  ldouble ph=xx[3];

  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[0][1]=2.*r/delta;
  dxdx[3][1]=a/delta;    

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> BL
int
dxdx_KS2BL(ldouble *xx, ldouble dxdx[][4])
{
  ldouble t=xx[0];
  ldouble r=xx[1];
  ldouble th=xx[2];
  ldouble ph=xx[3];

  ldouble a=BHSPIN;
  ldouble delta=r*r-2.*r+a*a;

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[0][1]=-2.*r/delta;
  dxdx[3][1]=-a/delta;    

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for KS -> MKS1
int
dxdx_KS2MKS1(ldouble *xx, ldouble dxdx[][4])
{
  ldouble KSx0=xx[0];
  ldouble KSx1=xx[1];
  ldouble KSx2=xx[2];
  ldouble KSx3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(KSx1-R0);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MKS1 -> KS
int
dxdx_MKS12KS(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MKS1COORDS)
  R0=MKS1R0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=exp(x1);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for CYL -> MCYL1
int
dxdx_CYL2MCYL1(ldouble *xx, ldouble dxdx[][4])
{
  ldouble CYLx0=xx[0];
  ldouble CYLx1=xx[1];
  ldouble CYLx2=xx[2];
  ldouble CYLx3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=1./(CYLx1-R0);

  return 0;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//calculates transformation matrices dxmu/dxnu
//for MCYL1 -> CYL
int
dxdx_MCYL12CYL(ldouble *xx, ldouble dxdx[][4])
{
  ldouble x0=xx[0];
  ldouble x1=xx[1];
  ldouble x2=xx[2];
  ldouble x3=xx[3];
  ldouble R0=0.;
#if(MYCOORDS==MCYL1COORDS)
  R0=MKS1R0;
#endif

  int i,j;
  for(i=0;i<4;i++)
    for(j=0;j<4;j++)
      dxdx[i][j]=delta(i,j);
  
  dxdx[1][1]=exp(x1);

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
	      printf("%10f ",g[i][j][k]);
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
	  printf("%10f ",g[j][k]);
	}
      printf("\n");
    }
  printf("\n");
  return -1;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//returns location of the horizone in BL
ldouble
r_horizon_BL(ldouble a)
{
  return 1.+sqrt(1-a*a);
}

//returns location of the co-rotating marginally bound orbit in BL
ldouble
r_mbound_BL(ldouble a)
{
  return 2.*(1.-a/2.+sqrt(1.-a));
}

//returns location of the photon orbit in BL
ldouble
r_photon_BL(ldouble a)
{
  return 2.*(1.-cosl(2./3.*acosl(-a)));
}

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************
//precalculates metric etc. and saves it to arrays
int
calc_metric()
{
  int ix,iy,iz;

  printf("Precalculating metrics... ");
  
  #pragma omp parallel for private(iy,iz) schedule (dynamic)
  for(ix=-NG;ix<NX+NG;ix++)
    {
      for(iy=-NG;iy<NY+NG;iy++)
	{
	  for(iz=-NG;iz<NZ+NG;iz++)
	    {
	      ldouble gloc[4][5];
	      ldouble Kr[4][4][4];
	      ldouble eup[4][4],elo[4][4];
	      ldouble tup[4][4],tlo[4][4];
	      ldouble xx[4];
	      int i,j,k;
	      //cell centers
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		    set_g(g,i,j,ix,iy,iz,gloc[i][j]);
	      for(j=0;j<3;j++)
		set_g(g,j,4,ix,iy,iz,calc_dlgdet(xx,j));
	      set_g(g,3,4,ix,iy,iz,calc_gdet(xx));

	      calc_ZAMOes(gloc,eup,elo,MYCOORDS);
	      calc_tetrades(gloc,tup,tlo,MYCOORDS);
	      
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    set_T(emuup,i,j,ix,iy,iz,eup[i][j]);
		    set_T(emulo,i,j,ix,iy,iz,elo[i][j]);
		    set_T(tmuup,i,j,ix,iy,iz,tup[i][j]);
		    set_T(tmulo,i,j,ix,iy,iz,tlo[i][j]);
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		    set_g(G,i,j,ix,iy,iz,gloc[i][j]);
	      
	      calc_Krzysie_at_center(ix,iy,iz,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKr(i,j,k,ix,iy,iz,Kr[i][j][k]);
	      	      
	      //x-faces
	      if(ix==-NG)
		{
		  xx[0]=0.;
		  xx[1]=get_xb(ix,0);
		  xx[2]=get_x(iy,1);
		  xx[3]=get_x(iz,2);
		  calc_g(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(gbx,i,j,ix,iy,iz,gloc[i][j],0);

		  calc_ZAMOes(gloc,eup,elo,MYCOORDS);
		  calc_tetrades(gloc,tup,tlo,MYCOORDS);

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      {
			/*
			set_Tb(emuupbx,i,j,ix,iy,iz,eup[i][j],0);
			set_Tb(emulobx,i,j,ix,iy,iz,elo[i][j],0);
			set_Tb(tmuupbx,i,j,ix,iy,iz,tup[i][j],0);
			set_Tb(tmulobx,i,j,ix,iy,iz,tlo[i][j],0);
			*/
		      }	      


		  calc_G(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(Gbx,i,j,ix,iy,iz,gloc[i][j],0);
		  for(j=0;j<3;j++)
		    set_gb(gbx,j,4,ix,iy,iz,calc_dlgdet(xx,j),0);
		  set_gb(gbx,3,4,ix,iy,iz,calc_gdet(xx),0);

		  /* Christoffels used only at centers
		  calc_Krzysie(xx,Kr);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      for(k=0;k<4;k++)
			set_gKrb(i,j,k,ix,iy,iz,Kr[i][j][k],0);
		  */


		}
	      xx[0]=0.;
	      xx[1]=get_xb(ix+1,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(gbx,i,j,ix+1,iy,iz,gloc[i][j],0);

	      calc_ZAMOes(gloc,eup,elo,MYCOORDS);
	      calc_tetrades(gloc,tup,tlo,MYCOORDS);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    /*
		    set_Tb(emuupbx,i,j,ix+1,iy,iz,eup[i][j],0);
		    set_Tb(emulobx,i,j,ix+1,iy,iz,elo[i][j],0);
		    set_Tb(tmuupbx,i,j,ix+1,iy,iz,tup[i][j],0);
		    set_Tb(tmulobx,i,j,ix+1,iy,iz,tlo[i][j],0);
		    */
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(Gbx,i,j,ix+1,iy,iz,gloc[i][j],0);
	      for(j=0;j<3;j++)
		set_gb(gbx,j,4,ix+1,iy,iz,calc_dlgdet(xx,j),0);
	      set_gb(gbx,3,4,ix+1,iy,iz,calc_gdet(xx),0);

	      /*
	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKrb(i,j,k,ix+1,iy,iz,Kr[i][j][k],0);
	      */

		  
	      //y-faces
	      if(iy==-NG)
		{
		  xx[0]=0.;
		  xx[1]=get_x(ix,0);
		  xx[2]=get_xb(iy,1);
		  xx[3]=get_x(iz,2);
		  calc_g(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(gby,i,j,ix,iy,iz,gloc[i][j],1);

		  calc_ZAMOes(gloc,eup,elo,MYCOORDS);
		  calc_tetrades(gloc,tup,tlo,MYCOORDS);

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      {
			/*
			set_Tb(emuupby,i,j,ix,iy,iz,eup[i][j],1);
			set_Tb(emuloby,i,j,ix,iy,iz,elo[i][j],1);
			set_Tb(tmuupby,i,j,ix,iy,iz,tup[i][j],1);
			set_Tb(tmuloby,i,j,ix,iy,iz,tlo[i][j],1);
			*/
		      }	      

		  calc_G(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(Gby,i,j,ix,iy,iz,gloc[i][j],1);
		  for(j=0;j<3;j++)
		    set_gb(gby,j,4,ix,iy,iz,calc_dlgdet(xx,j),1);
		  set_gb(gby,3,4,ix,iy,iz,calc_gdet(xx),1);

		  /*
		  calc_Krzysie(xx,Kr);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      for(k=0;k<4;k++)
			set_gKrb(i,j,k,ix,iy,iz,Kr[i][j][k],1);
		  */

		}
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_xb(iy+1,1);
	      xx[3]=get_x(iz,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(gby,i,j,ix,iy+1,iz,gloc[i][j],1);

	      calc_ZAMOes(gloc,eup,elo,MYCOORDS);
	      calc_tetrades(gloc,tup,tlo,MYCOORDS);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    /*
		    set_Tb(emuupby,i,j,ix,iy+1,iz,eup[i][j],1);
		    set_Tb(emuloby,i,j,ix,iy+1,iz,elo[i][j],1);
		    set_Tb(tmuupby,i,j,ix,iy+1,iz,tup[i][j],1);
		    set_Tb(tmuloby,i,j,ix,iy+1,iz,tlo[i][j],1);
		    */
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(Gby,i,j,ix,iy+1,iz,gloc[i][j],1);
	      for(j=0;j<3;j++)
		set_gb(gby,j,4,ix,iy+1,iz,calc_dlgdet(xx,j),1);
	      set_gb(gby,3,4,ix,iy+1,iz,calc_gdet(xx),1);
		  
	      /*
	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKrb(i,j,k,ix,iy+1,iz,Kr[i][j][k],1);
	      */

	      //z-faces
	      if(iz==-NG)
		{
		  xx[0]=0.;
		  xx[1]=get_x(ix,0);
		  xx[2]=get_x(iy,1);
		  xx[3]=get_xb(iz,2);
		  calc_g(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(gbz,i,j,ix,iy,iz,gloc[i][j],2);

		  calc_ZAMOes(gloc,eup,elo,MYCOORDS);
		  calc_tetrades(gloc,tup,tlo,MYCOORDS);

		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      {
			/*
			set_Tb(emuupbz,i,j,ix,iy,iz,eup[i][j],2);
			set_Tb(emulobz,i,j,ix,iy,iz,elo[i][j],2);
			set_Tb(tmuupbz,i,j,ix,iy,iz,tup[i][j],2);
			set_Tb(tmulobz,i,j,ix,iy,iz,tlo[i][j],2);
			*/
		      }	      

		  calc_G(xx,gloc);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      set_gb(Gbz,i,j,ix,iy,iz,gloc[i][j],2);
		  for(j=0;j<3;j++)
		    set_gb(gbz,j,4,ix,iy,iz,calc_dlgdet(xx,j),2);
		  set_gb(gbz,3,4,ix,iy,iz,calc_gdet(xx),2);

		  /*
		  calc_Krzysie(xx,Kr);
		  for(i=0;i<4;i++)
		    for(j=0;j<4;j++)
		      for(k=0;k<4;k++)
			set_gKrb(i,j,k,ix,iy,iz,Kr[i][j][k],2);
		  */

		}
	      xx[0]=0.;
	      xx[1]=get_x(ix,0);
	      xx[2]=get_x(iy,1);
	      xx[3]=get_xb(iz+1,2);
	      calc_g(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(gbz,i,j,ix,iy,iz+1,gloc[i][j],2);	  

	      calc_ZAMOes(gloc,eup,elo,MYCOORDS);
	      calc_tetrades(gloc,tup,tlo,MYCOORDS);

	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  {
		    /*
		    set_Tb(emuupbz,i,j,ix,iy,iz+1,eup[i][j],2);
		    set_Tb(emulobz,i,j,ix,iy,iz+1,elo[i][j],2);
		    set_Tb(tmuupbz,i,j,ix,iy,iz+1,tup[i][j],2);
		    set_Tb(tmulobz,i,j,ix,iy,iz+1,tlo[i][j],2);
		    */
		  }	      

	      calc_G(xx,gloc);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  set_gb(Gbz,i,j,ix,iy,iz+1,gloc[i][j],2);	  
	      for(j=0;j<3;j++)
		set_gb(gbz,j,4,ix,iy,iz+1,calc_dlgdet(xx,j),2);
	      set_gb(gbz,3,4,ix,iy,iz+1,calc_gdet(xx),2);

	      /*
	      calc_Krzysie(xx,Kr);
	      for(i=0;i<4;i++)
		for(j=0;j<4;j++)
		  for(k=0;k<4;k++)
		    set_gKrb(i,j,k,ix,iy,iz+1,Kr[i][j][k],2);
	      */

	      
	    }
	}
    }

  printf("done!\n");

  return 0;
}


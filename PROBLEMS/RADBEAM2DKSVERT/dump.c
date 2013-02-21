/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

v1=Tgas;
v2=Trad;
ldouble tau[3];
calc_tautot(pp,xxx,dx,tau);
v3=tau[0];
calc_tauabs(pp,xxx,dx,tau);
v4=tau[0];

ldouble nx,ny,nz,nlen,f;

  nx=Fx/E;
  ny=Fy/E;
  nz=Fz/E;

  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
 
   f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
v4=f;

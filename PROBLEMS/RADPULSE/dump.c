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
v4=get_u_scalar(aradx,ix,iy,iz);
v3=get_cflag(RADSOURCETYPEFLAG,ix,iy,iz);
  

//v4 filled by analytical_solution ()

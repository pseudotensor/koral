/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

ldouble cs  = sqrt(GAMMA*pp[UU]*GAMMAM1/pp[RHO]);

v5 = cs;

v6 = calc_lum(xxvecout[1],1)/calc_lumEdd()*(rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)*velGU2CGS(1.)*velGU2CGS(1.));


v1=Tgas;
v2=Trad;

ldouble tau[3],dxx[3]={1.,1.,1.};
calc_tautot(pp,xxx,dxx,tau);
v3=tau[0];
calc_tauabs(pp,xxx,dxx,tau);
v4=tau[0];


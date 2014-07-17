/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

ldouble cs  = sqrt(GAMMA*pp[UU]*GAMMAM1/pp[RHO]);

v5 = cs;
ldouble radlum,totlum;
calc_lum(xxvecout[1],3,&radlum,&totlum);
v6 = 1.e-40+radlum/calc_lumEdd()*(rhoGU2CGS(1.)*velGU2CGS(1.)*lenGU2CGS(1.)*lenGU2CGS(1.)*velGU2CGS(1.)*velGU2CGS(1.));

ldouble Rij[4][4];
calc_Rij(pp,&geom,Rij); 
v7=Rij[0][0]; //lab-frame energy density
//indices_2221(Rij,Rij,geom.gg);
//v7=-Rij[1][0]*sqrt(geom.gg[1][1]);//ortonormal flux

v1=Tgas;
v2=Trad;

ldouble tau[3],dxx[3]={1.,1.,1.};
calc_tautot(pp,xxx,dxx,tau);
v3=tau[0];
calc_tauabs(pp,xxx,dxx,tau);
v4=tau[0];


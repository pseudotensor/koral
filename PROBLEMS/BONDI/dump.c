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

//ldouble Rij[4][4];
//calc_Rij(pp,&geom,Rij);   
//v7=Rij[0][0]; //lab-frame energy density //(30)
//indices_2221(Rij,Rij,geom.gg);
//v7=-Rij[1][0]*sqrt(geom.gg[1][1]);//ortonormal flux

ldouble ehat,uconf[4];
calc_ff_Rtt(pp,&ehat,uconf,&geom);
ehat*=-1.;
v7=ehat;     
v8=calc_LTE_TfromE(ehat);

v1=Tgas;
v2=Trad;

ldouble tau[3],dxx[3]={1.,1.,1.};
calc_tautot(pp,&geomout,dxx,tau);
v3=tau[0];

calc_tauabs(pp,&geomout,dxx,tau);
v4=tau[0];



v9=pp[NF];//calc_ncompt_nphlab(pp,&geomout);


//temp
v8=get_u_scalar(cell_dt,ix,iy,iz);
v9=1./get_u_scalar(cell_tstepden,ix,iy,iz);

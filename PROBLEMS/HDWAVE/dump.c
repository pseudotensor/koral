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


ldouble om=1./CC*2.*Pi;
//hydro sound wave
#if (NWAVE==2)
ldouble myrho=RHOZERO*(1.+AAA*cos(KK*xx-om*t));
ldouble myuint=UINT*(1.+GAMMA*AAA*cos(KK*xx-om*t));
ldouble mycs=1./CC;
ldouble myvx=AAA*cos(KK*xx-om*t)*mycs;

v2=myrho;
v3=myuint;
v4=myvx;


#endif

//analytical profile of damped wave
//v3=RHO*(1+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx)));
//v4=ERAD*(1+DERE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DEIM/DERE*sin(OMRE*t-KK*xx)));



/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

//analytical solution is missing - look in anasol.c in previous versions to be placed as independent analytical_solution routine or something
v1=0.001;
v2=0.001;
v3=0.001;

ldouble podpierd=-(GG[0][0]-2.*ELL*GG[0][3]+ELL*ELL*GG[3][3]);
ldouble ulowert=-1./sqrt(podpierd);
if(ulowert<-1. || podpierd<0.)
  {
    v4=0.  ;
  }
 else
   {
     v4=ulowert;
   }

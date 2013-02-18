/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

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
v1=Tgas;

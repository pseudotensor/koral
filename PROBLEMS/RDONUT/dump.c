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
//tau in radius
//physical size of the cell
ldouble tautot[3];
calc_tautot(pp,xxvec,dx,tautot);
v2=tautot[0];
//prad/ptotal
ldouble prad=calc_LTE_EfromT(Trad)/3.;
ldouble pgas=GAMMAM1*calc_PEQ_ufromTrho(Tgas,rho);
				  v3=prad/(prad+pgas);

				  

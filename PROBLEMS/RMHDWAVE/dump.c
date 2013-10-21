/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

//analytical profile of damped wave

v1=RHOZERO+DRRE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DRIM/DRRE*sin(OMRE*t-KK*xx));
v1-=rho;
v2=UZERO+DURE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx)-DUIM/DURE*sin(OMRE*t-KK*xx));
v3=0.+DV1RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx))-DV1IM*exp(-OMIM*t)*sin(OMRE*t-KK*xx) ; 
v4=0.+DV2RE*exp(-OMIM*t)*(cos(OMRE*t-KK*xx))-DV2IM*exp(-OMIM*t)*sin(OMRE*t-KK*xx) ; 
v5=B2ZERO+DB2RE*exp(-OMIM*t)*cos(OMRE*t-KK*xx)-DB2IM*exp(-OMIM*t)*sin(OMRE*t-KK*xx); 
#ifdef RADIATION
v6=EEZERO+DEERE*exp(-OMIM*t)*cos(OMRE*t-KK*xx)-DEEIM*exp(-OMIM*t)*sin(OMRE*t-KK*xx);
v7=0.+DF1RE*exp(-OMIM*t)*cos(OMRE*t-KK*xx)-DF1IM*exp(-OMIM*t)*sin(OMRE*t-KK*xx) ; 
#endif





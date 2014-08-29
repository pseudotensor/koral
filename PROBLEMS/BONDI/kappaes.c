ldouble Tgas=T;
ldouble Trad=Tgas;

#ifdef NCOMPTONIZATION //adjust the Rosseland mean for the color temperature of radiation
Trad = calc_ncompt_Thatrad_full(pp,ggg);
#endif

ldouble kes=kappaCGS2GU(0.4)* rho;

//Klein-Nishina - Rybicki & Lightman
ldouble x=Trad/5.93e9;
ldouble lognum = log(1.+2.*x);
ldouble knfactor = 3./4.*((1.+x)/x/x/x*(2.*x*(1.+x)/(1.+2.*x)-lognum)+1./2./x*lognum - (1.+3.*x)/(1.+2.*x)/(1.+2.*x));

//override
knfactor=1.;

return kes*knfactor;


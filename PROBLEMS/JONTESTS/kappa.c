//ldouble calc_kappa(ldouble rho, ldouble T,ldouble x,ldouble y,ldouble z)

//ldouble pre=GAMMAM1*calc_PEQ_ufromTrho(T,rho);
//ldouble temp=(pre/rho);

ldouble temp=T;

/*
printf("arad: %e\n",SIGMA_RAD/4.);
printf("%e %e | %e\n",rho,temp,3.46764e-17 * rho * rho / pow(temp,4.)*sqrt(temp)); 
printf("KB: %e %e %e\n",K_BOLTZ,M_PROTON,MU_GAS);
getchar();
*/
return 1.*3.46764e-17 * rho * rho / pow(temp,4.)*sqrt(temp);


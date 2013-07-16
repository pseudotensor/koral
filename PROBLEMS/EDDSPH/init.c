//ix,iy,iz

ldouble mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

//pp[0] - rho
//pp[1] - uint
//pp[2]-pp[4] - vel_gas
//pp[5] - entropy
//pp[6] - Erf
//pp[7]-pp[9] - vel_rf

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);
ldouble xxsph[4];
coco_N(geom.xxvec,xxsph,MYCOORDS,SPHCOORDS);
ldouble rsph=xxsph[1];
/***********************************************/
//hydro atmosphere 
//set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,2);

#ifdef OPTTHIN
ldouble ALPHAP=(1.-FLUXBETA)/2.;
ldouble PSTAR=calc_PEQ_ufromTrho(TEMPSTAR,RHOSTAR)*GAMMAM1;
ldouble KKK = PSTAR / pow(RHOSTAR,GAMMA);
ldouble C0 = (GAMMA/GAMMAM1*PSTAR/RHOSTAR+1.)*pow(1.-2./RSTAR,ALPHAP);
ldouble pre = pow(KKK,-1./GAMMAM1)*pow(GAMMAM1/GAMMA*(C0*pow(rsph/(rsph-2.),ALPHAP)-1.),GAMMA/GAMMAM1);
ldouble rho = pow(pre/KKK,1./GAMMA);
ldouble temp=calc_PEQ_Tfromurho(pre/GAMMAM1,rho);

pp[RHO]=rho;
pp[UU]=calc_PEQ_ufromTrho(temp,pp[RHO]);
pp[VX]=pp[VY]=pp[VZ]=0.;
#endif


#ifdef OPTTHICK
ldouble rho = RHOSTAR;
ldouble temp=TEMPSTAR;

pp[RHO]=rho;
pp[UU]=calc_PEQ_ufromTrho(temp,pp[RHO]);
pp[VX]=pp[VY]=pp[VZ]=0.;

#endif

//rad atmosphere
#ifdef RADIATION

#ifdef OPTTHIN
//Abramowicz solution
ldouble alpha = asin(RSTAR/rsph*sqrt((1.-2./rsph)/(1.-2./RSTAR)));
ldouble I=FLUXBETA / KAPPA_ES_COEFF * pow((1.-2./RSTAR)/(1.-2./rsph),2.);
 //ldouble Fr = FLUXBETA / KAPPA_ES_COEFF / rsph / (rsph-2.);
ldouble Rtr = I/ RSTAR / (RSTAR-2.)*pow(sin(alpha),2.);
ldouble Rtt = 2.*I/ RSTAR / (RSTAR-2.)*(1.-cos(alpha));
//ortho-normal
pp[7]=Rtr; ////R^tr
pp[6]=Rtt; //R^tt
pp[8]=pp[9]=0.; //R^tinne
//converts upper row (R^tmu) from orthormal to coord. basis
prad_on2lab(pp,pp,&geom);
#endif

#ifdef OPTTHICK
//test
ldouble rhocgs=rhoGU2CGS(rho);
ldouble Tcgs=tempGU2CGS(temp);
ldouble ZZsun=1.;
ldouble kappaffcgs=6.4e22*rho/Tcgs/Tcgs/Tcgs/sqrt(Tcgs);
ldouble kappabfcgs=4.8e-24/1.67262158e-24/1.67262158e-24*rho/Tcgs/Tcgs/Tcgs/sqrt(Tcgs)*ZZsun;
ldouble kappaff=kappaCGS2GU(kappaffcgs)+kappaCGS2GU(kappabfcgs);

//Abramowicz solution
ldouble alpha = asin(RSTAR/rsph*sqrt((1.-2./rsph)/(1.-2./RSTAR)));
ldouble I=FLUXBETA / (kappaff+KAPPA_ES_COEFF) * pow((1.-2./RSTAR)/(1.-2./rsph),2.);
ldouble Rtr = I/ RSTAR / (RSTAR-2.)*pow(sin(alpha),2.);
ldouble Rtt = 2.*I/ RSTAR / (RSTAR-2.)*(1.-cos(alpha));
				       //ldouble Rtt = calc_LTE_EfromT(temp);
//Rtt from LTEQ
//ortho-normal
pp[7]=Rtr; ////R^tr
pp[6]=Rtt; //R^tt
pp[8]=pp[9]=0.; //R^tinne
//converts upper row (R^tmu) from orthormal to coord. basis
prad_on2lab(pp,pp,&geom);


				       //if(ix==1) {printf("%d | %e %e %e %e %e\n",ix,Rtr,Rtt,KKK,C0,rho,temp,kappaff);getchar();}

#endif
 
#endif

/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
//hd floors
check_floors_hd(pp,VELPRIM,&geom);
//to conserved
p2u(pp,uu,&geom);

/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);

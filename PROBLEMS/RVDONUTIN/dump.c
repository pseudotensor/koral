//v1 - 24
//..
//v7 - 30

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
ldouble tautot[3];
calc_tautot(pp,xxvec,dx,tautot);
v2=tautot[0];

//prad/ptotal
ldouble prad=calc_LTE_EfromT(Trad)/3.;
				  ldouble pgas=GAMMAM1*calc_PEQ_ufromTrho(Tgas,rho);
				  v3=prad/(prad+pgas);
				  ldouble PARAM=1.;
				  ldouble fdamptau=exp(-PARAM/tautot[0]/tautot[0]);
				  ldouble Rij[4][4];
				  calc_Rij(pp,&geom,Rij);
				  // skip boosts what reasonable when gas v<<1
				  boost22_lab2ff(Rij,Rij,pp,gg,GG);
				  trans22_cc2on(Rij,Rij,tup);

				  ldouble prad2=1./3.*Rij[0][0];//*fdamptau;
				  v4=prad2/prad;
				  v4=fdamptau;


				  ldouble nx,ny,nz,nlen,f;

				  nx=Fx/E;
				  ny=Fy/E;
				  nz=Fz/E;

				  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
 
				  f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  
				  

				  //				  v4=get_cflag(RADSOURCETYPEFLAG,ix,iy,iz);
				  v5=get_cflag(HDFIXUPFLAG,ix,iy,iz);
				  v6=get_cflag(RADFIXUPFLAG,ix,iy,iz);
				  v1=get_cflag(ENTROPYFLAG,ix,iy,iz);

				  v7=v1;

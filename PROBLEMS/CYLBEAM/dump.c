/*

int
fprint_profiles(ldouble t, ldouble totmass)
{
 (...)

*/ 

v1=Tgas;
//tau in radius
//physical size of the cell
ldouble tautot[3];
calc_tautot(pp,xxvec,dx,tautot);
v1=tautot[0];
//prad/ptotal
ldouble prad=calc_LTE_EfromT(Trad)/3.;
ldouble pgas=GAMMAM1*calc_PEQ_ufromTrho(Tgas,rho);
				  v3=prad/(prad+pgas);

				  
ldouble nx,ny,nz,nlen,f;

  nx=Fx/E;
  ny=Fy/E;
  nz=Fz/E;

  nlen=sqrt(nx*nx+ny*ny+nz*nz);
  
 
   f=(3.+4.*(nx*nx+ny*ny+nz*nz))/(5.+2.*sqrt(4.-3.*(nx*nx+ny*ny+nz*nz)));  


				  v2=get_cflag(RADSOURCETYPEFLAG,ix,iy,iz);
				  v3=get_cflag(HDFIXUPFLAG,ix,iy,iz);
				   v4=get_cflag(RADFIXUPFLAG,ix,iy,iz);

v3=f;

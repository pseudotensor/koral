
ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

/***********************************************/
//background
pp[RHO]=RHOZERO;
pp[UU]=calc_PEQ_ufromTrho(TEMPZERO,RHOZERO);
pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;

ldouble w1=exp(-((geom.xx-XBLOB1)*(geom.xx-XBLOB1) + geom.yy*geom.yy)/SIZEBLOB1/SIZEBLOB1);
ldouble w2=exp(-((geom.xx-XBLOB2)*(geom.xx-XBLOB2) + geom.yy*geom.yy)/SIZEBLOB2/SIZEBLOB2);


pp[RHO]*=(1. + BLOBMAG1*w1 + BLOBMAG2*w2);

pp[UU]=(1.*calc_PEQ_ufromTrho(TEMPZERO,pp[RHO]) + 
	w1*calc_PEQ_ufromTrho(TEMPBLOB1,pp[RHO]) +
	w2*calc_PEQ_ufromTrho(TEMPBLOB2,pp[RHO]))/(1.+w1+w2);


if(w1>0.01 && w2>0.01) pp[VX]=.5*(VELXBLOB1+VELXBLOB2);
 else if(w1>0.01) pp[VX]=VELXBLOB1;
 else if(w2>0.01) pp[VX]=VELXBLOB2;

if(w1>0.01 && w2>0.01) pp[VY]=.5*(VELYBLOB1+VELYBLOB2);
 else if(w1>0.01) pp[VY]=VELYBLOB1;
 else if(w2>0.01) pp[VY]=VELYBLOB2;

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

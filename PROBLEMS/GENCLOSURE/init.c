ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomCART;
fill_geometry_arb(ix,iy,iz,&geomCART,MINKCOORDS);

/***********************************************/
//background
pp[RHO]=RHOZERO;
pp[UU]=calc_PEQ_ufromTrho(TEMPZERO,RHOZERO);
pp[VX]=0.;
pp[VY]=0.;
pp[VZ]=0.;


ldouble w1=exp(-((geomCART.xx-XBLOB1)*(geomCART.xx-XBLOB1) + (geomCART.zz-YBLOB1)*(geomCART.zz-YBLOB1) )/SIZEBLOB1/SIZEBLOB1);
ldouble w2=exp(-((geomCART.xx-XBLOB2)*(geomCART.xx-XBLOB2) + (geomCART.zz-YBLOB2)*(geomCART.zz-YBLOB2) )/SIZEBLOB2/SIZEBLOB2);

//printf("%d %d > %f %f %f > %f %f\n",ix,iy,geomCART.xx,geomCART.yy,geomCART.zz,XBLOB1,YBLOB1);getch();


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

#ifdef RADIATION
ldouble temp=calc_PEQ_Tfromurho(pp[UU],pp[RHO]);

pp[EE0]=calc_LTE_EfromT(temp);
pp[FX0]=pp[FY0]=pp[FZ0]=0.;

//test
//pp[EE0]=1.;
//pp[FX0]=0.1;

#ifdef NCOMPTONIZATION
pp[NF0]=calc_NFfromE(pp[EE0]);
#endif

#endif

#ifdef TRACER
pp[TRA]=1.;
#endif

/***********************************************/
//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
//hd floors
//check_floors_hd(pp,VELPRIM,&geom);
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

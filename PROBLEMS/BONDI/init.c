

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);


ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE,ur;  
ldouble xx,yy,zz;
ldouble uu[NV];

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);
ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

ldouble pp[NV],T;

/************************/

int iix=ix;
#ifdef FLAT
iix=NX;
#endif

rho=get_u(pproblem1,RHO,iix,iy,iz);
uint=get_u(pproblem1,UU,iix,iy,iz);
ur=get_u(pproblem1,VX,iix,iy,iz);

#ifdef FLAT
ur=0.;
#endif

#ifdef INFLOW
rho*=1.e-10;
uint*=1.e-10;
#endif

//four-vel in BL
ldouble ucon[4]={0.,ur,0.,0.};
//rad. four-vel in BL
ldouble urfcon[4]={0.,ur,0.,0.}; 

pp[0]=rho;
pp[1]=uint;
pp[2]=ucon[1];
pp[3]=ucon[2];
pp[4]=ucon[3];
pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
E=get_u(pproblem1,EE0,iix,iy,iz); 
pp[6]=E;
pp[7]=urfcon[1];
pp[8]=urfcon[2];
pp[9]=urfcon[3]; 

#ifdef NCOMPTONIZATION
pp[NF0]=calc_NFfromE(pp[EE0]);
#endif
#endif	 

p2u(pp,uu,&geom);	 

/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }


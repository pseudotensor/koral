

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);


ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
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

ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut;

//at outern boundary
r=RMAX;
ur=-sqrt(2./r);
rho0=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
Tgas0=TGAS0;
            
//at given cell
r=geomBL.xx;
ur=-sqrt(2./r);    
printf("%d %e %e\n",ix,r,ur);

rho=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
Tgas=Tgas0*pow(rho/rho0,GAMMA-1.);      

uint=calc_PEQ_ufromTrho(Tgas,rho);

pgas=K_BOLTZ*rho*Tgas/MU_GAS/M_PROTON;
prad=PRADGAS*pgas;
E=prad*3.;

//four-vel in BL
ldouble ucon[4]={0.,ur,0.,0.};
conv_vels(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);

//rad. four-vel in BL
ldouble urfcon[4]={0.,0.,0.,0.};
conv_vels(urfcon,urfcon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);


pp[0]=rho;
pp[1]=uint;
pp[2]=ucon[1];
pp[3]=ucon[2];
pp[4]=ucon[3];
pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
pp[6]=E;
pp[7]=urfcon[1];
pp[8]=urfcon[2];
pp[9]=urfcon[3]; 
#endif	    
p2u(pp,uu,&geom);	 


/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }


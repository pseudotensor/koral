

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

/*
 ldouble rhos,Tgass,ur,Tgas,Trad,r,prad,pgas,ut,vx,Be,Kappa,urs;

 //at RBONDI
 r=RMAX;
 Be=(5.-3.*GAMMA)/(4.*(GAMMA-1.))/RBONDI;
 urs=-sqrt(Be/(5.-3.*GAMMA)*2.*GAMMAM1);
 rhos=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(RBONDI)*lenGU2CGS(RBONDI)*velGU2CGS(urs)));
 Kappa=GAMMAM1/pow(rhos,GAMMAM1)*(Be+1./RBONDI-0.5*urs*urs);
 pgas=Kappa*pow(rhos,GAMMA);
 Tgass=pgas/K_BOLTZ/rhos*MU_GAS*M_PROTON;


 //at given cell
 r=geomBL.xx;

 ur=urs*pow(RBONDI/r,1.5);
 rho=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
 Tgas=Tgass*pow(rho/rhos,GAMMA-1.);      
 uint=calc_PEQ_ufromTrho(Tgas,rho);
 pgas=K_BOLTZ*rho*Tgas/MU_GAS/M_PROTON;
 prad=PRADGAS*pgas;
 E=prad*3.;
*/

/* old
//at outern boundary
r=RMAX;
ur=-sqrt(2./r);
rho0=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
Tgas0=TGAS0;
            
//at given cell
r=geomBL.xx;
ur=-sqrt(2./r);    

rho=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
Tgas=Tgas0*pow(rho/rho0,GAMMA-1.);      

uint=calc_PEQ_ufromTrho(Tgas,rho);

pgas=K_BOLTZ*rho*Tgas/MU_GAS/M_PROTON;
prad=PRADGAS*pgas;
E=prad*3.;
*/

rho=get_u(pproblem1,RHO,ix,iy,iz);
uint=get_u(pproblem1,UU,ix,iy,iz);
ur=get_u(pproblem1,VX,ix,iy,iz);

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
E=get_u(pproblem1,EE0,ix,iy,iz);
pp[6]=E;
pp[7]=urfcon[1];
pp[8]=urfcon[2];
pp[9]=urfcon[3]; 
#endif	 

//transforming primitives from BL to MYCOORDS
trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,&geomBL,&geom);
   
p2u(pp,uu,&geom);	 


/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }


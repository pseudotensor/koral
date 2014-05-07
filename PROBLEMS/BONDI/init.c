

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

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
r=MAXX;
ur=-sqrtl(2./r);
rho0=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
Tgas0=TGAS0;
            
//at given cell
r=xx;
ur=-sqrtl(2./r);    
ut=sqrtl((-1.-ur*ur*gg[1][1])/gg[0][0]);
vx=ur/ut;  
rho=rhoCGS2GU(-MDOT*MDOTEDD/(4.*Pi*lenGU2CGS(r)*lenGU2CGS(r)*velGU2CGS(ur)));
Tgas=Tgas0*powl(rho/rho0,GAMMA-1.);      

uint=calc_PEQ_ufromTrho(Tgas,rho);

pgas=K_BOLTZ*rho*Tgas/MU_GAS/M_PROTON;
prad=PRADGAS*pgas;
E=prad*3.;


Fz=Fy=Fx=0.;
pp[0]=rho;
pp[1]=uint;
pp[2]=vx;

/*
  pp[0]=RHOAMB;
  pp[1]=calc_PEQ_ufromTrho(TAMB,RHOAMB);
  pp[3]=0.;
*/

pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(rho,uint);
#ifdef RADIATION
pp[6]=E;
pp[7]=Fx;
pp[8]=Fy;
pp[9]=Fz; 
prad_ff2lab(pp,pp,&geom);  

#endif	    
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

//if(isnan(get_u(p,5,ix,iy,iz))) {printf("pr: %d %d %d S: %Le\n",ix,iy,iz,0.);getchar();}

//mark initialy succesfull u2p_hot step
set_cflag(0,ix,iy,iz,0);

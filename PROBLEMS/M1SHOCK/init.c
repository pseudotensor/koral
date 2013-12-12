//KORAL - init.c

  struct geometry geom;
  fill_geometry(ix,iy,iz,&geom);

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];
ldouble gg[4][5],GG[4][5],tup[4][4],tlo[4][4];
ldouble pp[NV],T;

xx=get_x(ix,0);
yy=get_x(iy,1);
zz=get_x(iz,2);

pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
calc_tetrades(gg,tup,tlo,MYCOORDS);

/************************/
/************************/
/************************/
/* modify below */
      
ldouble rho0,Tgas0,ur,Tgas,Trad,r,rcm,prad,pgas,vx,ut,ux;

Fx=Fz=Fy=0.;
pp[0]=1.;
pp[1]=calc_PEQ_ufromTrho(TEMPIN,pp[RHO]);
ldouble ucon[4]={0.,0.,0.,0.};
conv_vels(ucon,ucon,VEL3,VELPRIM,gg,GG);
pp[2]=ucon[1]; 
pp[3]=ucon[2];
pp[4]=ucon[3]; 

pp[5]=calc_Sfromu(pp[RHO],pp[UU]);

pp[EE0]=calc_LTE_EfromT(TEMPIN);
pp[FX0]=Fx;
pp[FY0]=Fy;
pp[FZ0]=Fz; 

prad_ff2lab(pp,pp,&geom);
p2u(pp,uu,&geom);	

/* modify above */
/***********************************************/
/***********************************************/
/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
//mark initialy succesfull u2p_hot step
set_cflag(0,ix,iy,iz,0);

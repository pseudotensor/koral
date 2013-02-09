
ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];
ldouble gg[4][5],eup[4][4],elo[4][4],GG[4][5];
ldouble pp[NV],T,xxvec[4],xxvecBL[4];

get_xx(ix,iy,iz,xxvec);

coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);

xx=xxvecBL[1];
yy=xxvecBL[2];
zz=xxvecBL[3];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);


//KERR metric
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);

/************************/

//ambient
set_hdatmosphere(pp,xxvec,gg,GG,0);



//BL free-fall velocity
    ldouble ucon[4];
    ldouble r=xx;
    ucon[1]=-sqrt(2./r)*(1.-2./r);
    ucon[2]=ucon[3]=0.;

    conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);

    trans2_coco(xxvecBL,ucon,ucon,BLCOORDS,MYCOORDS);

    conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);

    pp[2]=ucon[1];
    pp[3]=ucon[2];
    pp[4]=ucon[3];

    
//calculating entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);

//testing if primitives make sense
check_floors_hd(pp,VELPRIM,gg,GG);
//end of floor section	      
	      
p2u(pp,uu,gg,GG);	 
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

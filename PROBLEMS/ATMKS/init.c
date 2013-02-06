
ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV];
ldouble gg[4][5],eup[4][4],elo[4][4],GG[4][5];
ldouble pp[NV],T,xxvec[4];

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvec,KSCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

/************************/

//ambient
set_hdatmosphere(pp,xxvec,gg,GG,0);

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

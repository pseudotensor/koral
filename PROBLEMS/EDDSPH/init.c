//ix,iy,iz

ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

//pp[0] - rho
//pp[1] - uint
//pp[2]-pp[4] - vel_gas
//pp[5] - entropy
//pp[6] - Erf
//pp[7]-pp[9] - vel_rf

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

/***********************************************/
//hydro atmosphere
set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,2);
//rad atmosphere
#ifdef RADIATION
set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,1);
#endif

/***********************************************/
//to cartesian
ldouble xxsph[4];
coco_N(geom.xxvec,xxsph,MYCOORDS,SPHCOORDS);

//distance from the center
ldouble dist = sqrt((xxsph[1]-BLOBR)*(xxsph[1]-BLOBR));

//increase rho
ldouble factor=(1.+BLOBP*exp(-dist*dist/BLOBW/BLOBW));
pp[0]*=factor;

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

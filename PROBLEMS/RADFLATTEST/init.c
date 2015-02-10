

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

ldouble pp[NV],uu[NV];

/************************/

if(geom.xx<0.3)
  pp[0]=LEFTRHO;
else
  pp[0]=RHOAMB;

pp[1]=UUAMB;
pp[2]=0.;
pp[3]=0.;
pp[4]=0.;
pp[5]=calc_Sfromu(pp[0],pp[1]);
#ifdef RADIATION
pp[6]=1.e-5;
pp[7]=0.;
pp[8]=0.;
pp[9]=0.;
#endif	 

p2u(pp,uu,&geom);	 


/***********************************************/

int iv;

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }


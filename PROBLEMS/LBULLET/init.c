
//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

int iv;
ldouble pp[NV], uu[NV];

for(iv=0; iv<NV; iv++) 
  { 
    pp[iv]=get_u(pproblem1,iv,ix,iy,iz);
  }

//to conserved
p2u(pp,uu,&geom);



/***********************************************/

for(iv=0;iv<NV;iv++)
  {
    set_u(u,iv,ix,iy,iz,uu[iv]);
    set_u(p,iv,ix,iy,iz,pp[iv]);
  }

//entropy
update_entropy(ix,iy,iz,0);
set_cflag(0,ix,iy,iz,0);

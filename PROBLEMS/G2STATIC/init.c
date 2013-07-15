
ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

/***********************************************/
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);
struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

/***********************************************/
//hydro atmosphere
//set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,2);

//Sgr A* atmosphere
#ifdef DONUT
int anret=donut_analytical_solution(pp,geomBL.xxvec,geomBL.gg,geomBL.GG);
if(anret<0) //atmosphere
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);   
  }
 else
  {
    //transforming primitives from BL to MYCOORDS
    trans_phd_coco(pp, pp, KERRCOORDS, MYCOORDS,geomBL.xxvec,geomBL.gg,geomBL.GG,geom.gg,geom.GG);     
  }
#else
set_sgradisk(pp,geom.xxvec,&geom,&geomBL);
#endif

/***********************************************/
//imposing cloud in rho and velocities
//cloud set up in prepinit.c and saved to pproblem

ldouble ucon[4];
ucon[1]=get_u(pproblem,VX,ix,iy,iz);
ucon[2]=get_u(pproblem,VY,ix,iy,iz);
ucon[3]=get_u(pproblem,VZ,ix,iy,iz);
ldouble atmrho=pp[0];
ldouble clrho=get_u(pproblem,RHO,ix,iy,iz);
pp[0]=atmrho+clrho;

pp[2]=(pp[2]*atmrho + ucon[1]*clrho ) / (atmrho + clrho);
pp[3]=(pp[3]*atmrho + ucon[2]*clrho ) / (atmrho + clrho);
pp[4]=(pp[4]*atmrho + ucon[3]*clrho ) / (atmrho + clrho);  

#ifdef TRACER
//tracer : fraction of cloud gas //test
pp[TRA]=clrho/pp[0];

if(pp[TRA]>MINTRACE)
  {
    pp[2]=ucon[1];
    pp[3]=ucon[2];
    pp[4]=ucon[3];    
  }
#endif
/***********************************************/

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


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
//cloud set up in prepinit.c and saved to pproblem1

ldouble ucon[4];
ucon[1]=get_u(pproblem1,VX,ix,iy,iz);
ucon[2]=get_u(pproblem1,VY,ix,iy,iz);
ucon[3]=get_u(pproblem1,VZ,ix,iy,iz);
ldouble diskrho=pp[0];
ldouble clrho=get_u(pproblem1,RHO,ix,iy,iz);
pp[0]=diskrho+clrho;


#ifdef TRACER

//tracer : fraction of cloud gas //test
ldouble tracer=clrho/(diskrho+clrho);

diskrho*=1.-step_function(tracer-DISKDAMPPAR1,DISKDAMPPAR2);

pp[0]=diskrho+clrho;
tracer=clrho/(diskrho+clrho);
pp[TRA]=tracer;

//to decrease cloud temperature
//pp[UU]= pp[UU]*(1. - tracer/(1/0.9));

#endif

pp[2]=(pp[2]*diskrho + ucon[1]*clrho ) / (diskrho + clrho);
pp[3]=(pp[3]*diskrho + ucon[2]*clrho ) / (diskrho + clrho);
pp[4]=(pp[4]*diskrho + ucon[3]*clrho ) / (diskrho + clrho);  

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

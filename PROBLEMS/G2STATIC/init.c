
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
//cloud
ldouble clix,cliy,cliz;
/*
clix=NX*0.66;
cliy=NY*0.3;
cliz=NZ;
*/
clix=NX*0.8;
cliy=NY/2;
cliz=NZ*.25;

ldouble clxx[4];
get_xx(clix,cliy,cliz,clxx);

//to cartesian
ldouble xxmink[4],clxxmink[4];
coco_N(geom.xxvec,xxmink,MYCOORDS,MINKCOORDS);
coco_N(clxx,clxxmink,MYCOORDS,MINKCOORDS);

//distance from the center
ldouble dist = sqrt((xxmink[1]-clxxmink[1])*(xxmink[1]-clxxmink[1])+
		    (xxmink[2]-clxxmink[2])*(xxmink[2]-clxxmink[2])+
		    (xxmink[3]-clxxmink[3])*(xxmink[3]-clxxmink[3]));


//increase rho
ldouble mag=CLMAG;
ldouble factor=(1.+mag*exp(-dist*dist/CLWIDTH/CLWIDTH));
ldouble atmrho = pp[0];
ldouble clrho = (factor-1.)*atmrho;
pp[0] =atmrho+clrho;

#ifdef TRACER
//tracer : fraction of cloud gas 
pp[TRA]=clrho/pp[0];
#endif

//velocity
ldouble OmKep = 1./sqrt(geomBL.xx*geomBL.xx*geomBL.xx);

//ldouble ucon[4]={0.,0,+OmKep,0.};
ldouble ucon[4]={0.,0,0,-OmKep};

conv_vels(ucon,ucon,VEL3,VEL4,geomBL.gg,geomBL.GG);
trans2_coco(geomBL.xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
conv_vels(ucon,ucon,VEL4,VELPRIM,geom.gg,geom.GG);

pp[2]=(pp[2]*atmrho + ucon[1]*clrho ) / (atmrho + clrho);
pp[3]=(pp[3]*atmrho + ucon[2]*clrho ) / (atmrho + clrho);
pp[4]=(pp[4]*atmrho + ucon[3]*clrho ) / (atmrho + clrho);    


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

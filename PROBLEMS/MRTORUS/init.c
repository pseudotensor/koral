
ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4],xxvecBL[4];
ldouble pp[NV],ppback[NV],T;

//geometries
get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvecBL,MYCOORDS,BLCOORDS);
xx=xxvecBL[1];
yy=xxvecBL[2];
zz=xxvecBL[3];

ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4],tlo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);

struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvecBL,ggBL,KERRCOORDS);
calc_G_arb(xxvecBL,GGBL,KERRCOORDS);

ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);

//donut formulae
ldouble podpierd=-(GGBL[0][0]-2.*ELL*GGBL[0][3]+ELL*ELL*GGBL[3][3]);
ldouble ut=-1./sqrt(podpierd);

ut/=UTPOT; //rescales rin
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

if(ut<-1 || podpierd<0. || xx<4.) //outside donut
  {
    //ambient
    set_hdatmosphere(pp,xxvec,gg,GG,0);
  }
 else //inside
   {
    //ambient
    set_hdatmosphere(ppback,xxvec,gg,GG,0);

    ldouble h=-1./ut;
    ldouble eps=(h-1.)/GAMMA;
    rho=powl(eps*(GAMMA-1.)/KKK,1./(GAMMA-1.));
    uint=rho*eps;
    uphi=-ELL*ut;
    uT=GGBL[0][0]*ut+GGBL[0][3]*uphi;
    uPhi=GGBL[3][3]*uphi+GGBL[0][3]*ut;
    Vphi=uPhi/uT;
    Vr=0.;

    //4-velocity in BL transformed to MYCOORDS
    ldouble ucon[4]={0.,-Vr,0.,Vphi};
    conv_vels(ucon,ucon,VEL3,VELPRIM,ggBL,GGBL);
   
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

    //transforming primitives from BL to MYCOORDS
    trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvecBL,&geomBL,&geom);
    
#ifdef MAGNFIELD 
    //MYCOORDS vector potential to calculate B's
    ldouble Acov[4];
    Acov[0]=Acov[1]=Acov[2]=0.;
    Acov[3]=my_max(pp[RHO]/4.e-20-0.02,0.)*sqrt(1.e-23);

    pp[B1]=0.;
    pp[B2]=0.;
    pp[B3]=Acov[3];
#endif

  }

#ifdef RADIATION
    //setting up radiative energy density to residual value - 
    //should immetiadelly increase sucking out the gas internal energy in opt.thick regions
    pp[EE0]=pp[UU]/1.e6;
    pp[FX0]=pp[FY0]=pp[FZ0]=0.;    
#endif

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
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

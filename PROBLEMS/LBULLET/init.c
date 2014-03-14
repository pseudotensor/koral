int init_dsandvels_limotorus(FTYPE r, FTYPE th, FTYPE a, FTYPE *rhoout, FTYPE *uuout, FTYPE *ell);

ldouble rho,mx,my,mz,m,E,uint,pgas,Fx,Fy,Fz,pLTE,ell;  
ldouble uu[NV], pp[NV],ppback[NV],T,uintorg;
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;

//geometries
struct geometry geom;
fill_geometry(ix,iy,iz,&geom);

struct geometry geomBL;
fill_geometry_arb(ix,iy,iz,&geomBL,KERRCOORDS);

ldouble r=geomBL.xx;
ldouble th=geomBL.yy;

init_dsandvels_limotorus(r, th, BHSPIN, &rho, &uint, &ell);
uintorg=uint;

if(rho<0.) //outside donut
  {
    //ambient
    set_hdatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,geom.xxvec,geom.gg,geom.GG,0);
#endif
  }
 else //inside donut
   {
    //ambient
    set_hdatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,geom.xxvec,geom.gg,geom.GG,0);
#endif

    uint=LT_KAPPA * pow(rho, LT_GAMMA) / (LT_GAMMA - 1.);
    pgas = GAMMAM1 * uint;
    ell*=-1.;

    ldouble ult,ulph,ucov[4],ucon[4];
    ulph = sqrt(-1./(geomBL.GG[0][0]/ell/ell + 2./ell*geomBL.GG[0][3] + geomBL.GG[3][3]));
    ult = ulph / ell;

    ucov[0]=ult;
    ucov[1]=0.;
    ucov[2]=0.;
    ucov[3]=ulph;
    
    indices_12(ucov,ucon,geomBL.GG);

    conv_vels_ut(ucon,ucon,VEL4,VELPRIM,geomBL.gg,geomBL.GG);
   
    pp[0]=my_max(rho,ppback[0]); 
    pp[1]=my_max(uint,ppback[1]);
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];

#ifdef MAGNFIELD//setting them zero not to break the following coordinate transformation
    pp[B1]=pp[B2]=pp[B3]=0.; 
#endif

   }

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

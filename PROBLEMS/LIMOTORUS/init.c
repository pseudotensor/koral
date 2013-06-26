
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

#ifdef BLOB
    pp[0]*=(1.+10.*exp(-(xx-14)*(xx-14)/.1));
#endif

#ifdef RADIATION
    set_radatmosphere(pp,xxvec,gg,GG,0);
#endif
  }
 else //inside
   {
    //ambient
    set_hdatmosphere(ppback,xxvec,gg,GG,0);
#ifdef RADIATION
    set_radatmosphere(ppback,xxvec,gg,GG,0);
#endif

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
   
    pp[2]=ucon[1]; 
    pp[3]=ucon[2];
    pp[4]=ucon[3];
    pp[0]=my_max(rho,ppback[0]); 

    //    pp[0]=ppback[0];

    pp[1]=my_max(uint,ppback[1]);


    ldouble P,aaa,bbb;
    P=GAMMAM1*uint;
    //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
    aaa=4.*SIGMA_RAD/3.;
    bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
    ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
    ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;

    E=calc_LTE_EfromT(T4);
    Fx=Fy=Fz=0.;
    uint=calc_PEQ_ufromTrho(T4,rho);

#ifdef HDDONUTASWITHRAD
    pp[1]=my_max(uint,ppback[1]);
#endif

#ifdef RADIATION
    pp[1]=my_max(uint,ppback[1]);
    pp[6]=my_max(E,ppback[6]);

    pp[7]=Fx;
    pp[8]=Fy;
    pp[9]=Fz;


 
    //now estimating flux in r,theta plane: F = -1/chi E,i in lab coordinates
    ldouble rho=pp[RHO];
    ldouble u=pp[1];  
    ldouble E=pp[6];  
    ldouble pr=(GAMMA-1.)*(u);
    ldouble T=pr*MU_GAS*M_PROTON/K_BOLTZ/rho;
    ldouble xx=get_x(ix,0);
    ldouble yy=get_x(iy,1);
    ldouble zz=get_x(iz,2);

    ldouble kappa=calc_kappa(rho,T,xx,yy,zz);
    ldouble chi=kappa+calc_kappaes(rho,T,xx,yy,zz);  
    
    ldouble xxvectemp[4]={xxvec[0],xxvec[1],xxvec[2],xxvec[3]};
    ldouble pptemp[NV],E1,E2,ggt[4][5],GGt[4][5];
    int anret,anretmin=0;

    //r dimension
    xxvectemp[1]=1.01*xxvecBL[1];
    xxvectemp[2]=1.0*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E1=pptemp[6];

    xxvectemp[1]=.99*xxvecBL[1];
    xxvectemp[2]=1.0*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E2=pptemp[6];

    Fx=(E2-E1)/(.02*xxvecBL[1]*sqrt(ggBL[1][1]))/chi/3.;

    //th dimension
    xxvectemp[1]=1.0*xxvecBL[1];
    xxvectemp[2]=1.01*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E1=pptemp[6];

    xxvectemp[1]=1.0*xxvecBL[1];
    xxvectemp[2]=0.99*xxvecBL[2];
    xxvectemp[3]=1.0*xxvecBL[3];
    calc_g_arb(xxvectemp,ggt,KERRCOORDS);
    calc_G_arb(xxvectemp,GGt,KERRCOORDS);

    anret=donut_analytical_solution(pptemp,xxvectemp,ggt,GGt);
    if(anret<0) anretmin=-1;
    E2=pptemp[6];

    Fy=(E2-E1)/(.02*xxvecBL[2]*sqrt(ggBL[2][2]))/chi/3.;

    //ph dimension - symmetry so not modified
    Fz=0.;

    if(anretmin<0) //one of the points outside the donut
      Fx=Fy=Fz=0.;
    else
      {
	ldouble Fl=sqrt(Fx*Fx+Fy*Fy+Fz*Fz);
	if(Fl>.99*E)
	  {
	    Fx=Fx/Fl*0.99*E;
	    Fy=Fy/Fl*0.99*E;
	    Fz=Fz/Fl*0.99*E;
	  }
      }

     //saving ff values to pp[]
     pp[7]=Fx;
     pp[8]=Fy;
     pp[9]=Fz;

#ifdef NOINITFLUX
     pp[7]=0.;
     pp[8]=0.;
     pp[9]=0.;
#endif


     //if(ix==NX-1 && iy==NY-1){print_Nvector(pp,NV);}

     //transforming from BL lab radiative primitives to code non-ortonormal primitives
     prad_ff2lab(pp,pp,&geomBL);
#endif
     //if(ix==NX-1 && iy==NY-1){print_Nvector(pp,NV);getchar();}
     

     //transforming primitives from BL to MYCOORDS
     trans_pall_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvecBL,ggBL,GGBL,gg,GG);
   }

//entropy
pp[5]=calc_Sfromu(pp[0],pp[1]);
//hd floors
//check_floors_hd(pp,VELPRIM,gg,GG);
//to conserved
p2u(pp,uu,&geom);

//test
//u2p_entropy(uu,pp,&geom);
//u2p_entropy_harm(uu,pp,&geom);

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


/*
  int
  set_initial_profile()
  {
  int ix,iy,iz;
  for(iz=0;iz<NZ;iz++)
  {
  for(iy=0;iy<NY;iy++)
  {
  for(ix=0;ix<NX;ix++)
  {
*/

ldouble rho,mx,my,mz,m,E,uint,Fx,Fy,Fz,pLTE;  
ldouble xx,yy,zz;
ldouble uu[NV],xxvec[4];

get_xx(ix,iy,iz,xxvec);
coco_N(xxvec,xxvec,MYCOORDS,BLCOORDS);
xx=xxvec[1];
yy=xxvec[2];
zz=xxvec[3];


ldouble gg[4][5],GG[4][5],eup[4][4],elo[4][4];
pick_g(ix,iy,iz,gg);
pick_G(ix,iy,iz,GG);
calc_ZAMOes(gg,eup,elo,MYCOORDS);

ldouble pp[NV],T;

//working in BL
ldouble ggBL[4][5],GGBL[4][5];
calc_g_arb(xxvec,ggBL,KERRCOORDS);
calc_G_arb(xxvec,GGBL,KERRCOORDS);
ldouble eupBL[4][4],eloBL[4][4];
ldouble tupBL[4][4],tloBL[4][4];
calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);

ldouble podpierd=-(GGBL[0][0]-2.*ELL*GGBL[0][3]+ELL*ELL*GGBL[3][3]);
ldouble ut=-1./sqrt(podpierd);

ut/=UTPOT; //rescales rin
ldouble Vphi,Vr;
ldouble D,W,eps,uT,uphi,uPhi;
if(ut<-1 || podpierd<0. || xx<3. || NODONUT || INFLOWING)
  {
    //ambient
    set_hdatmosphere(pp,xxvec,gg,GG,0);
#ifdef RADIATION
    set_radatmosphere(pp,xxvec,gg,GG,0);
#endif
  }
 else
   {
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
     conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
     trans2_coco(xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
     conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
   
     pp[2]=ucon[1]; 
     pp[3]=ucon[2];
     pp[4]=ucon[3];
     pp[0]=rho; pp[1]=uint; 

#ifdef RADIATION
     ldouble P,aaa,bbb;
     P=GAMMAM1*uint;
     //solving for T satisfying P=pgas+prad=bbb T + aaa T^4
     aaa=4.*SIGMA_RAD;
     bbb=K_BOLTZ*rho/MU_GAS/M_PROTON;
     ldouble naw1=cbrt(9*aaa*Power(bbb,2) - Sqrt(3)*Sqrt(27*Power(aaa,2)*Power(bbb,4) + 256*Power(aaa,3)*Power(P,3)));
     ldouble T4=-Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))/2. + Sqrt((4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 - naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa) + (2*bbb)/(aaa*Sqrt((-4*Power(0.6666666666666666,0.3333333333333333)*P)/naw1 + naw1/(Power(2,0.3333333333333333)*Power(3,0.6666666666666666)*aaa))))/2.;

     E=calc_LTE_EfromT(T4);
     Fx=Fy=Fz=0.;
     uint=calc_PEQ_ufromTrho(T4,rho);

     pp[1]=uint;
     pp[6]=E;
     pp[7]=Fx;
     pp[8]=Fy;
     pp[9]=Fz;

     //transforming BL ZAMO radiative primitives to BL non-ortonormal primitives
     ldouble eupBL[4][4],eloBL[4][4];
     ldouble tupBL[4][4],tloBL[4][4];
     calc_tetrades(ggBL,tupBL,tloBL,KERRCOORDS);
     calc_ZAMOes(ggBL,eupBL,eloBL,KERRCOORDS);
     prad_zamo2ff(pp,pp,ggBL,GGBL,eupBL);
     prad_ff2lab(pp,pp,ggBL,GGBL,tloBL);
     //transforming radiative primitives from BL to MYCOORDS
     trans_prad_coco(pp, pp, KERRCOORDS, MYCOORDS,xxvec,ggBL,GGBL,gg,GG);
#endif
   }


pp[5]=calc_Sfromu(pp[0],pp[1]);

//testing if interpolated primitives make sense
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

//if(isnan(get_u(p,5,ix,iy,iz))) {printf("pr: %d %d %d S: %Le\n",ix,iy,iz,0.);getchar();}

//mark initialy succesfull u2p_hot step
set_cflag(0,ix,iy,iz,0);

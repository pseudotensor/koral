
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

ldouble rho,mx,my,mz,m,E,uint,E0,Fx,Fy,Fz,pLTE;  
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

     //4-velocity in BL
     ldouble ucon[4]={0.,-Vr,0.,Vphi};
     conv_vels(ucon,ucon,VEL3,VEL4,ggBL,GGBL);
     trans2_coco(xxvec,ucon,ucon,BLCOORDS,MYCOORDS);
     conv_vels(ucon,ucon,VEL4,VELPRIM,gg,GG);
   
     pp[2]=ucon[1]; 
     pp[3]=ucon[2];
     pp[4]=ucon[3];

     //density etc.
     pp[0]=rho; pp[1]=uint; 
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
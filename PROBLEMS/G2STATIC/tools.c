//**********************************************************************
//**********************************************************************
//**********************************************************************
//sets rho, uint and velocities according to disk model from Sadowski+2012
int
set_sgradisk(ldouble *pp,ldouble *xx,void *ggg, void* gggBL)
{
  struct geometry *geom
    = (struct geometry *) ggg;

  struct geometry *geomBL
    = (struct geometry *) gggBL;

  // spherical coordinates
  ldouble xx2[4];
  coco_N(xx,xx2,MYCOORDS,BLCOORDS);
  ldouble r=xx2[1];
  ldouble th=xx2[2];

  //rotation here?

  //entropy = fixed at rin
  ldouble xxin[4]={0,get_x(0,0),get_x(NY/2,1),get_x(0,2)};
  coco_N(xxin,xxin,MYCOORDS,BLCOORDS);
  ldouble rin=xxin[1];

  rin=150.;
  ldouble rhoin = 2.02e5*pow(1.-pow((th-M_PI/2.)/(M_PI/2.),2.),1.69) * (150./rin);
  ldouble tempin = pow(10.,9.95+0.24*pow(fabs(th-M_PI/2.),2.93)) * (150./rin);
  //to code units
  rhoin = rhoCGS2GU(rhoin);
  tempin = tempCGS2GU(tempin);
  ldouble uin = calc_PEQ_ufromTrho(tempin,rhoin);
  ldouble K = (GAMMA-1.)*uin / pow(rhoin,GAMMA); //p = K rho^GAMMA

  //empirical fit in cgs
  ldouble rho = rhoin * (rin/r);

  //polytrope G=2
  rho = 1./2./K * 0.7 / r;

  ldouble temp = tempin * (150./r);
  ldouble vphi = pow(10.,9.15-0.24*pow(fabs(th-M_PI/2.),2.04)) * sqrt(150./r);

  //const.fraction of Keplerian:
  vphi=0.7*sqrt(1./r/r/r);

  ldouble chi = 0.1 + 0.31*pow(fabs(th-M_PI/2.),3.89); //pmag/pgas


  //to code units
  //rho = rhoCGS2GU(rho);
  temp = tempCGS2GU(temp);
  vphi = velCGS2GU(vphi);
  chi = chi;

  //temp following constant entropy: - actually should be satisfied already for GAMMA=2
  ldouble u = K * pow(rho,GAMMA) / (GAMMA-1.);
  temp = calc_PEQ_Tfromurho(u,rho); 
 
  pp[0] = rho;
  pp[1] = u;//calc_PEQ_ufromTrho(temp,rho);

  //to add extra magn-related pressure
  //pp[1] *= 1.+chi;

  ldouble ucon[4]={0.,0.,0.,vphi/r};

  conv_vels(ucon,ucon,VEL3,VEL4,geomBL->gg,geomBL->GG);
  trans2_coco(geomBL->xxvec,ucon,ucon,KERRCOORDS,MYCOORDS);
  conv_vels(ucon,ucon,VEL4,VELPRIM,geom->gg,geom->GG);

  pp[2]=ucon[1];
  pp[3]=ucon[2];
  pp[4]=ucon[3];

  return 0;

}

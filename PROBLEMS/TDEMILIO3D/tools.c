int find_globalindex(double r, double th, double ph, int gi[3])
{

  //in radius
    ldouble dist=BIG,distloc;
    int kbest=-1;
    ldouble rsph = r;
    int gix,giy,giz;
    dist=BIG;
    for(gix=0;gix<TNX;gix++)
      {
	int ix=gix-TOI;
	
	ldouble xx[4]={0.,0.5*(calc_xb(ix,0)+calc_xb(ix+1,0)),0.,0.};
	coco_N(xx,xx,MYCOORDS,BLCOORDS);
	ldouble rkoral=xx[1];


	//printf("%d %d %e %f\n",gix,ix,rkoral,rsph);
	
	distloc=fabs(rsph-rkoral);
	if(distloc<dist) {kbest=ix;dist=distloc;}
      }
    gi[0]=kbest;


    //in theta
    dist=BIG;
    kbest=-1;
   
    ldouble thsph = th;
    dist=BIG;
    for(giy=0;giy<TNY;giy++)
      {
	int iy=giy-TOJ;
	   
	ldouble xx[4]={0.,0.,0.5*(calc_xb(iy,1)+calc_xb(iy+1,1)),0.};
	coco_N(xx,xx,MYCOORDS,BLCOORDS);
	ldouble thkoral=xx[2];
	distloc=fabs(thsph-thkoral);
	if(distloc<dist) {kbest=iy;dist=distloc;}
      }
    gi[1]=kbest;

    //in phi
    dist=BIG;
    kbest=-1;
    ldouble phsph = ph;
    dist=BIG;
    for(giz=0;giz<TNZ;giz++)
      {
	int iz=giz-TOK;
	ldouble xx[4]={0.,0.,0.,0.5*(calc_xb(iz,2)+calc_xb(iz+1,2))};
	coco_N(xx,xx,MYCOORDS,BLCOORDS);
	ldouble phkoral=xx[3];
	distloc=fabs(phsph-phkoral);
	if(distloc<dist) {kbest=iz;dist=distloc;}
      }

    gi[2]=kbest;
   

    return 0;
}

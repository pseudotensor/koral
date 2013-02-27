//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

double minx,maxx,miny,maxy;


if(MYCOORDS==MKS1COORDS)
  {
    minx= (exp(get_xb( -1,0))+MKS1R0)*cos(get_xb(NZ+1,2));
    maxx= exp(get_xb(NX+1,0))+MKS1R0;
    miny=.5*sin(get_xb(-1,2))*(exp(get_xb(NX+1,0))+MKS1R0);
    maxy=.5*sin(get_xb(-1,2))*(exp(get_xb(NX+1,0))+MKS1R0)+( exp(get_xb(NX+1,0))+MKS1R0-((exp(get_xb( -1,0))+MKS1R0)*
											 cos(get_xb(NZ+1,2))));
  }
 else
   {
minx= get_xb( -1,0)*cos(get_xb(NZ+1,2));
maxx= get_xb(NX+1,0);
miny=.5*sin(get_xb(-1,2))*get_xb(NX+1,0);
maxy=.5*sin(get_xb(-1,2))*get_xb(NX+1,0)+( get_xb(NX+1,0)-get_xb( -1,0)*cos(get_xb(NZ+1,2)));
   }

  fprintf(fgnu,
	  "set view map\n"
	  "set pm3d\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 2\n"
	  //	  "set style arrow 1 head nofilled size screen 0.002,35 ls 1\n"



	  "set style arrow 1 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 35,3,9\n"



	  "unset surface\n"
	  "set term gif large size 600,500\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.8\n"
	  "set bmargin at screen .09\n"
	  "set tmargin at screen .92\n"
	  "set size .5,1\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"rad. energy density / flux \" offset 0,0\n"
#ifdef MINKOWSKI
	  "splot \"%s\" u 1:3:20 ti \"\"\n"
#else
	  //	  "set autoscale\n"
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):20 ti \"\" w l \n"
#endif

	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
#ifdef MINKOWSKI
	  "plot \"%s\" u 1:3:($21/(($21*$21+$22*$22+$23*$23)**.5)/%f):($23/(($21*$21+$22*$22+$23*$23)**.5)/%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname2,get_xb(0,0),get_xb(NX,0),get_xb(-NG,2),get_xb(NZ,2),fname,fname,
	  fname,30./(get_xb(NX,0)-get_xb(0,0)),30./(get_xb(NZ,2)-get_xb(0,2)),(int)(NX/20),(int)(NZ/20));
#else
//	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($21*cos($3)-($23)*sin($3))/(($21*cos($3)*$21*cos($3)+$23*sin($3)*$23*sin($3)+1.e-30*$20)**.5)/%f):(($23*cos($3)+($21)*sin($3))/(($21*cos($3)*$21*cos($3)+$23*sin($3)*$23*sin($3)+1.e-30*$20)**.5)/%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($21*cos($3)-$23*sin($3))/(%e)):(($23*cos($3)+$21*sin($3))/%e) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
  

,fname2,minx,maxx,miny,maxy,fname, fname,
#if(BEAMNO==1)
	    4.e-17,4.e-17,
#endif
#if(BEAMNO==2)
4.e-18,4.e-18,
#endif
#if(BEAMNO==3)
4.e-18,4.e-18,
#endif
(int)(NX/10),(int)(NZ/10)
);
#endif
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/

//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];


ldouble minx=0.;
ldouble miny=0.;
ldouble maxx=get_x(NX,0);
ldouble maxy=get_x(NY,0);


#ifdef myMKS1COORDS
minx= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxx= .66*(exp(get_xb(NX,0))+MKS1R0);
miny= -.02*(exp(get_xb(-NG,0))+MKS1R0);
maxy= .66*(exp(get_xb(NX,0))+MKS1R0);
#else
minx= -.02*get_xb(NX,0);
maxx= 1.02*get_xb(NX,0);
miny= -.02*get_xb(NX,0);
maxy= 1.02*get_xb(NX,0);
#endif

  fprintf(fgnu,
          "set table \"table.gp\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log cb\n"
	  "set cbrange [1.e-15:1e-11]\n"
	  "set cntrparam levels 20 \n"
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):20 ti \"\" w l ls 1\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

	  "set term gif large size 900,800\n"
	  "set output \"%s\"\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"

	  "set view map\n"
	  "set pm3d\n"
	  "unset surface\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 3 lc 3\n"
	  "set style line 11 lt 1 lw 2 lc 6\n"
	  "set style line 2 lt 1 lw 2 lc 2\n"
	  "set style line 3 lt 1 lw 2 lc 2\n"
	  "set style line 21 lt 3 lw 1 lc -1\n"
	  "set style line 1 lt 1 lw 2 lc 2\n"
	  "set style line 2 lt 1 lw 3 lc 3\n"
	 
	  "set style arrow 1 ls 1\n"
	  
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 23,28,3\n"
	  "set palette model RGB rgbformulae 7,8,9\n"
	  "set palette model RGB rgbformulae 35,3,9\n"
	  "set palette model RGB rgbformulae 6,3,21\n"
	 
	  "set autoscale\n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"

	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.85\n"
	  "set bmargin at screen .1\n"
	  "set tmargin at screen .95\n"
	  "set cblabel \"\"\n"
	  "set ylabel \"y\"\n"
	  "set title \"radiative E / flux\"\n"
	  "unset log cb\n"
	  "set log cb\n"
	  "set cbrange [1.e-13:1e-11]\n"
	  
	  "set palette model RGB rgbformulae 7,5,15\n"	 
	  "splot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):20 ti \"\" w l ls 1\n"

          "set isosam 10,10\n"
	  "set ylabel \"\"\n"
	  "unset tics\n"
	  "unset border\n"
	  "unset log cb\n"
	  "plot \"%s\" u (($1)*sin($2)):(($1)*cos($2)):(($21*sin($2)+$22*cos($2))/(($21*$21+$22*$22)**.5/%f*2)):((-$22*sin($2)+$21*cos($2))/(($21*$21+$22*$22)**.5/%f*2)) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  //	  "plot \"%s\" u 1:2:(($21)/(($21*$21+$22*$22)**.5/%f*2)):(($22)/(($21*$21+$22*$22)**.5/%f*2)) every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"table.gp\" w l ls 2\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"
	  ,fname
	  ,fname2,
	  minx,
	  maxx,
	  miny,
	  maxy,
	  
	  fname,
	  fname,1.,1.,NX/21+1,NY/21+1
	  );  
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/
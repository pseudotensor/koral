//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];

ldouble minx=0.;
ldouble miny=0.;
ldouble maxx=get_x(NX,0);
ldouble maxy=maxx;

#ifdef myMCYL1COORDS
maxx= 1.1*(exp(get_xb(NX,0))+MKS1R0);
maxy= maxx;
#endif


#ifdef FULLPHI
minx=-maxx;
miny=-maxx;
#endif

  fprintf(fgnu,
          "set table \"table.gp\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log z\n"
	  "set cntrparam levels 10 \n"
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):20 w l\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

	   "set term gif large size 800,700\n"
	   "set output \"%s\"\n"
	   "set view map\n"
	  "unset surface\n"
	  "set size 1,1\n"
	  "set origin 0,0\n"
	  "set multiplot\n"
	  "set pm3d\n"
	  "unset key\n"
	  "set style line 1 lt 1 lw 2 lc 2\n"
	  "set style line 2 lt 1 lw 3 lc 3\n"
	  "set style arrow 1 head nofilled size screen 0.002,35 ls 1\n"
	  "set palette model RGB rgbformulae 7,5,15\n"
	  "set palette model RGB rgbformulae 30,31,32\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  "set palette model RGB rgbformulae 35,3,9\n"
	  "set palette model RGB rgbformulae 21,22,23\n"
	  

	  "set lmargin at screen 0.1\n"
	  "set rmargin at screen 0.8\n"
	  "set bmargin at screen .12\n"
	  "set tmargin at screen .95\n"
	  "set autoscale \n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"
	  "set log cb\n"
	  "set cbrange [0.1:10.]\n"
	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
#ifdef MULTIRADFLUID
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):($20+$24+$28+$32) w l ti \"\"\n"
#else
	  "splot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):20 w l ti \"\"\n"
#endif

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  //"plot \"table.gp\" w l ls 2\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"

	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
#ifdef MULTIRADFLUID
	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):((($21+$25+$29+$33)*cos($3)-(($23+$27+$31+$35))*sin($3))*%f)"
	  ":((($23+$27+$31+$35)*cos($3)+(($21+$25+$29+$33))*sin($3))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname,fname2,minx,maxx,miny,maxy,fname,fname,
	  1.,1.,(int)(NX/15),(int)(NZ/20));
#else
	  "plot \"%s\" u (($1)*cos($3)):(($1)*sin($3)):(($21*cos($3)-($23)*sin($3))*%f)"
	  ":(($23*cos($3)+($21)*sin($3))*%f) every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname,fname2,minx,maxx,miny,maxy,fname,fname,
	  1.,1.,(int)(NX/15),(int)(NZ/20));
#endif	    
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/

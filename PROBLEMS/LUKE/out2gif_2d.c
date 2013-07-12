//int
//convert_out2gif_2d(char *fname,char *fname2,int niter,ldouble t)
//{
//  FILE *fgnu=fopen("plot.gp","w");
//  char bufor[50];
  fprintf(fgnu,
          

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
	  "set style arrow 1 ls 1\n"
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

	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"rho\" offset 0,-1\n"
	  "splot \"%s\" u 1:2:($15) w l ti \"\"\n"

	 
	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
	  "plot \"%s\" u 1:2:(($16)"
	  "/((($16)*($16)+($17)*($17))**.5)/%f):"
	  "(($17)/((($16)*($16)+($17)*($17))**.5)/%f)"
	  "every %d:%d w vectors arrowstyle 1 ti \"\"\n"
	  ,fname2,get_xb(-NG,0),get_xb(NX+NG,0),get_xb(-NG,1),get_xb(NY+NG,1),fname,fname,
	  30./(get_xb(NX,0)-get_xb(0,0)),30./(get_xb(NY,1)-get_xb(0,1)),(int)(NX/20),(int)(NY/20));
	    
	    

/*
  fprintf(fgnu,"\n");
  fclose(fgnu);   
  
  int i=system("gnuplot plot.gp ");
  return 0;
*/

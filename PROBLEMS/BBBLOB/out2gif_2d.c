char bufloc[100];
sprintf(bufloc,"table.gp.%d",PROCID);
 fprintf(fgnu,
          
	  "set table \"%s\"\n"
	  "set contour base\n"
	  "unset surface\n"
	  "set log z\n"
	  "set cntrparam levels discrete 1.1,3,6,10,30,60,100 \n"
	  "splot \"%s\" u 1:2:14 w l\n"
	  "unset dgrid3d\n"
	  "unset log z\n"
	  "unset table\n"
	  "unset contour\n"
	  "unset surface\n"

#ifndef RADIATION
	  "set term gif large size 800,700\n"
#else
	  "set term gif large size 1600,700\n"
#endif
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

#ifndef RADIATION
	  "set lmargin at screen 0.07\n"
	  "set rmargin at screen 0.85\n"
#else
	  "set lmargin at screen 0.05\n"
	  "set rmargin at screen 0.425\n"
#endif
	  "set bmargin at screen .12\n"
	  "set tmargin at screen .95\n"
	  "set autoscale \n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"
	  "set cbrange [1:2e10]\n"
	  //"set log cb\n"

	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"hydro\" offset 0,-1\n"
	  "splot \"%s\" u 1:2:($24) w l ti \"\"\n"
	 
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

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"%s\" w l ls 2\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"

#ifdef RADIATION
	  "set lmargin at screen 0.525\n"
	  "set rmargin at screen 0.90\n"
	  "set bmargin at screen .12\n"
	  "set tmargin at screen .95\n"
	  "set autoscale \n"
	  "set xrange [%f:%f]\n"
	  "set yrange [%f:%f]\n"
	  "set cbrange [1:2e10]\n"
	  //"set log cb\n"

	  "set xlabel \"x\"\n"
	  "set ylabel \"y\"\n"
	  "set cblabel \"\"\n"
	  "set title \"rad\" offset 0,-1\n"
	  "splot \"%s\" u 1:2:($25) w l ti \"\"\n"
	 
	  "unset pm3d\n"
	  "set isosam 10,10\n"
	  "set format x \"\"\n"
	  "set format y \"\"\n" 
	  "set xlabel \"\"\n"
	  "set ylabel \"\"\n"
	  "set title \"\" offset 0,-1\n"
	  "plot \"%s\" u 1:2:(($21)"
	  "/((($21)*($21)+($22)*($22))**.5)/%f):"
	  "(($22)/((($21)*($21)+($22)*($22))**.5)/%f)"
	  "every %d:%d w vectors arrowstyle 1 ti \"\"\n"

	  "unset tics\n"
	  "unset border\n"
	  "unset pm3d\n"
	  "unset surface\n"
	  "plot \"%s\" w l ls 2\n"
	  "set pm3d\n"
	  "set tics\n"
	  "set border\n"
#endif

	 ,bufloc,fname,fname2,get_xb(0,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1)
	 ,fname,fname,50./(get_xb(NX,0)-get_xb(0,0)),50./(get_xb(NY,1)-get_xb(0,1)),(int)(NX/20),(int)(NY/20),bufloc
#ifdef RADIATION
	  ,get_xb(0,0),get_xb(NX,0),get_xb(0,1),get_xb(NY,1)
	 ,fname,fname,50./(get_xb(NX,0)-get_xb(0,0)),50./(get_xb(NY,1)-get_xb(0,1)),(int)(NX/20),(int)(NY/20),bufloc
#endif
);
	    
	

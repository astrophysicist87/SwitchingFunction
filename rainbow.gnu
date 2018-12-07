FILES = system("ls -1 QCDPD*.dat")
TEMPS = system("ls -1 QCDPD*.dat | sed 's/QCDPD_circle_R_//g' | sed 's/.dat//g'")
plot for [i=1:words(FILES)] word(FILES, i) u 3:9 w l lc rgb hsv2rgb((word(TEMPS,i)-100.0)/200,1,1) lw 2

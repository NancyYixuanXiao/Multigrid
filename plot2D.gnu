#reset;
L=99;


#set isosample 5, 5
#set contour base
#set view map
#set dgrid3d;
#unset surface
#set cntrparam level incremental -1.1, 0.1, 1
#splot "T.txt" using 1:2:3 with lines




plot 'parallel_data.dat' using 1:2:3 with image

pause mouse;

!--------------------------------------------------------
!                             Kferimilevelplot
!--------------------------------------------------------
subroutine Kferimilevelplot (  )
use moduleone
!============  start   K_ferimilevelplot.plt'=======

open(unit =UNIT202, file = FILENAME202, status="unknown", action="write" , DELIM='none', iostat=ierror)

write(UNIT202,'(a)')  " set terminal pdf " 
write(UNIT202,'(a)')  " set output 'Ks_points0.pdf'   "
write(UNIT202,'(a)')  " set border linewidth 1.5    " 
write(UNIT202,'(a)')  " set pointsize 0.20  " 
write(UNIT202,'(a)')  " set style line 1 lc rgb ""red"" pt 7 ps 0.20   # circle    " 
write(UNIT202,'(a)')  " unset key  " 
write(UNIT202,'(a)')  " set tics scale 0.75  " 
write(UNIT202,'(a)')  " set xtics 1 " 
write(UNIT202,'(a)')  " set ytics 1  "
write(UNIT202,'(a)')  " set xlabel 'Kx'  "
write(UNIT202,'(a)')  " set ylabel 'Ky'  "
write(UNIT202,'(a)')  " m = ""./K_ferimilevel.dat""  "
write(UNIT202,'(a)')  "  plot m using 1:2 notitle w p ls 1 "
write(UNIT202,'(a)')  " quit " 
close(UNIT202)

call system('gnuplot -persist  K_ferimilevelplot.plt')
!~ call system('del *.plt')

end subroutine






!--------------------------------------------------------
!                             showgap1plot
!--------------------------------------------------------
subroutine showgap1plot (Y0, Y1, np, ua, uz  )
use moduleone
use moduletwo
implicit real (Kind=DP) (a-h, o-z)

real (KInd=DP), intent(in) :: Y0, Y1, np, ua, uz
character (len=7) :: material_name

!~ y0=0.2
!~ y1=1.2
!~ np=4.0
!~ ua=0.812
!~ uz=0.901


If( SetTransitionMetalDichalcogenides==1 ) then
material_name="MoS_{2}"
else If( SetTransitionMetalDichalcogenides==2 ) then
material_name="WS_{2}"
end if

open(unit =UNIT13, file = FILENAME13, status="unknown", action="write" , DELIM='none', iostat=ierror) !!  showgap1.plt

!~ write(UNIT13,'(a)')  "set term postscript eps enhanced color 24  " 
!~ write(UNIT13,'(a)')  "set output 'gap1.dat.eps'   "

write(UNIT13,'(a)')  "set terminal pdfcairo enhanced  " 
write(UNIT13,'(a)')  "set output 'gap1_dat.pdf'   "
write(UNIT13,'(a)')  "set style data linespoints   " 
write(UNIT13,'(a)')  "set style line 1 lc rgb ""blue"" pt 7 ps 0.5   # circle    " 
write(UNIT13,'(a)')  "set style line 2 lc rgb ""red"" pt 7 ps 0.5   # circle    " 
write(UNIT13,'(a)')  "set style line 3 lc rgb '#0060ad' pt 9 ps 0.5   # triangle    " 

write(UNIT13,  '(  3a, f6.2, a)'   ) "set key title"" ", material_name, "---- n'=", np,&
  &"{/Symbol \264} 10^{13} cm^{-2}"" at screen 0.8, 0.45  "

write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 1 from  187,", Y0, ", 0  to  187,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 2 from  296,", Y0, ", 0  to  296,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 3 from  513,", Y0, ", 0  to  513,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT13,'(a, f6.2, a, f6.2, a)')  "set arrow 4 from  730,", Y0, ", 0  to  730,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT13,'(a, f8.4, a, f8.4, a)')  "set arrow 5 from  0,", ua, ", 0  to 839,", ua, &
  &", 0 nohead back nofilled linecolor rgbcolor ""violet""  lw 1.5 dashtype 2  " 
write(UNIT13,'(a, f8.4, a, f8.4, a)')  "set arrow 6 from  0,", uz, ", 0  to 839,", uz, &
  &", 0 nohead back nofilled linecolor rgbcolor ""violet""  lw 1.5  " 

write(UNIT13,'(a)')  "set autoscale x   "
write(UNIT13,'(a)')  "set autoscale y   "
write(UNIT13,'(  2(a, f6.2), a )')  "set yrange [", y0, ":", y1, "]" 

write(UNIT13,'(a)')  "set xtics ( ""{/Symbol G}"" 0, ""M"" 187, ""K"" 296, &
  &""{/Symbol G}"" 513, ""K'"" 730, ""M"" 839)   font "",26""    "
write(UNIT13,'(a)')  "set ylabel ""{/Symbol=32 e }  (eV) ""   font "",28""  "
write(UNIT13,'(a)')  "plot ""gap1.dat""  using 1:2 notitle w p pt 8 ps 0.5 "
write(UNIT13,'(a)')  "quit " 
close(UNIT13)

call system('gnuplot -persist  showgap1.plt')
!~ call system('del *.plt')

end subroutine

!--------------------------------------------------------
!                             plot Bandstructure0.dat
!--------------------------------------------------------
subroutine showBandstructure0plot (Y0, Y1)
use moduleone
use moduletwo
implicit real (Kind=DP) (a-h, o-z)

real (KInd=DP), intent(in) :: Y0, Y1
character (len=7) :: material_name
!~ y0=0.2
!~ y1=1.2

If( SetTransitionMetalDichalcogenides==1 ) then
material_name="MoS_{2}"
else If( SetTransitionMetalDichalcogenides==2 ) then
material_name="WS_{2}"
end if

open(unit =UNIT14, file = FILENAME14, status="unknown", action="write" , DELIM='none', iostat=ierror) !!  'showBandstructure0.plt'

!~ write(UNIT13,'(a)')  "set term postscript eps enhanced color 24  " 
!~ write(UNIT1,'(a)')  "set output '***.eps'   "
write(UNIT14,'(a)')  "set terminal pdfcairo enhanced  " 
write(UNIT14,'(a)')  "set output 'Bandstructure0.pdf'   "
write(UNIT14,'(a)')  "set style data linespoints   " 
write(UNIT14,'(a)')  "set style line 1 lc rgb ""blue"" pt 7 ps 0.5   # circle    " 
write(UNIT14,'(a)')  "set style line 2 lc rgb ""red"" pt 7 ps 0.5   # circle    " 
write(UNIT14,'(a)')  "set style line 3 lc rgb '#0060ad' pt 9 ps 0.5   # triangle    " 

write(UNIT14,  '(  2a  )'   ) "set key title"" ", material_name

write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 1 from  187,", Y0, ", 0  to  187,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 2 from  296,", Y0, ", 0  to  296,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 3 from  513,", Y0, ", 0  to  513,", Y1, &
  & ", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT14,'(a, f6.2, a, f6.2, a)')  "set arrow 4 from  730,", Y0, ", 0  to  730,", Y1, &
  &", 0 nohead back nofilled linecolor rgbcolor ""black""  linewidth 1.500  " 
write(UNIT14,'(   a    )'                )  "set arrow 5 from  0, 0, 0  to 839, 0, 0 &
  &nohead back nofilled linecolor rgbcolor ""violet""  lw 1.5 dashtype 2  " 
write(UNIT14,'(a)')  "set autoscale x   "
write(UNIT14,'(a)')  "set autoscale y   "
write(UNIT14,'(  2(a, f6.2), a )')  "set yrange [", y0, ":", y1, "]" 

write(UNIT14,'(a)')  "set xtics ( ""{/Symbol G}"" 0, ""M"" 187, ""K"" 296, &
  &""{/Symbol G}"" 513, ""K'"" 730, ""M"" 839)   font "",26""    "
write(UNIT14,'(a)')  "set ylabel ""{/Symbol=32 e }  (eV) ""   font "",28""  "
write(UNIT14,'(a)')  "plot ""Bandstructure0.dat""  using 1:2 notitle w p pt 8 ps 0.5"    !! w p pt 8 ps 0.5
write(UNIT14,'(a)')  "quit " 
close(UNIT14)

call system('gnuplot -persist  showBandstructure0.plt')
!~ call system('del *.plt')

end subroutine

!--------------------------------------------------------
!                             
!--------------------------------------------------------







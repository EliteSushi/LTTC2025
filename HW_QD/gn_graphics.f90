subroutine initgnu(nframes, dt)
  implicit none
  integer :: nframes
  real(8) :: dt
  character(9) :: fname

  fname = 'plot.gn'
  open(1, file=fname, status='replace')

  ! Gnuplot script setup
  !write(1, *) "set key outside"
  write(1,*) "set terminal gif animate delay 10 size 600,800"
  write(1,*) "set output 'evo.gif'"
  write(1,*) "stats 'datafile' nooutput"
  write(1, *) "set xrange [-1:1]"
  write(1, *) "set yrange [0:200]"
  write(1, *) "set xlabel 'x'"
  write(1, *) "set ylabel 'Potential / {|Ψ|}^{2}'"
  write(1, *) "set title 'Wavefunction and Potential Evolution'"

  ! Animation loop using Gnuplot's multi-block indexing
  write(1, *) "do for [t=0:", nframes, "] {"
  write(1, *) '    set key title sprintf("Time: %.2f fs", t * ', dt, ')'
  write(1, *) "    plot 'psi.dat' index t using 1:2 with lines title 'Potential' lw 4, \"
  write(1, *) "         '' index t using 1:3 with lines title sprintf('Psi', t) lw 4"
  write(1, *) "}"

  close(1)
end subroutine initgnu

!-------------------------------!
subroutine graphgnu(dx, npoints, psi, t)
  !-------------------------------!
  use parameters
  implicit none

  integer :: i, j, npoints
  real(8) :: x, dx, psi(npoints), t, potential,potscale,cutoff, psiwr

  if (potentialtype=='harmonic') then
      potscale=20.d0
      cutoff=16.d0
  elseif (potentialtype=='doublewell') then
      potscale=40.d0
      cutoff=8.d0
  endif

  write(12,*) ""
  write(12,*) ""
  write (12, *) "# x          pot          psi       (t = ", t, ")"
  do i=-npoints/2+1,npoints/2
      x=dble(i)*dx
      if (potentialtype=='harmonic') then
        potential=0.5d0*mass*angfreq**2*x**2*au2kcalmol*potscale
      elseif (potentialtype=='doublewell') then
        potential=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)*potscale
      endif
      if (i>0) then
        j=i
      else     
        j=i+npoints
      endif
      if (abs(x)<=1.d0) then
          psiwr = 100.d0*psi(j)
          write(12, '(F10.5, 1X, F10.5, 1X, F10.5)') x, potential, psiwr
    endif
  end do

end subroutine graphgnu
!-----------------------!
subroutine gnucoeff(ncoeff)
  implicit none
  integer :: ncoeff, i
  character(9) :: fname

  fname = 'coeff.gn'
  open(1, file=fname, status='replace')

  write(1,*) "set terminal png size 800,600"
  write(1,*) "set output 'coeff.png'"
  write(1,*) "stats 'coeff.dat' nooutput"
  write(1, *) "set xlabel 't (fs)'"
  write(1, *) "set ylabel 'Real part of Coefficient'"
  write(1, *) "set title 'Wavefunction Coefficient Time Evolution'"

  ! Animation loop using Gnuplot's multi-block indexing
  write(1, *) "p 'coeff.dat' u 1:2 w l lw 3 t 'c0', \"
  do i=2, ncoeff-1
    write(1, '(A, I0, A, I0, A)') " 'coeff.dat' u 1:", i+1, " w l lw 3 t 'c", i-1, "', \"
  end do
  write(1, '(A, I0, A, I0, A)') " 'coeff.dat' u 1:", i+1, " w l lw 3 t 'c", i-1, "'"

  close(1)
end subroutine gnucoeff
!-----------------!
 module parameters
!-----------------!

 real(8), parameter :: length=5.12d0   ! length of the box (in Bohr)
 real(8), parameter :: mass = 1822.88839d0 != 1 atomic mass unit (=1 g/mol)
 real(8), parameter :: pi=3.141592653589793d0
 real(8), parameter :: au2kcalmol=627.509d0
 real(8), parameter :: fs2au=41.341373336561d0
 real(8) :: angfreq, barrier
 character(10) :: potentialtype

 end module parameters
!---------------------!

!-----------------!
module functions2use
!-----------------!
   implicit none
   contains

   !---------------------------------!
   recursive function hermite(n,y) result(h) !Recursive function to calculate the hermite polynomials
      implicit none

      integer, intent(in) :: n
      real(8), intent(in) :: y
      real(8) :: h

      if (n<0) then !Check for proper usage
         write(*,*) "Error in hermite: n must be <=0"
         h=0.d0
      endif

      if (n==0) then !If statements to end recursion
         h=1.d0
      elseif (n==1) then
         h=2.d0*y
      else
         h=2.d0*y*hermite(n-1,y)-2.d0*(n-1)*hermite(n-2,y) !Recursion
      endif
   end function hermite
   !----------------------!

   !---------------------------------!
   recursive function factorial(n) result(fact) !Recursive function to calculate the factorial
      implicit none

      integer :: n
      integer :: fact

      if (n<0) then
         write(*,*) "Error in factorial: n must be greater than or equal to 0"
         fact=0
      endif
      
      if (n==0) then
         fact=1
      else
         fact=n*factorial(n-1)
      endif
   end function factorial
   !----------------------!
end module functions2use
!-----------------!

!-----------------!
 program propagate
!-----------------!
 use parameters
 use functions2use
 implicit none

 integer :: npoints,ntime,snapshot,i,ncoeff,iostat
 real(8) :: alpha,dt,t,dx,x0
 real(8), allocatable :: pot(:),kin(:),psisquare(:),coeff(:)
 complex(8), allocatable :: psi(:),psi0(:),exppot(:),expkin(:)


 open(unit=10,file='wavepacket') !, status='old', action='read', iostat=iostat)
   read(10,*) npoints               !Number of lattice points
   read(10,*) x0                    !Initial position
   read(10,*) alpha                 !Governs the initial width of the wave packet
   read(10,*) dt                    !Propagation time step
   read(10,*) ntime                 !Number of propagation steps
   read(10,*) snapshot              !snapshot frequency
   read(10,*) ncoeff                !Number of coefficients
   allocate(coeff(ncoeff))
   coeff=0.d0
   do i = 1, ncoeff                 !Read ncoeff coefficients or until EOF 
      read(10,*,iostat=iostat) coeff(i)
      if (iostat /= 0) exit 
   end do
 close(10)

 open(unit=11,file='potential')
   read(11,*) potentialtype         !harmonic or double well potential
   read(11,*) angfreq               !Angular frequency for harmonic potential
   read(11,*) barrier               !Height of barrier in double well potential (in kcal/mol)
 close(11)

 dt=dt*fs2au                        !convert femtoseconds to atomic units
 angfreq=angfreq/fs2au              !convert femtoseconds to atomic units

 allocate(psi(npoints),psi0(npoints))
 allocate(pot(npoints),exppot(npoints))
 allocate(kin(npoints),expkin(npoints))
 allocate(psisquare(npoints))

 dx=length/dble(npoints)

 call initpsi(npoints,dx,x0,psi0,ncoeff,coeff)             !Obtain initial wavepacket psi0
 call fourier(0,npoints,psi0)                        !Initialize the FFT
 call operators(npoints,dx,dt,pot,kin,exppot,expkin) !Calculate the kinetic and potential operators

 open(12, file="psi.dat", status='replace')
 open(14, file="coeff.dat", status='replace')  
 psi=psi0                                            !Set the wavepacket psi at t=0 equal psi0
 do i=0,ntime                                        !Start propagation
    t=i*dt
    if (i>0) then
       psi=psi*exppot                                !Multiply psi with exp(-i*dt*potential)
       call fourier(1,npoints,psi)                   !Forward FFT to momentum space
       psi=psi*expkin                                !Multiply psi with the exp(-i*dt*kinetic operator)
       call fourier(-1,npoints,psi)                  !Backward FFT to position space
    endif
    if (mod(i,snapshot)==0) then                     !Take a snapshot if the remainder of i/snapshot equals 0
       !call initgraph(i/snapshot,t)                  !Initialize graph
       psisquare=(abs(psi))**2
       !call graphpot(dx,npoints)                     !Plot the potential
       !call graphpsi(dx,npoints,psisquare)           !Plot |psi|^2
       call graphgnu(dx,npoints,psisquare,t)          !Plot potential, |psi|^2
       !call getcs(npoints,dx,x0,psi,ncoeff,coeff,t) !Get coeffs
    endif
 end do                                              !End propagation

 close(12)
 call initgnu(ntime/snapshot, dt*dble(snapshot)/fs2au)                                  !Initialize gnuplot
 call gnucoeff(ncoeff)
 call system('gnuplot plot.gn')                      !Execute gnuplot
 !call system('gnuplot coeff.gn')                      !Execute gnuplot

 deallocate(psi,psi0)
 deallocate(pot,exppot)
 deallocate(kin,expkin)
 deallocate(psisquare)

 end program propagate
!---------------------!

!------------------------------------------------!
 subroutine initpsi(npoints,dx,x0,psi0,ncoeff,coeff)      !Modified initpsi to initiate wavefunction expanded in the harmonic eigenfunctions
 use parameters
 use functions2use
 implicit none

 integer :: i,j,l,npoints,ncoeff
 real(8) :: x,x0,dx,mo,sqmopi,sqmo,A,coeff(ncoeff), N
 complex(8) :: psi0(npoints), res
 
 mo = mass*angfreq !Constants that will be needed in the loop
 sqmo = sqrt(mo) 
 sqmopi = sqmo/sqrt(pi)
 psi0 = (0.d0,0.d0) 
do l = 0, ncoeff-1 !Loop through the user specified eigenfunctions
   if (coeff(l+1)==0) cycle !if Coefficient = 0 skip
   A = sqrt(1.d0/(2.d0**l*dble(gamma(real(l+1))))*sqmopi) !Constant of the Eigenfunction
   do i=-npoints/2+1,npoints/2
      x=dble(i)*dx
      if (i>0) then
         j=i
      else     
         j=i+npoints
      endif
      res = coeff(l+1)*A*exp(-mo*(x-x0)**2/2.d0)*hermite(l,(x-x0)*sqmo) !coefficient*Eigenfucntion
      psi0(j) = psi0(j) + res
   end do
enddo

 N = sqrt(sum(coeff**2))!Normalize eigenfunction
 psi0 = psi0 / N 
 end subroutine initpsi
!----------------------!

!------------------------------------------------------!
 subroutine operators(npoints,dx,dt,pot,kin,exppot,expkin)
!------------------------------------------------------!
 use parameters
 implicit none

 integer :: i,j,npoints
 real(8) :: x,p,b,dt,dx,dp,pot(npoints),kin(npoints)
 complex(8) :: exppot(npoints),expkin(npoints)

 dp=2.d0*pi/length
 do i=-npoints/2+1,npoints/2
    x=dble(i)*dx
    p=dble(i-1)*dp
    if (i>0) then
       j=i
    else
       j=i+npoints
    endif
    if (potentialtype=='harmonic') then
       pot(j)=0.5d0*mass*angfreq**2*x**2
    elseif (potentialtype=='doublewell') then
       pot(j)=barrier*(16.d0*x**4 - 8.d0*x**2 + 1.d0)/au2kcalmol
    endif
    kin(j)=0.5d0*p**2/mass
    exppot(j)=exp(-dt*(0,1)*pot(j))
    expkin(j)=exp(-dt*(0,1)*kin(j))
 end do

 end subroutine operators
!------------------------!

!-----------------------------------!
 subroutine fourier(dir,npoints,psi)
!-----------------------------------!
 implicit none

 integer :: i,npoints,dir
 real(8) :: nr
 complex(8) :: psi(npoints)
 real(8), allocatable, save :: wsave(:)       

 if (dir==1) then
    call dcfftf(npoints,psi,wsave)
    nr=1.d0/dble(npoints)
    do i=1,npoints
       psi(i)=psi(i)*nr
    end do
 elseif (dir==-1) then
    call dcfftb(npoints,psi,wsave)
 elseif (dir==0) then
    if (allocated(wsave)) deallocate(wsave)
    allocate(wsave(4*npoints+20))
    call dcffti(npoints,wsave)
 endif

 end subroutine fourier
!----------------------!

!------------------------------------------------!
 subroutine getcs(npoints,dx,x0,psi,ncoeff,coeff,t)      !Modified initpsi to initiate wavefunction expanded in the harmonic eigenfunctions
!------------------------------------------------!
 use parameters
 use functions2use
 implicit none

 integer :: i,j,l,npoints,ncoeff
 real(8) :: x,x0,dx,mo,sqmopi,sqmo,A,coeff(ncoeff),t
 complex(8) :: psi(npoints), res
 
mo = mass*angfreq
 sqmo = sqrt(mass*angfreq)
 sqmopi = sqmo/sqrt(pi)
 coeff = 0.d0
   do l = 0, ncoeff-1
      A = sqrt(1.d0/(2.d0**l*dble(gamma(real(l+1))))*sqmopi)
      do i=-npoints/2+1,npoints/2
         x=dble(i)*dx
         if (i>0) then
            j=i
         else     
            j=i+npoints
         endif
         res = A*exp(-mo*(x-x0)**2/2.d0)*hermite(l,(x-x0)*sqmo)
         coeff(l+1) = coeff(l+1) + real(conjg(psi(j)))*res*dx
      end do
   enddo
   write(14,*) t/fs2au, coeff
 end subroutine getcs
!----------------------!
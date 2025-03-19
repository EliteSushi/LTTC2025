module utilities

public::lecture_xyz             ! Read and save positions in an .xyz file
public::get_distances           ! Calculate square of distances between atoms
public::get_forces              ! Calculate forces from the positions
public::get_pot_ener            ! Calculate energy from the positions
public::output_structure        ! Write output structures 
public::output_energy           ! Write output energies
public::get_kin_ener            ! Calculate kinetic energy from velocities
public::get_T                   ! Calculate T from kinetic energy
public::thermostat               ! Apply thermostat
public::get_com                 ! Calculate Center of Mass

contains

subroutine lecture_xyz (filename,tableau,nbr_atomes)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Record frames from an .xyz file   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

implicit none

character*100, intent(in)  :: filename
integer,       intent(out) :: nbr_atomes

real*8, allocatable, dimension (:,:),intent(out) :: tableau

character*10 :: var_poub
character*10 :: sep
integer      :: OK
integer      :: i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(1,file=filename,iostat=OK)

read(1,*,iostat=OK), nbr_atomes
read(1,*,iostat=OK), sep

allocate(tableau(nbr_atomes,3))

do i=1,nbr_atomes
  read(1,*,iostat=OK), var_poub, tableau(i,1), tableau(i,2), tableau(i,3)
enddo

close(1)

endsubroutine lecture_xyz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate square of distances between atoms !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_distances (positions,distances,nbr_atomes)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: positions
real*8, allocatable, dimension (:,:),intent(inout) :: distances

integer, intent(in)  :: nbr_atomes
integer              :: i,j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!do i=1,nbr_atomes
!  do j=1,nbr_atomes
!   distances(i,j)=      (positions(i,1)-positions(j,1))**2 + &
!					&   (positions(i,2)-positions(j,2))**2 + &
!					&	(positions(i,3)-positions(j,3))**2
!  enddo
!enddo

do i=1,nbr_atomes-1
  do j=i+1,nbr_atomes
   distances(i,j)= (positions(i,1)-positions(j,1))**2 + &
					&        (positions(i,2)-positions(j,2))**2 + &
					&	       (positions(i,3)-positions(j,3))**2
   distances(j,i)=distances(i,j)
  enddo
enddo

endsubroutine get_distances

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ forces               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_forces (forces,distances,positions,nbr_atomes,sigma,eps)

implicit none

real*8, allocatable, dimension (:,:),intent(inout) :: forces
real*8, allocatable, dimension (:,:),intent(in)    :: distances
real*8, allocatable, dimension (:,:),intent(in)    :: positions

real*8,intent(in)  :: sigma,eps
real*8             :: sigma_12,sigma_6
real*8             :: dist2
real*8             :: dr4,dr8,dr14
real*8             :: part_6,part_12
integer,intent(in) :: nbr_atomes
integer            :: i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,nbr_atomes
  do j=1,3
    forces(i,j)=0.0d0
  enddo
enddo

sigma_6=sigma**6
sigma_12=sigma**12

do i=1,nbr_atomes-1
  do k=i+1,nbr_atomes
	do j=1,3
	  dist2=(positions(i,j)-positions(k,j)) ! this will not comply with the PBC
	  dr4=distances(i,k)*distances(i,k)
	  dr8=dr4*dr4
	  dr14=dr8*dr4*distances(i,k)
	  part_6  = sigma_6*6.0d0*(-(dist2)/dr8)
	  part_12 = sigma_12*6.d0*(-(2.0d0*dist2)/dr14) 
	  forces(i,j)=forces(i,j)-4.0d0*eps*(part_12-part_6)
	  forces(k,j)=forces(k,j)+4.0d0*eps*(part_12-part_6)
	enddo
  enddo
enddo

endsubroutine get_forces

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ Energy               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_pot_ener (distances,nbr_atomes,sigma,eps,energy)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: distances

real*8,intent(in)   :: sigma,eps
real*8              :: sigma_6,sigma_12
real*8              :: part_6,part_12
real*8              :: dr6,dr12
real*8,intent(inout):: energy
integer,intent(in)  :: nbr_atomes
integer             :: i,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

energy=0.0d0

sigma_6=sigma**6
sigma_12=sigma**12

do i=1,nbr_atomes-1
    do k=i+1,nbr_atomes
	  dr6=distances(i,k)*distances(i,k)*distances(i,k)
	  dr12=dr6*dr6
	  part_6  = sigma_6/dr6
	  part_12 = sigma_12/dr12
      energy=energy+4*eps*(part_12-part_6)
	enddo
enddo

endsubroutine get_pot_ener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output gemetry                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_structure (index,positions,nbr_atomes,energy_pot)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: positions

real*8,intent(in)   :: energy_pot
integer,intent(in)  :: index
integer,intent(in)  :: nbr_atomes
integer             :: i,j
real*8, parameter   :: bohr_to_ang=0.5291772108d0 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

write(index,*), nbr_atomes
write(index,*), energy_pot
do i=1,nbr_atomes
  write(index,*), "Ar", positions(i,1)*bohr_to_ang, positions(i,2)*bohr_to_ang, positions(i,3)*bohr_to_ang
enddo

endsubroutine output_structure

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Ouput Energies                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_energy (index,index_loop,timestep,energy_pot,energy_kin,temp)

implicit none

real*8, parameter   :: fs_to_timeau=41.34137314d0
real*8, parameter   :: Hartree_to_kcal=627.50947d0
real*8,intent(in)   :: energy_pot,energy_kin,temp
integer,intent(in)  :: index,index_loop
real*8,intent(in)   :: timestep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (index_loop == 1) then
  write(index, '(A10, 2x, A14, 2x, A16, 2x, A16)') "Timestep", "Temp", "Potential Energy", "Kinetic Energy"
  write(index, '(f10.2, 2x, f14.10, 2x, f16.10, 2x, f16.10)') timestep * index_loop * 1.0d0 / fs_to_timeau, temp, &
                                                             energy_pot * Hartree_to_kcal, energy_kin * Hartree_to_kcal
else
  write(index, '(f10.2, 2x, f14.10, 2x, f16.10, 2x, f16.10)') timestep * index_loop * 1.0d0 / fs_to_timeau, temp, &
                                                             energy_pot * Hartree_to_kcal, energy_kin * Hartree_to_kcal
endif

endsubroutine output_energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Kinetic Energy          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_kin_ener (velo,nbr_atomes,kin_ener,mass)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: velo

real*8,intent(inout):: kin_ener
real*8,intent(in)   :: mass
integer,intent(in)  :: nbr_atomes
integer             :: i,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

kin_ener=0.0d0

do i=1,nbr_atomes
 do k=1,3
   kin_ener=kin_ener+mass*velo(i,k)*velo(i,k)
 enddo
enddo

kin_ener=0.5d0*kin_ener

endsubroutine get_kin_ener

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate Instantaneous T         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_T (kin_ener,temp,nbr_atomes)

implicit none

real*8,intent(inout) :: temp
real*8,intent(in)    :: kin_ener
integer,intent(in)   :: nbr_atomes
real*8, parameter    :: kb=3.166811d-6   ! Hartree/K

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

temp=2.0d0*kin_ener/(kb*3.d0*dble(nbr_atomes))

endsubroutine get_T


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Apply Thermostat                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine thermostat(velo,thermo,temp,temp_0,val,nbr_atomes,timestep)

implicit none

real*8, allocatable, dimension (:,:),intent(inout)    :: velo

real*8               :: lambda     ! Scaling parameter of the thermostat
real*8, intent(in)   :: temp,temp_0,val,timestep
integer, intent(in)  :: thermo,nbr_atomes              
integer              :: i,j  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of definition Variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(thermo==1) then   ! Rescaling Thermostat

  lambda=sqrt(temp_0/temp)

  do i=1,nbr_atomes
    do j=1,3
      velo(i,j)=velo(i,j)*lambda
    enddo
  enddo

elseif (thermo==2) then ! Berendsen Thermostat


  lambda=1.0d0+(timestep/val)*((temp_0/temp)-1.0d0)
  lambda=sqrt(lambda)
  
  do i=1,nbr_atomes
    do j=1,3
      velo(i,j)=velo(i,j)*lambda
    enddo
  enddo

endif

endsubroutine thermostat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Gauss Distribution               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real*8 function gaussian_distr()
implicit none
integer, save :: iset=0
real*8 :: fac,rsq,var_1,var_2,x,y
real*8, save :: gset
call random_seed()
if (iset==0) then
! first time to calculate rsq
call random_number(x)
call random_number(y)
var_1 = 2.d0*x-1.d0
var_2 = 2.d0*y-1.d0
rsq = var_1**2+var_2**2
do while (rsq.ge.1..or.rsq.eq.0.) ! repeat the operation until rsq is between 0 and 1
call random_number(x)
call random_number(y)
var_1 = 2.d0*x-1.d0
var_2 = 2.d0*y-1.d0
rsq = var_1**2+var_2**2
enddo
fac = sqrt(-2.d0*log(rsq)/rsq)
gset = var_1*fac
gaussian_distr = var_2*fac
iset = 1
else
gaussian_distr = gset
iset = 0
endif
end function gaussian_distr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Spherical confinement potential (scp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine apply_scp(forces,pot_ener,positions,nbr_atomes,R_scp)
  implicit none

  integer,intent(in)                                    :: nbr_atomes
  real*8, intent(in)                                    :: R_scp
  real*8, allocatable, dimension (:,:),intent(inout)       :: positions
  real*8, allocatable, dimension (:,:),intent(inout)    :: forces
  real*8,intent(inout)                                  :: pot_ener
  
  real*8, dimension(3)  :: com
  real*8  :: d_com
  integer :: i,j

  !!!!!! Get center of Mass !!!!!
  com = 0.0d0
  do i=1,nbr_atomes
    do j=1,3
      com(j)=com(j)+positions(i,j) !if we had diffrent masses we would have to "mass-avarage"
    enddo
  enddo
  com=com/dble(nbr_atomes) 

  !!!!!! Apply SCP !!!!!
  do i=1,nbr_atomes
    d_com = 0.0d0
    do j=1,3
      positions(i,j) = positions(i,j) - com(j) ! Center the system
      d_com = d_com + (positions(i,j))**2 ! Distance to origin (COM)
    enddo
    d_com = sqrt(d_com)-R_scp
    if (d_com .gt. 0.0d0) then
      pot_ener = pot_ener + 0.008d0*(d_com)**4 ! update potential energy
      do j=1,3
        forces(i,j) = forces(i,j) -0.032d0*(d_com)**3*(positions(i,j))/(d_com+R_scp) !update forces calculateed previously
      enddo
    end if
  enddo
end subroutine apply_SCP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bolztmann distribution for initial velocities !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initial_velocities(velocities,nbr_atomes,temp,mass)
  implicit none

  real*8, dimension (:,:),intent(inout) :: velocities 
  real*8, intent(in) :: temp,mass
  integer, intent(in) :: nbr_atomes
  real*8 :: kb=3.166811d-6
  integer :: i,j

  do i=1,nbr_atomes
    do j=1,3
      velocities(i,j)=0.0d0
      velocities(i,j)=sqrt(kb*temp/mass)*gaussian_distr() ! Initial velocities following the Boltzmann distribution
    enddo
  enddo

end subroutine initial_velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate distance in PBC (closest image of evry pair)!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_distances_PBC(positions,distances,nbr_atomes,L)

implicit none

real*8, allocatable, dimension (:,:),intent(in)    :: positions
real*8, allocatable, dimension (:,:),intent(inout) :: distances
real*8, intent(in)  :: L
integer, intent(in)  :: nbr_atomes

real*8              :: dist
integer              :: i,j,k

distances=0.0d0
do i=1,nbr_atomes-1
  do j=i+1,nbr_atomes
    do k=1,3
      dist = positions(i,k)-positions(j,k)
      distances(i,j)=distances(i,j) + (dist - L*dble(nint(dist/L)))**2 !We use the nearest image of each particle pair
    enddo 
   distances(j,i)=distances(i,j)
  enddo
enddo

endsubroutine get_distances_PBC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Calculate LJ distances and forces  in PBC!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_distances_forces_PBC(positions,distances,forces,nbr_atomes,L,sigma,eps)

implicit none

real*8, allocatable, dimension (:,:),intent(inout) :: positions
real*8, allocatable, dimension (:,:),intent(inout) :: distances
real*8, allocatable, dimension (:,:),intent(inout) :: forces

real*8, intent(in)   :: L, sigma, eps
integer, intent(in)  :: nbr_atomes

real*8             :: dist
real*8             :: sigma_12,sigma_6
real*8             :: dist2
real*8             :: dr4,dr8,dr14
real*8             :: part_6,part_12
integer            :: i,j,k

forces=0.0d0
distances=0.0d0

do i=1,nbr_atomes-1 !calculate distances**2
  do j=i+1,nbr_atomes
    do k=1,3
      dist = positions(i,k)-positions(j,k)
      distances(i,j)=distances(i,j) + (dist - L*dble(nint(dist/L)))**2 !We use the nearest image of each particle pair
    enddo 
   distances(j,i)=distances(i,j)
  enddo
enddo

sigma_6=sigma**6*6.0d0 !calculate forces
sigma_12=sigma**12*6.0d0

do i=1,nbr_atomes-1
  do k=i+1,nbr_atomes
    do j=1,3
      dist2=positions(i,j)-positions(k,j) ! Need dist component for the F direction
      dist2=dist2-L*dble(nint(dist2/L)) ! closest image for PBC
      dr4=distances(i,k)*distances(i,k)
      dr8=dr4*dr4
      dr14=dr8*dr4*distances(i,k)
      part_6  = sigma_6*(-(dist2)/dr8)
      part_12 = sigma_12*(-(2.0d0*dist2)/dr14) 
      forces(i,j)=forces(i,j)-4.0d0*eps*(part_12-part_6)
      forces(k,j)=forces(k,j)+4.0d0*eps*(part_12-part_6)
    enddo
  enddo
enddo
endsubroutine get_distances_forces_PBC

subroutine wrap_atoms_PBC(positions, nbr_atomes, L)
  implicit none
  real*8, allocatable, dimension (:,:),intent(inout) :: positions
  real*8, intent(in)  :: L
  integer, intent(in)  :: nbr_atomes
  integer            :: i,j

  do i=1,nbr_atomes
    do j=1,3
      positions(i,j) = positions(i,j) - L*dble(nint(positions(i,j)/L))
    enddo
  enddo
end subroutine wrap_atoms_PBC

end module utilities

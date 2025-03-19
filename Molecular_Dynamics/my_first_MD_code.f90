program my_first_MD_Code

use utilities

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Variable Déclaration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real*8, parameter  :: pi=3.141592653589
real*8, parameter  :: kb=3.166811563d-6 ! Boltzmann constant in Hartree/K
real*8, parameter  :: kcal_to_Hartree=0.0015936015d0
real*8, parameter  :: Hartree_to_kcal=627.50947d0
real*8, parameter  :: ang_to_bohr=1.8897261d0
real*8, parameter  :: bohr_to_ang=0.52917721d0
real*8, parameter  :: fs_to_timeau=41.34137314d0
real*8, parameter  :: amu_to_au=1822.8884850d0
real*8, parameter  :: au_to_amu=0.00054857990943d0

real*8             :: mass=40.0d0*amu_to_au ! Mass of the particle
real*8             :: timestep              ! MD Timestep in fs
real*8             :: energy_pot            ! Potential energy of the system
real*8             :: energy_kin            ! Kinetic energy of the system
real*8             :: sigma=2.0d0*ang_to_bohr     ! Definition of the LJ potential
real*8             :: eps=1.0d0*kcal_to_Hartree   ! Definition of the LJ potential
real*8             :: temp                        ! Instantaneous T°
real*8             :: temp_0=100.0d0              ! Simulation T°
real*8             :: thermo_val=100000.0d0*fs_to_timeau  ! Strenght of the Thermostat
real*8             :: R_scp=10.0d0*ang_to_bohr              ! scp Radius
real*8             :: box_L                                 ! Box size

integer            :: max_md_Step=100       ! Maximum number of MD steps
integer            :: str_md_out=1          ! Number of steps between writing structure
integer            :: i,j,k,l,m             ! Indexes
integer            :: thermo=0              ! To define the thermostat 
integer            :: nbr_atomes            ! Number of atom in the simulation
integer            :: save_freq=1           ! Frequency of saving the data

logical           :: scp=.FALSE.           ! Activate the scp
logical           :: pbc=.FALSE.           ! Activate the PBC

real*8, allocatable, dimension(:,:) :: distances
real*8, allocatable, dimension(:,:) :: positions
real*8, allocatable, dimension(:,:) :: forces
real*8, allocatable, dimension(:,:) :: velocities

character(len=100) :: string_lec        ! To read input data
character(len=100) :: filename="input.xyz" ! Input file name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of Variable Declaration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (COMMAND_ARGUMENT_COUNT().eq.0) then
		write(*,*) "This a kind of Manual for your program!"
        write(*,*) "OPTIONS: "
        write(*,*) " -FILE           Filename (default is input.xyz)"
        write(*,*) " -TIMESTEP       MD Timestep in fs"
        write(*,*) " -MAX_MD_STEPS   Maximum number of MD steps (default is 100)"
        write(*,*) " -EPS            Energy definition of the LJ potential in kcal.mol-1"
        write(*,*) " -SIGMA          Distance definition of the LJ potential in A"
		    write(*,*) " -MASS           Mass in a.m.u. of the particle"
		    write(*,*) " -TEMP           Temperature of the Simulation in K"
        write(*,*) " -THERMOSTAT     Define the thermostats : SCALE or BEREN, if BEREN, &
		                                 &the coupling parameter should follow BEREN"
		    write(*,*) " --SAVEFREQ      Frequency of saving the data"
        write(*,*) " -SCP            Activate spherical confinement potential by providing &
                                     &the threshold radius (in A). Cannot be used with -PBC"
        write(*,*) " -PBC            Activate Periodic Boundry Conditions by providing &
                                     &the box length (in A). Cannot be used with -SCP"                             
        stop
endif

scp = .FALSE. ! SCP and PBC are deafaulted to inactive.
pbc = .FALSE.

do i = 1, COMMAND_ARGUMENT_COUNT()

    call getarg(i,string_lec)
    if (index(string_lec,'-FILE').NE.0)then
      call getarg(i+1,filename)
    endif
    if (index(string_lec,'-TIMESTEP').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) timestep
      timestep=timestep*fs_to_timeau
    endif
    if (index(string_lec,'-MAX_MD_STEPS').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) max_md_Step
    endif
    if (index(string_lec,'-EPS').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) eps
      eps=eps*kcal_to_Hartree
    endif
    if (index(string_lec,'-SIGMA').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) sigma
      sigma=sigma*ang_to_bohr
    endif
    if (index(string_lec,'-TEMP').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) temp_0
    endif
    if (index(string_lec,'-MASS').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) mass
      mass=mass*amu_to_au
    endif
    if (index(string_lec,'-THERMOSTAT').NE.0)then
        call getarg(i+1,string_lec)
        if (string_lec=='SCALE') then
          thermo=1
        elseif (string_lec=='BEREN') then
          thermo=2
          call getarg(i+2,string_lec)
          read(string_lec,*) thermo_val
          thermo_val=thermo_val*fs_to_timeau
        else
          print*, 'Thermostat not well defined'
        stop
      endif
    endif
    if (index(string_lec,'-SCP').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) R_scp
      R_scp=R_scp*ang_to_bohr
      scp=.TRUE.
    endif
    if (index(string_lec,'-PBC').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) box_L
      box_L=box_L*ang_to_bohr
      pbc=.TRUE.
    endif
    if (index(string_lec,'-SAVEFREQ').NE.0)then
      call getarg(i+1,string_lec)
      read(string_lec,*) save_freq
    endif
enddo

!!! Just a check on thermo_val

if (thermo_val.lt.timestep) then
  write(*,*) "the coupling parameter in Berendsen should be larger or equal to timestep"
  stop
endif

!!! Make sure that PBC and SCP are not activated at the same time
if (pbc .and. scp) then
  write(*,*) "PBC and SCP cannot be used together!"
  stop
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Opening of files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

open(10,file='md.xyz')
open(20,file='ENERGIES')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Lecture Initial positions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call lecture_xyz(filename,positions,nbr_atomes) !type xyz

do i=1,nbr_atomes
  do j=1,3
    positions(i,j)=positions(i,j)*ang_to_bohr
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End Lecture Initial positions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate Initial Velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(velocities(nbr_atomes,3))

call initial_velocities(velocities,nbr_atomes,temp_0,mass)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End Generate Initial Velocities
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Generate Initial Forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(forces(nbr_atomes,3))
allocate(distances(nbr_atomes,nbr_atomes))

do i=1,nbr_atomes
  do j=1,3
    forces(i,j)=0.0d0
  enddo
enddo

do i=1,nbr_atomes
  do j=1,nbr_atomes
    distances(i,j)=0.0d0
  enddo
enddo

if (pbc .eqv. .TRUE.) then  !Apply PBC if activated
  !call get_distances_PBC(positions,distances,nbr_atomes,Box_L)
  call get_distances_forces_PBC(positions,distances,forces,nbr_atomes,Box_L,sigma,eps)
else
  call get_distances(positions,distances,nbr_atomes)
  call get_forces(forces,distances,positions,nbr_atomes,sigma,eps)
  if (scp .eqv. .TRUE.) then ! Apply SCP if activated
    call apply_scp(forces,energy_pot,positions,nbr_atomes,R_scp)
  endif
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End Generate Initial Forces
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!! Master MD Loop !!!!!!!!!!!!!!!!! 
do m=1, max_md_Step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Velocity Verlet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Velocities to t+0.5dt
do i=1,nbr_atomes
  do j=1,3
    velocities(i,j)=velocities(i,j)+0.5*timestep*forces(i,j)/mass
  enddo
enddo

do i=1,nbr_atomes
  do j=1,3
    positions(i,j)=positions(i,j)+timestep*velocities(i,j)
  enddo
enddo

!Apply PBC if activated and wrap atoms to unit cell
if (pbc .eqv. .TRUE.) then
  call wrap_atoms_PBC(positions, nbr_atomes, Box_L)  
  call get_distances_forces_PBC(positions,distances,forces,nbr_atomes,Box_L,sigma,eps)
else !call get_forces(forces,distances,positions,nbr_atomes,sigma,eps)
  call get_distances(positions,distances,nbr_atomes)
  call get_forces(forces,distances,positions,nbr_atomes,sigma,eps)
endif
! Apply scp if activated and get pot ener as scp modifies it
call get_pot_ener (distances,nbr_atomes,sigma,eps,energy_pot)
if (scp .eqv. .TRUE.) then 
  call apply_scp(forces,energy_pot,positions,nbr_atomes,R_scp)
endif

do i=1,nbr_atomes
  do j=1,3
    velocities(i,j)=velocities(i,j)+0.5*timestep*forces(i,j)/mass
  enddo
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!End of Velocity Verlet
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Apply Thermostat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_kin_ener (velocities,nbr_atomes,energy_kin,mass)
call get_T (energy_kin,temp,nbr_atomes)
call thermostat(velocities,thermo,temp,temp_0,thermo_val,nbr_atomes,timestep)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Output data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_kin_ener (velocities,nbr_atomes,energy_kin,mass)
call get_T (energy_kin,temp,nbr_atomes)

if ((mod(m,save_freq) .eq. 0) .or. (m .eq. 1)) then
 call output_structure (10,positions,nbr_atomes,energy_pot)
 call output_energy (20,m,timestep,energy_pot,energy_kin,temp)
endif

!!!!!!!!!!!!!!!!! Master MD Loop !!!!!!!!!!!!!!!!! 
enddo

close(10)
close(20)

end program my_first_MD_Code

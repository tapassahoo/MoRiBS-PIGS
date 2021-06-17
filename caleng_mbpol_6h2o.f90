  subroutine mbxeng(com,Eulang,Vpot)
  implicit none

  double precision :: angHOH,dOH,pi,xH,zH,ang1
  double precision :: Vpot
  double precision, allocatable :: coord(:)
  double precision, allocatable :: box(:)
  double precision, allocatable :: com(:)
  double precision, allocatable :: Eulang(:)
  double precision, allocatable :: ROwf(:)
  double precision, allocatable :: RH1wf(:)
  double precision, allocatable :: RH2wf(:)
  double precision, allocatable :: com_1(:)
  double precision, allocatable :: Eulang_1(:)
  double precision, allocatable :: RO_1_sf(:)
  double precision, allocatable :: RH1_1_sf(:)
  double precision, allocatable :: RH2_1_sf(:)
  double precision, DIMENSION(3,3) :: rotmat_1

  character(len=5), allocatable :: at_name(:)
  character(len=5), allocatable :: monomers(:)
  character(len=20) :: xyz, json_file
  
  integer, allocatable :: nats(:)
  character :: atom

  integer :: n_at,i,j,nmon,k,ii,jj,j1

  external initialize_system
  external get_energy
  external finalize_system

  ! need to define here the numebr of monomers, and number of atoms in each one
  ! In the case of a water cluster of 256 molecules, the variables should be:
  
  ! nmon = 256
  ! allocate(nats(nmon),monomers(nmon))

  ! do i=1,nmon
  !   nats(i) = 3
  !   monomers(i) = "h2o"
  ! enddo

  ! This will work for the test case, which contains 2 co2 monomers
  ! at the beggining, and then 2 water molecules
  nmon = 6
  allocate(nats(nmon),monomers(nmon),box(9))

  ! h2o
  do i=1,nmon
    nats(i) = 3
    monomers(i) = "h2o"
  enddo

  ! Set json file
  json_file = 'mbx_gas_phase.json'//CHAR(0)
  ! Open file called input.xyz. 
  ! Contains the coordinates in XYZ format.
  ! Assumes halide at the end (because the C++ code does that)
  !xyz='input.xyz'
  !open(unit = 55, file = xyz, status = 'old', action = 'read')

  ! Number of atoms
  !read(55,*) n_at

  n_at=18
  ! Allocation of the coord vector and the atom name vector,
  ! as so as the vector that will store the gradients
  allocate(coord(n_at*3),at_name(n_at))

  ! Reading comment file from XYZ. Not saving it
  !read(55,*)

  ! Reading and saving the names and coordinates
  ! Coordinates will be XYZXYZXYZ... OHHOHH...OHHX
  !do i=1,n_at
  !  read(55,*) at_name(i), (coord(3*(i-1)+j),j=1,3)
  !enddo

  ! Done reading. Closing file.
  !close(55)

  at_name(1)='O'
  at_name(2)='H'
  at_name(3)='H'
  at_name(4)='O'
  at_name(5)='H'
  at_name(6)='H'
  at_name(7)='O'
  at_name(8)='H'
  at_name(9)='H'
  at_name(10)='O'
  at_name(11)='H'
  at_name(12)='H'
  at_name(13)='O'
  at_name(14)='H'
  at_name(15)='H'
  at_name(16)='O'
  at_name(17)='H'
  at_name(18)='H'

  allocate(ROwf(3),RH1wf(3),RH2wf(3))

  angHOH=107.4d0
  dOH=0.9419d0
  pi=4.0d0*datan(1.0d0)
  ang1=(angHOH*pi)/180.0d0
  zH=ROwf(3)-dsqrt(0.5*dOH*dOH*(1.0+dcos(ang1))) 
  xH=dsqrt(dOH*dOH-(ROwf(3)-zH)*(ROwf(3)-zH))
  RH1wf(1)=xH
  RH1wf(2)=0.0d0
  RH1wf(3)=zH
  RH2wf(1)=-RH1wf(1)
  RH2wf(2)=RH1wf(2)
  RH2wf(3)=RH1wf(3)

  do i=1,nmon
    do j=1,3
        ii=j+(i-1)*3
        com_1(j)=com(ii)
        Eulang_1(j)=Eulang(ii)
    enddo

    call matpre(Eulang_1, rotmat_1)
    do j=1,3
        RO_1_sf(j)=0.d0
    enddo
    call rottrn(rotmat_1, ROwf, RO_1_sf, com_1)
    do j=1,3
        jj=j+(i-1)*9
        coord(jj)=RO_1_sf(j)
    enddo

    do j=1,3
        RH1_1_sf(j)=0.d0
    enddo
    call rottrn(rotmat_1, RH1wf, RH1_1_sf, com_1)
    do j=4,6
        jj=j+(i-1)*9
        j1=j-3
        coord(jj)=RH1_1_sf(j1)
    enddo

    do j=1,3
        RH2_1_sf(j)=0.d0
    enddo
    call rottrn(rotmat_1, RH2wf, RH2_1_sf, com_1)
    do j=7,9
        jj=j+(i-1)*9
        j1=j-6
        coord(jj)=RH2_1_sf(j1)
    enddo
  enddo

  ! MAKE SURE YOU HAVE LEN=5 IN THE CHAR ARRAY!!!!!!

  do i =1,n_at
    at_name(i)=trim(at_name(i))//CHAR(0)
  enddo
  do i=1,nmon
    monomers(i) = trim(monomers(i))//CHAR(0)
  enddo
    
  ! initialize_system
  ! @coord is a double 1D array with the coordinates of all atoms.
  !        It needs xyzxyzxyzxyz... Basically reading an xyz file
  !        and has length 3*n_at
  ! @nats is an integer 1D array of length nmon with the number of atoms
  !       in each monomer.
  ! @at_name is a char 2D array. The second dimention needs to be 5. 
  !          The first dimention needs to be of length of number of atoms
  !          It contains the name of the atoms as they appear in the 
  !          periodic table
  ! @monomers is a char 2D array. The second dimention needs to be 5.
  !           The first dimention needs to be of length of number of mons
  !           It contains the names of the monomers (h2o for water,
  !           li, na, k, rb, cs, f, cl, br, i for the alkali metal ions
  !           and the halides)
  ! @nmon is an integer with the total number of monomers 

  ! get_energy
  ! @coord is a double 1D array with the coordinates of all atoms.
  !        It needs xyzxyzxyzxyz... Basically reading an xyz file
  !        and has length 3*n_at
  ! @n_at is an integer with the number of atoms
  ! @Vpot is an output double variable that will contain the energy 
  !       after the call to the energy function

  ! get_energy_g 
  ! @grads is an output double 1D array that will contain the gradients 
  !        (XYZ) of each atom as gx1 gy1 gz1 gx2 gy2 gz2 ...
 
  ! First, initialize system
  call initialize_system(coord, nats, at_name, monomers, nmon, json_file)

  ! Energy call no gradients no pbc
  call get_energy(coord, n_at, Vpot)
  !write(*,*) "Energy = " , Vpot

  ! Don't forget to free the system
  call finalize_system()


END SUBROUTINE


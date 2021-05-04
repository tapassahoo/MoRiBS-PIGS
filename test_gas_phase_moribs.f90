program test
  implicit none

  double precision :: Vpot
  double precision, allocatable :: coord(:)
  double precision, allocatable :: box(:)

  character(len=5), allocatable :: at_name(:)
  character(len=5), allocatable :: monomers(:)
  character(len=20) :: xyz, json_file
  
  integer, allocatable :: nats(:)
  character :: atom

  integer :: n_at,i,j,nmon,k

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

  coord(1)=0.803889
  coord(2)=0.381762 
  coord(3)= -1.685143
  coord(4)=0.362572
  coord(5)=-0.448201 
  coord(6)=-1.556674
  coord(7)=1.668734
  coord(8)=0.275528
  coord(9)=-1.301550
  coord(10)=0.666169
  coord(11)=-0.420958
  coord(12)=1.707749
  coord(13)=0.236843
  coord(14)=0.404385
  coord(15)=1.523931
  coord(16)=0.226003
  coord(17)=-1.053183
  coord(18)=1.153395
  coord(19)=2.996112
  coord(20)=0.001740
  coord(21)=0.125207
  coord(22)=2.356345
  coord(23)=-0.159970
  coord(24)=0.813642
  coord(25)=3.662033
  coord(26)=-0.660038
  coord(27)=0.206711
  coord(28)=-0.847903
  coord(29)=-1.777751
  coord(30)=-0.469278
  coord(31)=-1.654759
  coord(32)=-1.281222
  coord(33)=-0.344427
  coord(34)=-1.091666
  coord(35)=-2.653858
  coord(36)=-0.718356
  coord(37)=-2.898828
  coord(38)=0.065636
  coord(39)=0.089967
  coord(40)=-3.306527
  coord(41)=0.037245
  coord(42)=0.940083
  coord(43)=-2.312757
  coord(44)=0.817025
  coord(45)=0.097526
  coord(46)=-0.655160
  coord(47)=1.814997
  coord(48)=0.176741
  coord(49)=-0.134384
  coord(50)=1.449649
  coord(51)=-0.543456
  coord(52)=-0.526672
  coord(53)=2.749233
  coord(54)=0.167243

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
  write(*,*) "Energy = " , Vpot

  ! Don't forget to free the system
  call finalize_system()


end program


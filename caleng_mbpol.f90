!program test
  !implicit none

  !double precision :: Vpot1
  !double precision, allocatable :: com(:)
  !double precision, allocatable :: Eulang(:)
  !integer :: nmon, ndim, nats

  !nmon=2
  !ndim=3
  !nats=3

  !allocate(com(nmon*ndim),Eulang(nmon*ndim))

  !com(1)=0.0
  !com(2)=0.0
  !com(3)=0.0
  !com(4)=0.0
  !com(5)=0.0
  !com(6)=3.0

  !Eulang(1)=0.0
  !Eulang(2)=1.0
  !Eulang(3)=0.0
  !Eulang(4)=0.0
  !Eulang(5)=1.0
  !Eulang(6)=0.0

  !call mbxeng(com,Eulang,Vpot1,nmon,ndim)
  !write(6,*) Vpot1
!end program test
!
  subroutine mbxeng(com,Eulang,Vpot1)
  implicit none

  integer :: nmon,ndim
  double precision, dimension(6) :: com
  double precision, dimension(6) :: Eulang
  double precision :: Vpot1

  double precision, dimension(3) :: ROwf
  double precision, dimension(3) :: RH1wf
  double precision, dimension(3) :: RH2wf
  double precision, dimension(3) :: RO_sf
  double precision, dimension(3) :: RH1_sf
  double precision, dimension(3) :: RH2_sf
  double precision, dimension(3) :: Eulang1
  double precision, dimension(3) :: com1
  double precision, DIMENSION(3,3) :: rotmat

  double precision :: Vpot

  double precision, dimension(18) :: coord
  character(len=5), dimension(6) :: at_name
  character(len=5), dimension(2) :: monomers
  double precision, dimension(9) :: box
  character(len=20) :: xyz, json_file
  
  integer, dimension(2) :: nats
  character :: atom

  integer :: n_at,i,j,k,ii,jj,j1

  external initialize_system
  external get_energy
  external finalize_system

  nmon=2
  ndim=3
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
  !nmon = 2

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

  n_at=nmon*3
  ! Allocation of the coord vector and the atom name vector,
  ! as so as the vector that will store the gradients

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


  ROwf(1)=0.0
  ROwf(2)=0.0
  ROwf(3)=0.062404
  RH1wf(1)=-0.759103
  RH1wf(2)=0.0
  RH1wf(3)=-0.495212
  RH2wf(1)=-RH1wf(1)
  RH2wf(2)=RH1wf(2)
  RH2wf(3)=RH1wf(3)

  do i=1,nmon
    do j=1,ndim
        ii=j+(i-1)*ndim
        com1(j)=com(ii)
        Eulang1(j)=Eulang(ii)
    enddo

    call matpre1(Eulang1, rotmat)
    do j=1,3
        RO_sf(j)=0.d0
    enddo
    call rottrn1(rotmat, ROwf, RO_sf, com1)
    do j=1,3
        jj=j+(i-1)*9
        coord(jj)=RO_sf(j)
    enddo

   do j=1,3
       RH1_sf(j)=0.d0
    enddo
    call rottrn1(rotmat, RH1wf, RH1_sf, com1)
    do j=4,6
        jj=j+(i-1)*9
        j1=j-3
        coord(jj)=RH1_sf(j1)
    enddo

    do j=1,3
        RH2_sf(j)=0.d0
    enddo
    call rottrn1(rotmat, RH2wf, RH2_sf, com1)
    do j=7,9
        jj=j+(i-1)*9
        j1=j-6
        coord(jj)=RH2_sf(j1)
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
  Vpot1=Vpot

END SUBROUTINE mbxeng

    subroutine matpre1(Eulang,rotmat)
    implicit none

    double precision, dimension(3), intent(in)::Eulang
    double precision, dimension(3,3), intent(out)::rotmat

    double precision :: phi, theta, chi, cp, sp, ct, st, ck, sk

    phi=Eulang(1)
    theta=Eulang(2)
    chi=Eulang(3)

    cp=cos(phi)
    sp=sin(phi)
    ct=cos(theta)
    st=sin(theta)
    ck=cos(chi)
    sk=sin(chi)

    rotmat(1,1)=cp*ct*ck-sp*sk
    rotmat(1,2)=-cp*ct*sk-sp*ck
    rotmat(1,3)=cp*st
    rotmat(2,1)=sp*ct*ck+cp*sk
    rotmat(2,2)=-sp*ct*sk+cp*ck
    rotmat(2,3)=sp*st
    rotmat(3,1)=-st*ck
    rotmat(3,2)=st*sk
    rotmat(3,3)=ct

    end subroutine matpre1

    subroutine rottrn1(rotmat,rwf,rsf,rcom)
    implicit none
    double precision, dimension(3), intent(in) :: rwf,rcom
    double precision, dimension(3,3), intent(in) :: rotmat(3,3)
    double precision, dimension(3), intent(out) :: rsf
    integer :: i,j

    do i=1,3
      rsf(i)=rcom(i)
      do j=1,3
        rsf(i)=rsf(i)+rotmat(i,j)*rwf(j)
      enddo
    enddo
   end subroutine rottrn1

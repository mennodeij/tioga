


program grid2d

  use gridtype
  implicit none
  include "mpif.h"

  integer :: ier, myproc, numprocs
  integer :: i, j, k, o, v, ng
  
  type dblArray
    real*8, allocatable :: data(:)
  end type
  
  type intArray
    integer*4, allocatable :: data(:)
  end type
  
  real*8 :: xyz(3), x, y, z
  integer :: cell(8)
    
  type(dblArray) :: coord(2), nodeResolution(2), cellResolution(2)
  type(intArray) :: connectivity(2), iblank(2), obcnode(2), wbcnode(2), iblank_c(2)
  
  integer :: nvertx(2), nverty(2), nvertz(2)
  integer :: ncellx(2), ncelly(2), ncellz(2)
  integer :: nverts(2), ncells(2), nwbc(2), nobc(2)
  
    
  real*8 :: sizex(2), sizey(2), sizez(2)
  real*8 :: stepx(2), stepy(2), stepz(2)
  
  ! scale 
  REAL*8 :: scaling(3,3) =  RESHAPE( (/ 0.5d0, 0.0d0, 0.0d0,  &
                                        0.0d0, 0.5d0, 0.0d0,  &
                                        0.0d0, 0.0d0, 1.0d0 /), &
                                      (/ 3, 3 /) &
                                    )
  ! rotate over 15.93 degrees
  REAL*8 :: rotation(3,3) = RESHAPE( (/ 0.96159773301087720159134301585929d0, -0.27446274768780867507596949035934d0, 0.0d0,  &
                                        0.27446274768780867507596949035934d0,  0.96159773301087720159134301585929d0, 0.0d0,  &
                                        0.0d0, 0.0d0, 1.0d0 /), &
                                      (/ 3, 3 /) &
                                    )
  REAL*8 :: translation(3) 
  
  
  call mpi_init(ier)
  call mpi_comm_rank(mpi_comm_world, myproc, ier)
  call mpi_comm_size(mpi_comm_world, numprocs, ier)
  
  call tioga_init_f90(mpi_comm_world)
  call mpi_barrier(mpi_comm_world,ier)
  
  if (numprocs /= 1) then
    print *, "Only one process supported."
    call exit(-1)
  endif
  
  ncellx(1) = 32
  ncelly(1) = 16
  ncellz(1) = 1
  
  ncellx(2) = 64
  ncelly(2) = 32
  ncellz(2) = 1
  
  sizex(1) = 2.0d0 
  sizey(1) = 1.0d0
  sizez(1) = 0.1d0
  
  sizex(2) = 2.0d0 
  sizey(2) = 1.0d0
  sizez(2) = 0.1d0
  
  do ng=1,2
    stepx(ng) = sizex(ng)/ncellx(ng)
    stepy(ng) = sizey(ng)/ncelly(ng)
    stepz(ng) = sizez(ng)/ncellz(ng)
  end do
  
  nvertx = 1+ncellx
  nverty = 1+ncelly
  nvertz = 1+ncellz
  
  do ng=1,2
    nverts(ng) = nvertx(ng)*nverty(ng)*nvertz(ng)
    ncells(ng) = ncellx(ng)*ncelly(ng)*ncellz(ng)
    
    allocate(coord(ng)%data(3*nverts(ng)), stat=ier)
    allocate(iblank(ng)%data(nverts(ng)), stat=ier)
    allocate(connectivity(ng)%data(8*ncells(ng)), stat=ier)
    allocate(iblank_c(ng)%data(ncells(ng)), stat=ier)
    ! node and cell resolution is used to indicate which nodes/cells solve
    ! the PDE's, see https://github.com/jsitaraman/tioga/issues/44 for more info
    allocate(nodeResolution(ng)%data(nverts(ng)), stat=ier)
    allocate(cellResolution(ng)%data(ncells(ng)), stat=ier)
    
    if (ng == 1) then 
      nwbc(ng) = 0
      nobc(ng) = nvertz(ng)*(2*(ncellx(ng)+ncelly(ng)))
      ! grid 1 is the foreground grid: use smaller resolution values
      nodeResolution(ng)%data = 1.0d0
      cellResolution(ng)%data = 1.0d0
    else 
      nwbc(ng) = 0
      nobc(ng) = 0
      ! grid 2 is the background grid: use larger resolution values
      nodeResolution(ng)%data = 10.0d0
      cellResolution(ng)%data = 10.0d0
    endif
    
    allocate(wbcnode(ng)%data(nwbc(ng)))
    allocate(obcnode(ng)%data(nobc(ng)))

    o = 0
    v = 0
    z = 0.0d0
    y = 0.0d0
    x = 0.0d0
    
    do k=1,nvertz(ng)
      do j=1,nverty(ng)
        do i=1,nvertx(ng)
          v = v+1
          coord(ng)%data(3*v-2) = x
          coord(ng)%data(3*v-1) = y
          coord(ng)%data(3*v-0) = z
          x = x + stepx(ng)
          
          ! tag the overset boundary nodes
          if (ng .eq. 1) then
            if (i .eq. 1 .or. i .eq. nvertx(ng)) then
              o = o + 1
              obcnode(ng)%data(o) = v
            else if (j .eq. 1 .or. j .eq. nverty(ng)) then
              o = o + 1
              obcnode(ng)%data(o) = v
            endif
          endif
        end do
        y = y + stepy(ng)
        x = 0.0d0
      end do
      z = z + stepz(ng)
      y = 0.0d0
    end do

    
    v = 0
    do k=1,ncellz(ng)
      do j=1,ncelly(ng)
        do i=1,ncellx(ng)
          cell(1) = i+0+(j-1)*nvertx(ng)
          cell(2) = i+1+(j-1)*nvertx(ng)
          cell(3) = i+1+(j-0)*nvertx(ng)
          cell(4) = i+0+(j-0)*nvertx(ng)
          
          cell(5) = cell(1)+nvertx(ng)*nverty(ng)
          cell(6) = cell(2)+nvertx(ng)*nverty(ng)
          cell(7) = cell(3)+nvertx(ng)*nverty(ng)
          cell(8) = cell(4)+nvertx(ng)*nverty(ng)
          v = v+1
          connectivity(ng)%data(8*v-7:8*v) = cell(1:8)
        end do
      end do
    end do
  end do
  
  ! transform the foreground grid
  ng = 1
  translation(1) = sizex(1)/4.0d0+0.03563
  translation(2) = sizey(1)/4.0d0+0.02987
  translation(3) = 0.0d0
  do v=1,nverts(1)
    xyz(1:3) = coord(1)%data(3*v-2:3*v) 
    xyz = MATMUL(scaling, xyz)
    coord(ng)%data(3*v-2:3*v) = MATMUL(rotation, xyz) + translation
  end do
  
  do ng=1,2
    call tioga_registergrid_data_mb(ng, ng, &
                                    nverts(ng), &
                                    coord(ng)%data, &
                                    iblank(ng)%data, &
                                    nwbc(ng), &
                                    nobc(ng), &
                                    wbcnode(ng)%data, &
                                    obcnode(ng)%data, &
                                    1, &
                                    8, &
                                    ncells(ng), &
                                    connectivity(ng)%data)
    call tioga_setcelliblank_multi(ng, iblank_c(ng)%data)
    call tioga_setresolutions_multi(ng, nodeResolution(ng)%data, cellResolution(ng)%data)
  end do
  
  call tioga_setnfringe(2) ! value of 2 (or 3) does not seem to work; always 1 layerof fringe cells
  call tioga_preprocess_grids
  call tioga_performconnectivity
  
  ! check that all outer cells are fringe on the foreground grid
  
  v = 0
  do k=1,ncellz(1)
    do j=1,ncelly(1)
      do i=1,ncellx(1)
        v = v + 1
        if (i .eq. 1 .or. i .eq. ncellx(1)) then
          if (iblank_c(1)%data(v) /= -1) print*, "cell ", v, " not fringe"
        elseif (j .eq. 1 .or. j .eq. ncelly(1)) then
          if (iblank_c(1)%data(v) /= -1) print*, "cell ", v, " not fringe"
        else
          if (iblank_c(1)%data(v) /= +1) print*, "cell ", v, " not field"
        endif
      end do
    end do
  end do

  call tioga_writeoutputfiles(0,"row")

  call tioga_delete
  call mpi_finalize(ier)

end program
          
          
          
          
          
          
          
          
          
          
          
          
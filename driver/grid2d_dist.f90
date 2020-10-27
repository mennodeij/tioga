


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
  
  type boolArray
    logical, allocatable :: data(:)
  end type
  
  real*8 :: xyz(3), x, y, z
  integer :: cell(8)
    
  type(dblArray) :: coord(2), nodeResolution(2), cellResolution(2)
  type(intArray) :: connectivity(2), iblank(2), obcnode(2), wbcnode(2), iblank_c(2)
  type(intArray):: vertOnProc(2), cellOnProc(2), iblankOnProc(2), connectivityOnProc(2)
  type(dblArray):: coordOnProc(2), cellResolutionOnProc(2), nodeResolutionOnProc(2)
  type(intArray):: obcOnProc(2), wbcOnProc(2)
  
  type(boolArray):: isHoleNode(2), isFringeNode(2), isInNode(2)
  type(boolArray):: isWhatNode
  
  integer :: nCellsPerProc(2), nVertsOnProc(2), nwbcOnProc(2), nobcOnProc(2), totalCells
  integer :: nFringeVerts,nTotalFringeVerts,nHoleVerts,nTotalHoleVerts,nInVerts,nTotalInVerts
  
  integer :: nvertx(2), nverty(2), nvertz(2)
  integer :: ncellx(2), ncelly(2), ncellz(2)
  integer :: nverts(2), ncells(2), nwbc(2), nobc(2)

  real*8 :: sizex(2), sizey(2), sizez(2)
  real*8 :: stepx(2), stepy(2), stepz(2)
  
  logical :: grouped

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
  
  ! if (numprocs /= 1) then
  !   print *, "Only one process supported."
  !   call exit(-1)
  ! endif
  
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

  ! all process now have the whole grid - distribute between processes 
  ! grouped: in batches 
  ! not grouped: in 'pathological way' => each nth cell goes to the nth process.
  
  grouped = .FALSE.
  nCellsPerProc = 0
  do ng=1,2
    ! determine which cells go where
    if (grouped) then
      ! each process gets a contiguous list of cells
      nCellsPerProc(ng) = ncells(ng) / numprocs
      if (MOD(ncells(ng), numprocs) /= 0 .AND. myproc == numprocs-1) then
        nCellsPerProc(ng) = nCellsPerProc(ng) + MOD(ncells(ng), numprocs)
      endif
      
      ALLOCATE(vertOnProc(ng)%data(nverts(ng)))
      ALLOCATE(cellOnProc(ng)%data(ncells(ng)))

      vertOnProc(ng)%data= -1
      cellOnProc(ng)%data= -1
      
      j = 1+(myproc)*(ncells(ng) / numprocs)
      k = j + nCellsPerProc(ng) -1

      do i=j,k
        cellOnProc(ng)%data(i) = 1
        do j=1,8
          v = connectivity(ng)%data((i-1)*8+j)
          vertOnProc(ng)%data(v) = 1
        end do
       end do
      
    else
      ! each process gets the n'th cell
      do i=1,ncells(ng)
        if (MOD(i, numprocs) == myproc) nCellsPerProc(ng) = nCellsPerProc(ng) + 1
      end do

      ALLOCATE(vertOnProc(ng)%data(nverts(ng)))
      ALLOCATE(cellOnProc(ng)%data(ncells(ng)))

      vertOnProc(ng)%data= -1
      cellOnProc(ng)%data= -1

      do i=1,ncells(ng)
        if (MOD(i, numprocs) == myproc) then 
          cellOnProc(ng)%data(i) = 1
          do j=1,8
            v = connectivity(ng)%data((i-1)*8+j)
            vertOnProc(ng)%data(v) = 1
          end do
        endif
      end do
    endif

    nVertsOnProc(ng) = COUNT(vertOnProc(ng)%data == 1)

    ALLOCATE(coordOnProc(ng)%data(3*nVertsOnProc(ng)))
    ALLOCATE(iblankOnProc(ng)%data(nVertsOnProc(ng)))
    ALLOCATE(connectivityOnProc(ng)%data(8*nCellsPerProc(ng)))
    
    coordOnProc(ng)%data = -100
    j = 0
    do i=1,nverts(ng)
      if (vertOnProc(ng)%data(i) == 1) then ! 1 => present on proc
        j = j+1
        coordOnProc(ng)%data(3*j-2:3*j) = coord(ng)%data(3*i-2:3*i)
        vertOnProc(ng)%data(i) = j ! keep for lookup
      endif
    end do
    
    j = 0
    do i=1,ncells(ng)
      if (cellOnProc(ng)%data(i) == 1) then ! 1 => present on proc
        do k=1,8
          v = connectivity(ng)%data(8*(i-1)+k)
          connectivityOnProc(ng)%data(8*j+k) = vertOnProc(ng)%data(v)
        end do
        j = j+1
      endif
    end do
    
    nwbcOnProc(ng) = 0
    do i=1,nwbc(ng)
      if (vertOnProc(ng)%data(wbcnode(ng)%data(i)) > 0) then
        nwbcOnProc(ng) = nwbcOnProc(ng) + 1
      endif
    end do
    
    nobcOnProc(ng) = 0
    do i=1,nobc(ng)
      if(vertOnProc(ng)%data(obcnode(ng)%data(i)) > 0) then
        nobcOnProc(ng) = nobcOnProc(ng) + 1
      endif
    end do
    
    ALLOCATE(obcOnProc(ng)%data(nobcOnProc(ng)))
    ALLOCATE(wbcOnProc(ng)%data(nwbcOnProc(ng)))
    
    j = 0
    do i=1,nwbc(ng)
      if (vertOnProc(ng)%data(wbcnode(ng)%data(i)) > 0) then
        j = j+1
        wbcOnProc(ng)%data(j) = vertOnProc(ng)%data(wbcnode(ng)%data(i))
      endif
    end do
    
    j = 0
    do i=1,nobc(ng)
      if(vertOnProc(ng)%data(obcnode(ng)%data(i)) > 0) then
        j = j+1
        obcOnProc(ng)%data(j) = vertOnProc(ng)%data(obcnode(ng)%data(i))
      endif
    end do
    
    call tioga_registergrid_data_mb(ng, ng, &
                                    nVertsOnProc(ng), &
                                    coordOnProc(ng)%data, &
                                    iblankOnProc(ng)%data, &
                                    nwbcOnProc(ng),  &
                                    nobcOnProc(ng), &
                                    wbcOnProc(ng)%data, &
                                    obcOnProc(ng)%data, &
                                    1, &
                                    8, &
                                    nCellsPerProc(ng), &
                                    connectivityOnProc(ng)%data)
    ALLOCATE(nodeResolutionOnProc(ng)%data(nVertsOnProc(ng)))
    ALLOCATE(cellResolutionOnProc(ng)%data(nCellsPerProc(ng)))

    if (nobc(ng) > 0) then
      nodeResolutionOnProc(ng)%data = 1.0d0
      cellResolutionOnProc(ng)%data = 1.0d0
    else
      nodeResolutionOnProc(ng)%data = 10.0d0
      cellResolutionOnProc(ng)%data = 10.0d0
    endif

    call tioga_setresolutions_multi(ng, nodeResolutionOnProc(ng)%data, cellResolutionOnProc(ng)%data)
  end do ! ng 

  ! do ng=1,2
  !   call tioga_registergrid_data_mb(ng, ng, &
  !                                   nverts(ng), &
  !                                   coord(ng)%data, &
  !                                   iblank(ng)%data, &
  !                                   nwbc(ng), &
  !                                   nobc(ng), &
  !                                   wbcnode(ng)%data, &
  !                                   obcnode(ng)%data, &
  !                                   1, &
  !                                   8, &
  !                                   ncells(ng), &
  !                                   connectivity(ng)%data)
  !   call tioga_setcelliblank_multi(ng, iblank_c(ng)%data)
  !   call tioga_setresolutions_multi(ng, nodeResolution(ng)%data, cellResolution(ng)%data)
  ! end do
  
  call tioga_setnfringe(1) ! only works for the background grid
  call tioga_setmexclude(3)
  call tioga_preprocess_grids
  
  !call mpi_finalize(ier)
  !call exit(0)
  
  call tioga_performconnectivity
  call tioga_reduce_fringes
  
  nFringeVerts = 0
  nHoleVerts = 0
  nInVerts = 0
  
  
  do ng=1,2
    ALLOCATE(isFringeNode(ng)%data(nverts(ng)))
    ALLOCATE(isHoleNode(ng)%data(nverts(ng)))
    ALLOCATE(isInNode(ng)%data(nverts(ng)))
    ALLOCATE(isWhatNode%data(nverts(ng)))
    isFringeNode(ng)%data = .FALSE.
    isHoleNode(ng)%data = .FALSE.
    isInNode(ng)%data = .FALSE.
    
    
    do i=1,nverts(ng)
      if (vertOnProc(ng)%data(i) > 0) then
        if (iblankOnProc(ng)%data(vertOnProc(ng)%data(i)) == -1) isFringeNode(ng)%data(i) = .TRUE.
        if (iblankOnProc(ng)%data(vertOnProc(ng)%data(i)) ==  0) isHoleNode(ng)%data(i) = .TRUE.
        if (iblankOnProc(ng)%data(vertOnProc(ng)%data(i)) == +1) isInNode(ng)%data(i) = .TRUE. 
      endif
    enddo

    isWhatNode%data = .FALSE.
    CALL MPI_REDUCE(isFringeNode(ng)%data, isWhatNode%data, nverts(ng), MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, ier)
    nTotalFringeVerts = COUNT(isWhatNode%data)

    isWhatNode%data = .FALSE.
    CALL MPI_REDUCE(isHoleNode(ng)%data, isWhatNode%data, nverts(ng), MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, ier)
    nTotalHoleVerts = COUNT(isWhatNode%data)

    isWhatNode%data = .FALSE.
    CALL MPI_REDUCE(isInNode(ng)%data, isWhatNode%data, nverts(ng), MPI_LOGICAL, MPI_LOR, 0, MPI_COMM_WORLD, ier)
    nTotalInVerts = COUNT(isWhatNode%data)

    if (myproc == 0) print *, "Grid ", ng, " has ", nTotalInVerts, " IN, ", nTotalFringeVerts, " FRINGE and ", nTotalHoleVerts, " HOLE nodes."
    DEALLOCATE(isWhatNode%data)
  end do

  call tioga_writeoutputfiles(0,"row")

  call tioga_delete
  call mpi_finalize(ier)
  
  

  do ng=1,2
    deallocate(coord(ng)%data)
    deallocate(iblank(ng)%data)
    deallocate(connectivity(ng)%data)
    deallocate(iblank_c(ng)%data)
    deallocate(nodeResolution(ng)%data)
    deallocate(cellResolution(ng)%data)
    deallocate(wbcnode(ng)%data)
    deallocate(obcnode(ng)%data)
    
    deallocate(vertOnProc(ng)%data)
    deallocate(cellOnProc(ng)%data)
    deallocate(coordOnProc(ng)%data)
    deallocate(iblankOnProc(ng)%data)
    deallocate(connectivityOnProc(ng)%data)
    deallocate(wbcOnProc(ng)%data)
    deallocate(obcOnProc(ng)%data)
    deallocate(cellResolutionOnProc(ng)%data)
    deallocate(nodeResolutionOnProc(ng)%data)
    
    deallocate(isFringeNode(ng)%data)
    deallocate(isHoleNode(ng)%data)
    deallocate(isInNode(ng)%data)
  end do

end program


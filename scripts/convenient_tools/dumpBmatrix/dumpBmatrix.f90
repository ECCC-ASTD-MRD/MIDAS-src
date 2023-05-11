
program dumpBmatrix
  implicit none
  integer               :: extractDate, nlev_T, nlev_M,Vcode, ip1_sfc, ip1_T_2m
  integer               :: ip1_M_10m,varCount,nkgdim, nBands
  integer               :: nulmat, ierr
  integer               :: indexBand, indexI,indexJ
  integer               :: levelIndex, controlVectorIndex, varIndex
  integer, allocatable  :: ip1_T(:), ip1_M(:), controlVectorIp1(:)
  integer, external     :: fnom, fclos
  real(4)               :: latitude, longitude
  real(8), allocatable  :: Bmatrix(:,:)

  character(len=4), allocatable :: varList(:), controlVectorVar(:)
  character(len=256)            :: bmatFileName

  
  if ( command_argument_count() /= 1 ) then
    write(*,*) '<!> wrong argument number <!>'
    write(*,*) '    SYNOPSIS:'
    write(*,*) '       ./dumpBmatrix <bmatFileName>'
    call exit(1)
  end if

  call get_command_argument(1,bmatFileName)

  nulmat = 0
  ierr = fnom(nulmat,trim(bmatFileName),'FTN+SEQ+UNF',0)
  read(nulmat) extractDate, nlev_T, nlev_M, Vcode, &
      ip1_sfc, ip1_T_2m, ip1_M_10m, varCount, nkgdim, nBands
  write(*,*) extractDate
  write(*,*) nlev_T, nlev_M, Vcode
  write(*,*) ip1_sfc, ip1_T_2m, ip1_M_10m
  write(*,*) varCount, nkgdim, nBands
  allocate(ip1_T(nlev_T),ip1_M(nlev_M),varList(varCount))
  allocate(Bmatrix(nkgdim,nkgdim))
  read(nulmat) ip1_T(:), ip1_M(:), varList(:)
  write(*,*) ip1_T
  write(*,*) ip1_M
  write(*,*) varList


  allocate(controlVectorVar(nkgdim))
  allocate(controlVectorIp1(nkgdim))
  controlVectorIndex = 0
  do varIndex = 1, varCount
    select case(trim(varList(varIndex)))
    case('TT','HU')
      do levelIndex= 1, nlev_T
        controlVectorIndex =  controlVectorIndex + 1
        controlVectorVar(controlVectorIndex) = varList(varIndex)
        controlVectorIp1(controlVectorIndex) = ip1_T(levelIndex)
      end do
    case('UU','VV')
      do levelIndex= 1, nlev_M
        controlVectorIndex =  controlVectorIndex + 1
        controlVectorVar(controlVectorIndex) = varList(varIndex)
        controlVectorIp1(controlVectorIndex) = ip1_M(levelIndex)
      end do
    case('P0','TG')
      controlVectorIndex =  controlVectorIndex + 1
      controlVectorVar(controlVectorIndex) = varList(varIndex)
      controlVectorIp1(controlVectorIndex) = ip1_sfc
    case default
      write(*,*) "Unknown variable ", varList(varIndex)
      call exit(1)
    end select
  end do
  
  do indexBand = 1, nBands 
    read(nulmat) latitude, longitude,  Bmatrix(:,:)
    write(*,*) '--- band ', indexBand, latitude, longitude
    do indexI = 1, nkgdim
      do indexJ =indexI, nkgdim
        write(*,'(3(i3,1x),2(A4,1x),2(i12,1x),e23.15)') indexBand, indexI,indexJ, &
            controlVectorVar(indexI),controlVectorVar(indexJ), &
            controlVectorIp1(indexI),controlVectorIp1(indexJ), &
            Bmatrix(indexI,indexJ)
      end do
    end do
  end do

  ierr=fclos(nulmat)
  deallocate(ip1_T,ip1_M,varList)
  deallocate(Bmatrix)
  deallocate(controlVectorVar)
  deallocate(controlVectorIp1)
  
end program dumpBmatrix

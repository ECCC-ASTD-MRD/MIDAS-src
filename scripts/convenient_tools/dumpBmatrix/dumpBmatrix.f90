!--------------------------------------- LICENCE BEGIN -----------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------


program dumpBmatrix
  implicit none
  integer               :: extractDate, nlev_T, nlev_M,Vcode, ip1_sfc, ip1_T_2m
  integer               :: ip1_M_10m,varCount,nkgdim, nBands
  integer               :: nulmat, ierr
  integer               :: indexBand, indexI,indexJ
  integer, allocatable  :: ip1_T(:), ip1_M(:)
  integer, external     :: fnom, fclos
  real(4)               :: latitude, longitude
  real(8), allocatable  :: Bmatrix(:,:)

  character(len=4), allocatable :: varList(:)
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
  
  do indexBand = 1, nBands 
    read(nulmat) latitude, longitude,  Bmatrix(:,:)
    write(*,*) '--- band ', indexBand, latitude, longitude
    do indexI = 1, nkgdim
      do indexJ =indexI, nkgdim
        write(*,*) indexBand, indexI,indexJ,Bmatrix(indexI,indexJ)
      end do
    end do
  end do

  ierr=fclos(nulmat)
  deallocate(ip1_T,ip1_M,varList)
  deallocate(Bmatrix)
  
end program dumpBmatrix

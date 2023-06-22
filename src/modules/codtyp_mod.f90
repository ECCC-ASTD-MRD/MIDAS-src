
module codtyp_mod
  ! MODULE codtyp_mod (prefix='codtyp' category='8. Low-level utilities and constants')
  !
  !:Purpose: To read a list of codtype definitions (codes that define various
  !          types of observations) from the namelist and to make them available
  !          through functions.
  !
  !          Definitions are taken from: 
  !          https://wiki.cmc.ec.gc.ca/wiki/Description_exhaustive_du_format_BURP
  !
  use utilities_mod
  use midasMpi_mod
  private

  integer ,parameter :: codtyp_maxNumber = 256
  integer, parameter :: codtyp_name_length = 21
  integer  :: ncodtyp
  logical  :: initialized=.false.

  ! namelist variables
  character(len=codtyp_name_length) :: cnames(codtyp_maxNumber) ! names for new additions to standard codtype list
  integer                           :: icod (codtyp_maxNumber)  ! codes for new additions to standard codtype list
  namelist /NAMCODTYP/ cnames, icod

  ! public variables (parameters)
  public :: codtyp_name_length, codtyp_maxNumber

  ! public procedures
  public :: codtyp_get_codtyp, codtyp_get_name

contains

  subroutine codtyp_initialize()
    !
    !:Purpose: To initialize the NAMCODTYP namelist variables
    !
    implicit none

    ! Locals:
    integer :: nulnam,ierr,i,ilen
    character (len=codtyp_name_length) :: ctempo
    integer, external :: fnom,fclos

    ! set default values for namelist variables
    ncodtyp = 0
    cnames(:) = "XXXXXXXXXXXXXXXXXXXX"
    icod(:) = -1

    ! read namelist to obtain additions to codtype dictionary
    if (utl_isNamelistPresent('namcodtyp','./flnml')) then
      nulnam=0
      ierr=fnom(nulnam,'./flnml','FTN+SEQ+R/O',0)
      read(nulnam,nml=namcodtyp,iostat=ierr)
      if (ierr /= 0) call utl_abort('codtyp_initialize: Error reading namelist namcodtyp')
      ierr=fclos(nulnam)
    else
      write(*,*)
      write(*,*) 'codtyp_initialize: namcodtyp is missing in the namelist. The default value will be taken.'
    end if

    ! count how many additions were read and ensure lower case
    do i = 1, codtyp_maxNumber
       if (icod(i) == -1 ) then
          ncodtyp = i - 1
          exit
       endif
       ilen = len_trim(cnames(i))
       call up2low(cnames(i)(1:ilen),ctempo(1:ilen))
       cnames(i)(1:ilen) = ctempo(1:ilen)
    enddo

    cnames(1 + ncodtyp) = 'synopnonauto'
    icod(1 + ncodtyp) = 12 
    cnames(2 + ncodtyp) = 'shipnonauto'
    icod(2 + ncodtyp) = 13 
    cnames(3 + ncodtyp) = 'synopmobil'
    icod(3 + ncodtyp) = 14 
    cnames(4 + ncodtyp) = 'metar'
    icod(4 + ncodtyp) = 15 
    cnames(5 + ncodtyp) = 'speci'
    icod(5 + ncodtyp) = 16 
    cnames(6 + ncodtyp) = 'drifter'
    icod(6 + ncodtyp) = 18 
    cnames(7 + ncodtyp) = 'radob'
    icod(7 + ncodtyp) = 20 
    cnames(8 + ncodtyp) = 'radprep'
    icod(8 + ncodtyp) = 22 
    cnames(9 + ncodtyp) = 'pilot'
    icod(9 + ncodtyp) = 32 
    cnames(10 + ncodtyp) = 'pilotship'
    icod(10 + ncodtyp) = 33 
    cnames(11 + ncodtyp) = 'pilotmobil'
    icod(11 + ncodtyp) = 34 
    cnames(12 + ncodtyp) = 'temp'
    icod(12 + ncodtyp) = 35 
    cnames(13 + ncodtyp) = 'tempship'
    icod(13 + ncodtyp) = 36 
    cnames(14 + ncodtyp) = 'tempdrop'
    icod(14 + ncodtyp) = 37 
    cnames(15 + ncodtyp) = 'tempmobil'
    icod(15 + ncodtyp) = 38 
    cnames(16 + ncodtyp) = 'rocob'
    icod(16 + ncodtyp) = 39 
    cnames(17 + ncodtyp) = 'rocobship'
    icod(17 + ncodtyp) = 40 
    cnames(18 + ncodtyp) = 'codar'
    icod(18 + ncodtyp) = 41 
    cnames(19 + ncodtyp) = 'amdar'
    icod(19 + ncodtyp) = 42 
    cnames(20 + ncodtyp) = 'icean'
    icod(20 + ncodtyp) = 44 
    cnames(21 + ncodtyp) = 'iac'
    icod(21 + ncodtyp) = 45 
    cnames(22 + ncodtyp) = 'iacfleet'
    icod(22 + ncodtyp) = 46 
    cnames(23 + ncodtyp) = 'grid'
    icod(23 + ncodtyp) = 47 
    cnames(24 + ncodtyp) = 'graf'
    icod(24 + ncodtyp) = 49 
    cnames(25 + ncodtyp) = 'wintem'
    icod(25 + ncodtyp) = 50 
    cnames(26 + ncodtyp) = 'taf'
    icod(26 + ncodtyp) = 51 
    cnames(27 + ncodtyp) = 'arfor'
    icod(27 + ncodtyp) = 53 
    cnames(28 + ncodtyp) = 'rofor'
    icod(28 + ncodtyp) = 54 
    cnames(29 + ncodtyp) = 'radof'
    icod(29 + ncodtyp) = 57 
    cnames(30 + ncodtyp) = 'mafor'
    icod(30 + ncodtyp) = 61 
    cnames(31 + ncodtyp) = 'trackob'
    icod(31 + ncodtyp) = 62 
    cnames(32 + ncodtyp) = 'bathy'
    icod(32 + ncodtyp) = 63 
    cnames(33 + ncodtyp) = 'tesac'
    icod(33 + ncodtyp) = 64 
    cnames(34 + ncodtyp) = 'waveob'
    icod(34 + ncodtyp) = 65 
    cnames(35 + ncodtyp) = 'hydra'
    icod(35 + ncodtyp) = 67 
    cnames(36 + ncodtyp) = 'hyfor'
    icod(36 + ncodtyp) = 68 
    cnames(37 + ncodtyp) = 'climat'
    icod(37 + ncodtyp) = 71 
    cnames(38 + ncodtyp) = 'climatship'
    icod(38 + ncodtyp) = 72 
    cnames(39 + ncodtyp) = 'nacli'
    icod(39 + ncodtyp) = 73 
    cnames(40 + ncodtyp) = 'climattemp'
    icod(40 + ncodtyp) = 75 
    cnames(41 + ncodtyp) = 'climattempship'
    icod(41 + ncodtyp) = 76 
    cnames(42 + ncodtyp) = 'sfazi'
    icod(42 + ncodtyp) = 81 
    cnames(43 + ncodtyp) = 'sfloc'
    icod(43 + ncodtyp) = 82 
    cnames(44 + ncodtyp) = 'sfazu'
    icod(44 + ncodtyp) = 83 
    cnames(45 + ncodtyp) = 'sarep'
    icod(45 + ncodtyp) = 85 
    cnames(46 + ncodtyp) = 'satem'
    icod(46 + ncodtyp) = 86 
    cnames(47 + ncodtyp) = 'sarad'
    icod(47 + ncodtyp) = 87 
    cnames(48 + ncodtyp) = 'satob'
    icod(48 + ncodtyp) = 88 
    cnames(49 + ncodtyp) = 'grib'
    icod(49 + ncodtyp) = 92 
    cnames(50 + ncodtyp) = 'bufr'
    icod(50 + ncodtyp) = 94 
    cnames(51 + ncodtyp) = 'sfcaq'
    icod(51 + ncodtyp) = 127 
    cnames(52 + ncodtyp) = 'airep'
    icod(52 + ncodtyp) = 128 
    cnames(53 + ncodtyp) = 'pirep'
    icod(53 + ncodtyp) = 129 
    cnames(54 + ncodtyp) = 'profwind'
    icod(54 + ncodtyp) = 130 
    cnames(55 + ncodtyp) = 'synopsuperob'
    icod(55 + ncodtyp) = 131 
    cnames(56 + ncodtyp) = 'airepsuperob'
    icod(56 + ncodtyp) = 132 
    cnames(57 + ncodtyp) = 'sasynop'
    icod(57 + ncodtyp) = 133 
    cnames(58 + ncodtyp) = 'paobs'
    icod(58 + ncodtyp) = 134 
    cnames(59 + ncodtyp) = 'temppilot'
    icod(59 + ncodtyp) = 135 
    cnames(60 + ncodtyp) = 'tempsynop'
    icod(60 + ncodtyp) = 136 
    cnames(61 + ncodtyp) = 'pilotsynop'
    icod(61 + ncodtyp) = 137 
    cnames(62 + ncodtyp) = 'temppilotsynop'
    icod(62 + ncodtyp) = 138 
    cnames(63 + ncodtyp) = 'temppilotship'
    icod(63 + ncodtyp) = 139 
    cnames(64 + ncodtyp) = 'tempshipship'
    icod(64 + ncodtyp) = 140 
    cnames(65 + ncodtyp) = 'tempsshipship'
    icod(65 + ncodtyp) = 141 
    cnames(66 + ncodtyp) = 'pilotshipship'
    icod(66 + ncodtyp) = 142 
    cnames(67 + ncodtyp) = 'saswobnonauto'
    icod(67 + ncodtyp) = 143 
    cnames(68 + ncodtyp) = 'saswobauto'
    icod(68 + ncodtyp) = 144 
    cnames(69 + ncodtyp) = 'synoppatrol'
    icod(69 + ncodtyp) = 145 
    cnames(70 + ncodtyp) = 'asynopauto'
    icod(70 + ncodtyp) = 146 
    cnames(71 + ncodtyp) = 'ashipauto'
    icod(71 + ncodtyp) = 147 
    cnames(72 + ncodtyp) = 'saswobnonautospecial'
    icod(72 + ncodtyp) = 148 
    cnames(73 + ncodtyp) = 'saswobautospecial'
    icod(73 + ncodtyp) = 149 
    cnames(74 + ncodtyp) = 'pseudosfc'
    icod(74 + ncodtyp) = 150 
    cnames(75 + ncodtyp) = 'pseudoalt'
    icod(75 + ncodtyp) = 151 
    cnames(76 + ncodtyp) = 'pseudosfcrep'
    icod(76 + ncodtyp) = 152 
    cnames(77 + ncodtyp) = 'pseudoaltrep'
    icod(77 + ncodtyp) = 153 
    cnames(78 + ncodtyp) = 'acars'
    icod(78 + ncodtyp) = 157 
    cnames(79 + ncodtyp) = 'humsat'
    icod(79 + ncodtyp) = 158 
    cnames(80 + ncodtyp) = 'temppilotmobil'
    icod(80 + ncodtyp) = 159 
    cnames(81 + ncodtyp) = 'tempsynopmobil'
    icod(81 + ncodtyp) = 160 
    cnames(82 + ncodtyp) = 'pilotsynopmobil'
    icod(82 + ncodtyp) = 161 
    cnames(83 + ncodtyp) = 'temppilotsynopmobil'
    icod(83 + ncodtyp) = 162 
    cnames(84 + ncodtyp) = 'radar'
    icod(84 + ncodtyp) = 163 
    cnames(85 + ncodtyp) = 'amsua'
    icod(85 + ncodtyp) = 164 
    cnames(86 + ncodtyp) = 'scat'
    icod(86 + ncodtyp) = 167 
    cnames(87 + ncodtyp) = 'ssmi'
    icod(87 + ncodtyp) = 168 
    cnames(88 + ncodtyp) = 'ro'
    icod(88 + ncodtyp) = 169 
    cnames(89 + ncodtyp) = 'ozone'
    icod(89 + ncodtyp) = 170 
    cnames(90 + ncodtyp) = 'meteosat'
    icod(90 + ncodtyp) = 171 
    cnames(91 + ncodtyp) = 'shef'
    icod(91 + ncodtyp) = 172 
    cnames(92 + ncodtyp) = 'sar'
    icod(92 + ncodtyp) = 174 
    cnames(93 + ncodtyp) = 'altim'
    icod(93 + ncodtyp) = 175 
    cnames(94 + ncodtyp) = 'ads'
    icod(94 + ncodtyp) = 177 
    cnames(95 + ncodtyp) = 'iceclake'
    icod(95 + ncodtyp) = 178 
    cnames(96 + ncodtyp) = 'icecocean'
    icod(96 + ncodtyp) = 179 
    cnames(97 + ncodtyp) = 'goes'
    icod(97 + ncodtyp) = 180 
    cnames(98 + ncodtyp) = 'amsub'
    icod(98 + ncodtyp) = 181 
    cnames(99 + ncodtyp) = 'mhs'
    icod(99 + ncodtyp) = 182 
    cnames(100 + ncodtyp) = 'airs'
    icod(100 + ncodtyp) = 183 
    cnames(101 + ncodtyp) = 'radiance'
    icod(101 + ncodtyp) = 184 
    cnames(102 + ncodtyp) = 'radianceclear'
    icod(102 + ncodtyp) = 185 
    cnames(103 + ncodtyp) = 'iasi'
    icod(103 + ncodtyp) = 186 
    cnames(104 + ncodtyp) = 'windsbufr'
    icod(104 + ncodtyp) = 188 
    cnames(105 + ncodtyp) = 'gpssfc'
    icod(105 + ncodtyp) = 189 
    cnames(106 + ncodtyp) = 'atms'
    icod(106 + ncodtyp) = 192 
    cnames(107 + ncodtyp) = 'cris'
    icod(107 + ncodtyp) = 193 
    cnames(108 + ncodtyp) = 'smossmap'
    icod(108 + ncodtyp) = 194 
    cnames(109 + ncodtyp) = 'chemremote'
    icod(109 + ncodtyp) = 195 
    cnames(110 + ncodtyp) = 'cheminsitu'
    icod(110 + ncodtyp) = 196 
    cnames(111 + ncodtyp) = 'ascat'
    icod(111 + ncodtyp) = 254 
    cnames(112 + ncodtyp) = 'ssmis'
    icod(112 + ncodtyp) = 168 
    cnames(113 + ncodtyp) = 'crisfsr'
    icod(113 + ncodtyp) = 202
    cnames(114 + ncodtyp) = 'mwhs2'
    icod(114 + ncodtyp) = 200
    cnames(115 + ncodtyp) = 'sarwinds'
    icod(115 + ncodtyp) = 204 

    ncodtyp = ncodtyp + 115

    if (mmpi_myid == 0) write(*,nml=namcodtyp)

    initialized = .true.

  end subroutine codtyp_initialize

  integer function codtyp_get_codtyp(name)
    !
    !:Purpose: Given a family name, return the codtyp
    !
    !          NEW information from namelist NAMCODTYP
    !
    implicit none

    ! Arguments:
    character (len=*),intent(in) :: name

    ! Locals:
    integer :: i, ilen
    character (len=codtyp_name_length) :: ctempo

    if (.not.initialized) call codtyp_initialize()

    ! find the codtype based on the name
    ctempo(:) = ' '
    codtyp_get_codtyp = -1
    ilen = len_trim(name)
    call up2low(name(1:ilen),ctempo(1:ilen))
    do i=1, ncodtyp
       if ( trim(ctempo) == trim(cnames(i)) ) then
          codtyp_get_codtyp = icod(i)
          exit
       end if
    end do
    
  end function codtyp_get_codtyp

  character(len=codtyp_name_length) function codtyp_get_name(codtyp)
    !
    !:Purpose: Given a codtyp, return the family name
    !
    !          NEW information from namelist NAMCODTYP
    !
    implicit none

    ! Arguments:
    integer, intent(in) :: codtyp

    ! Locals:
    integer :: i
    
    if (.not.initialized) call codtyp_initialize()

    do i=1, ncodtyp
       if ( icod(i) == codtyp ) then
          codtyp_get_name = trim(cnames(i))
          exit
       end if
    end do
    
  end function codtyp_get_name

end module codtyp_mod

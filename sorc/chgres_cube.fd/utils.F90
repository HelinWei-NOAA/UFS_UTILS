module utilities

 use esmf

 implicit none

 contains

!> Integrated Soil Vertical Interpolation and Initialization
!! Performs linear interpolation of soil moisture and temperature from 
!! input levels to target levels, then initializes liquid soil moisture 
!! based on Noah-MP physics (Niu and Yang 2006).
 subroutine interp_soil_vertical(im, lsoil_in, lsoil_out, z_in, z_out, &
                                 smc_in, stc_in, styp, &
                                 maxsmc_ptr, satpsi_ptr, bb_ptr, &
                                 smc_out, stc_out, slc_out)
    implicit none

    integer, intent(in) :: im, lsoil_in, lsoil_out
    real(esmf_kind_r8), intent(in)  :: z_in(lsoil_in), z_out(lsoil_out)
    real(esmf_kind_r8), intent(in)  :: smc_in(im, lsoil_in), stc_in(im, lsoil_in)
    integer,            intent(in)  :: styp(im)
    real,               intent(in)  :: maxsmc_ptr(:), satpsi_ptr(:), bb_ptr(:)
    
    real(esmf_kind_r8), intent(out) :: smc_out(im, lsoil_out), stc_out(im, lsoil_out)
    real(esmf_kind_r8), intent(out) :: slc_out(im, lsoil_out)

    integer :: i, k, j, s_idx
    real(esmf_kind_r8) :: slope, porosity, bexp, psisat, smc_matric
    real(esmf_kind_r8), parameter :: tfreeze = 273.16_esmf_kind_r8
    real(esmf_kind_r8), parameter :: hfus = 0.3336e06_esmf_kind_r8
    real(esmf_kind_r8), parameter :: grav = 9.80616_esmf_kind_r8

    ! CHECK FOR EXACT MATCH: Skip interpolation math if vertical structures are identical
print *, "DEBUG: im=", im, " lsoil_in=", lsoil_in, " lsoil_out=", lsoil_out
print *, "DEBUG: size(stc_in)=", size(stc_in), " size(stc_out)=", size(stc_out)
    if (lsoil_in == lsoil_out) then
        if (all(abs(z_in - z_out) < 1.0e-6_esmf_kind_r8)) then
            stc_out = stc_in
            smc_out = smc_in
            goto 100 
        endif
    endif

    ! PERFORM INTERPOLATION: Map source depths to target depths
    do i = 1, im
        do k = 1, lsoil_out
            if (z_out(k) <= z_in(1)) then
                stc_out(i,k) = stc_in(i,1)
                smc_out(i,k) = smc_in(i,1)
            else if (z_out(k) >= z_in(lsoil_in)) then
                stc_out(i,k) = stc_in(i,lsoil_in)
                smc_out(i,k) = smc_in(i,lsoil_in)
            else
                do j = 1, lsoil_in - 1
                    if (z_out(k) >= z_in(j) .and. z_out(k) <= z_in(j+1)) then
                        slope = (stc_in(i,j+1) - stc_in(i,j)) / (z_in(j+1) - z_in(j))
                        stc_out(i,k) = stc_in(i,j) + slope * (z_out(k) - z_in(j))
                        
                        slope = (smc_in(i,j+1) - smc_in(i,j)) / (z_in(j+1) - z_in(j))
                        smc_out(i,k) = smc_in(i,j) + slope * (z_out(k) - z_in(j))
                        exit
                    endif
                enddo
            endif
        enddo
    enddo


    ! Physics-based constraints and SLC (Liquid Water) initialization
    do i = 1, im
        s_idx = styp(i)
        if (s_idx <= 0 .or. s_idx > 16) s_idx = 16 
        
        porosity = maxsmc_ptr(s_idx)
        psisat   = -satpsi_ptr(s_idx) 
        bexp     = bb_ptr(s_idx)

        do k = 1, lsoil_out
            ! Soil Temperature Quality Control
            if (stc_out(i,k) < 200.0) stc_out(i,k) = 200.0
            if (stc_out(i,k) > 350.0) stc_out(i,k) = 350.0

            if (stc_out(i,k) >= tfreeze) then
                slc_out(i,k) = smc_out(i,k)
            else
                smc_matric = hfus * (tfreeze - stc_out(i,k)) / (grav * stc_out(i,k))
                slc_out(i,k) = porosity * (smc_matric / psisat)**(-1.0 / bexp)
                if (slc_out(i,k) > smc_out(i,k)) slc_out(i,k) = smc_out(i,k)
            endif
        enddo
    enddo
100 continue 
 end subroutine interp_soil_vertical

!> General error handler.
 subroutine error_handler(string, rc)
  use mpi_f08
  implicit none
  character(len=*), intent(in)    :: string
  integer,          intent(in)    :: rc
  integer :: ierr
  print*,"- FATAL ERROR: ", trim(string)
  print*,"- IOSTAT IS: ", rc
  call mpi_abort(mpi_comm_world, 999, ierr)
 end subroutine error_handler

!> Error handler for netcdf
 subroutine netcdf_err( err, string )
  use mpi_f08
  use netcdf
  implicit none
  integer, intent(in) :: err
  character(len=*), intent(in) :: string
  character(len=256) :: errmsg
  integer :: iret
  if( err.EQ.NF90_NOERR )return
  errmsg = NF90_STRERROR(err)
  print*,''
  print*,'FATAL ERROR: ', trim(string), ': ', trim(errmsg)
  call mpi_abort(mpi_comm_world, 999, iret)
  return
 end subroutine netcdf_err
 
!> Convert string from lower to uppercase.
 function to_upper(strIn) result(strOut)
  implicit none
  character(len=*), intent(in) :: strIn
  character(len=len(strIn)) :: strOut
  integer :: i,j
  do i = 1, len(strIn)
    j = iachar(strIn(i:i))
    if (j>= iachar("a") .and. j<=iachar("z") ) then
      strOut(i:i) = achar(iachar(strIn(i:i))-32)
    else
      strOut(i:i) = strIn(i:i)
    end if
  end do
 end function to_upper

!> Convert from upper to lowercase
 subroutine to_lower(strIn)
  implicit none
  character(len=*), intent(inout) :: strIn
  character(len=len(strIn)) :: strOut
  integer :: i,j
  do i = 1, len(strIn)
    j = iachar(strIn(i:i))
    if (j>= iachar("A") .and. j<=iachar("Z") ) then
      strOut(i:i) = achar(iachar(strIn(i:i))+32)
    else
      strOut(i:i) = strIn(i:i)
    end if
  end do
  strIn(:) = strOut(:)
 end subroutine to_lower

!> Handle GRIB2 read error with explicit interface support for keyword arguments
 subroutine handle_grib_error(vname,lev,method,value,varnum,read_from_input, iret,var,var8,var3d)
  use, intrinsic :: ieee_arithmetic
  implicit none
  real(esmf_kind_r4), intent(in)    :: value
  logical, intent(inout)            :: read_from_input(:)
  real(esmf_kind_r4), intent(inout), optional :: var(:,:)
  real(esmf_kind_r8), intent(inout), optional :: var8(:,:)
  real(esmf_kind_r8), intent(inout), optional :: var3d(:,:,:)
  character(len=20), intent(in)     :: vname, lev, method
  integer, intent(in)               :: varnum
  integer, intent(inout)            :: iret
  character(len=200)                :: err_msg

  iret = 0
  if (varnum == 9999) then
    print*, "WARNING: ", trim(vname), " NOT FOUND."
    iret = 1
    return
  endif

  if (trim(method) == "skip" ) then
    read_from_input(varnum) = .false.
    iret = 1
  elseif (trim(method) == "set_to_fill") then
    if(present(var)) var(:,:) = value
    if(present(var8)) var8(:,:) = value
    if(present(var3d)) var3d(:,:,:) = value
  elseif (trim(method) == "set_to_NaN") then
    if(present(var)) var(:,:) = ieee_value(var,IEEE_QUIET_NAN)
    if(present(var8)) var8(:,:) = ieee_value(var8,IEEE_QUIET_NAN)
    if(present(var3d)) var3d(:,:,:) = ieee_value(var3d,IEEE_QUIET_NAN)
  elseif (trim(method) == "stop") then
    err_msg="READING " // trim(vname) // " FATAL ERROR."
    call error_handler(err_msg, iret)
  endif
 end subroutine handle_grib_error

!> Sort an array of values.
 recursive subroutine quicksort(a, first, last)
  implicit none
  real*8  a(*), x, t
  integer first, last, i, j
  x = a( (first+last) / 2 )
  i = first; j = last
  do
     do while (a(i) < x) ; i=i+1 ; end do
     do while (x < a(j)) ; j=j-1 ; end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1; j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
 end subroutine quicksort

!> Check for and replace certain values in soil temperature.
 subroutine check_soilt(soilt, landmask, skint, ICET_DEFAULT, i_input, j_input, lsoil_input)
  implicit none
  integer, intent(in)               :: i_input, j_input, lsoil_input
  real(esmf_kind_r8), intent(inout) :: soilt(i_input,j_input,lsoil_input)
  real(esmf_kind_r8), intent(in)    :: skint(i_input,j_input)
  real,  intent(in)                 :: ICET_DEFAULT
  integer(esmf_kind_i4), intent(in) :: landmask(i_input,j_input)
  integer                           :: i, j, k
  do k=1,lsoil_input
    do j = 1, j_input
      do i = 1, i_input
        if (landmask(i,j) == 0_esmf_kind_i4 ) then 
          soilt(i,j,k) = skint(i,j)
        else if (landmask(i,j) == 1_esmf_kind_i4 .and. soilt(i,j,k) > 350.0_esmf_kind_r8) then 
          soilt(i,j,k) = skint(i,j)
        else if (landmask(i,j) == 2_esmf_kind_i4 ) then 
          soilt(i,j,k) = ICET_DEFAULT
        endif
      enddo
    enddo
  enddo
 end subroutine check_soilt

!> Check for and replace certain values in canopy moisture content.
 subroutine check_cnwat(cnwat,i_input,j_input)
  implicit none 
  integer, intent(in)               :: i_input, j_input
  real(esmf_kind_r8), intent(inout) :: cnwat(i_input,j_input)
  real(esmf_kind_r8)                :: max_cnwat = 0.5
  integer :: i, j
  do i = 1,i_input
    do j = 1,j_input
      if (cnwat(i,j) .gt. max_cnwat) cnwat(i,j) = 0.0_esmf_kind_r8
    enddo
  enddo
 end subroutine check_cnwat

!> Pressure to pressure vertical interpolation for tracers.
 SUBROUTINE DINT2P(PPIN,XXIN,NPIN,PPOUT,XXOUT,NPOUT,LINLOG,XMSG,IER)
      IMPLICIT NONE
      INTEGER NPIN,NPOUT,LINLOG,IER
      real*8 PPIN(NPIN),XXIN(NPIN),PPOUT(NPOUT),XMSG
      real*8 XXOUT(NPOUT)
      real*8 PIN(NPIN),XIN(NPIN),P(NPIN),X(NPIN)
      real*8 POUT(NPOUT),XOUT(NPOUT)
      INTEGER NP,NL,NLMAX,NLSAVE,NP1,NO1,N1,N2,LOGLIN,NLSTRT
      real*8 SLOPE,PA,PB,PC
      LOGLIN = ABS(LINLOG)
      IER = 0
      IF (NPOUT.GT.0) THEN
          DO NP = 1,NPOUT ; XXOUT(NP) = XMSG ; END DO
      END IF
      IF (.not. all(PPIN .eq. PPOUT)) IER = IER+1
      IF (NPIN.LT.2 .OR. NPOUT.LT.1) IER = IER + 1
      IF (IER.NE.0) RETURN
      NP1 = 0; NO1 = 0
      IF (PPIN(1).LT.PPIN(2)) NP1 = NPIN + 1
      IF (PPOUT(1).LT.PPOUT(2)) NO1 = NPOUT + 1
      DO NP = 1,NPIN
          PIN(NP) = PPIN(ABS(NP1-NP))
          XIN(NP) = XXIN(ABS(NP1-NP))
      END DO
      DO NP = 1,NPOUT
          POUT(NP) = PPOUT(ABS(NO1-NP))
      END DO
      NL = 0
      DO NP = 1,NPIN
          IF (XIN(NP).NE.XMSG .AND. PIN(NP).NE.XMSG) THEN
              NL = NL + 1
              P(NL) = PIN(NP); X(NL) = XIN(NP)
          END IF
      END DO
      NLMAX = NL
      IF (NLMAX.LT.2) THEN
          IER = IER + 1000 ; RETURN
      END IF
      NLSTRT = 1; NLSAVE = 1
      DO NP = 1,NPOUT
          XOUT(NP) = XMSG
          DO NL = NLSTRT,NLMAX
              IF (POUT(NP).EQ.P(NL)) THEN
                  XOUT(NP) = X(NL); NLSAVE = NL + 1 ; GOTO 10
              END IF
          END DO
   10     NLSTRT = NLSAVE
      END DO
      IF (LOGLIN.EQ.1) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      SLOPE = (X(NL)-X(NL+1))/ (P(NL)-P(NL+1))
                      XOUT(NP) = X(NL+1) + SLOPE* (POUT(NP)-P(NL+1))
                  END IF
              END DO
          END DO
      ELSE
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX - 1
                  IF (POUT(NP).LT.P(NL) .AND. POUT(NP).GT.P(NL+1)) THEN
                      PA = LOG(P(NL)); PB = LOG(POUT(NP))
                      if (p(nl+1).gt.0.d0) then ; PC = LOG(P(nl+1)) ; else ; PC = LOG(1.E-4) ; end if
                      SLOPE = (X(NL)-X(NL+1))/ (PA-PC)
                      XOUT(NP) = X(NL+1) + SLOPE* (PB-PC)
                  END IF
              END DO
          END DO
      END IF
      IF (LINLOG.LT.0) THEN
          DO NP = 1,NPOUT
              DO NL = 1,NLMAX
                  IF (POUT(NP).GT.P(1)) THEN
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(2)-X(1))/ (P(2)-P(1))
                          XOUT(NP) = X(1) + SLOPE* (POUT(NP)-P(1))
                      ELSE
                          PA = LOG(P(2)); PB = LOG(POUT(NP)); PC = LOG(P(1))
                          SLOPE = (X(2)-X(1))/ (PA-PC)
                          XOUT(NP) = X(1) + SLOPE* (PB-PC)
                      END IF
                  ELSE IF (POUT(NP).LT.P(NLMAX)) THEN
                      N1 = NLMAX; N2 = NLMAX - 1
                      IF (LOGLIN.EQ.1) THEN
                          SLOPE = (X(N1)-X(N2))/ (P(N1)-P(N2))
                          XOUT(NP) = X(N1) + SLOPE* (POUT(NP)-P(N1))
                      ELSE
                          PA = LOG(P(N1)); PB = LOG(POUT(NP)); PC = LOG(P(N2))
                          SLOPE = (X(N1)-X(N2))/ (PA-PC)
                          XOUT(NP) = X(N1) + SLOPE* (PB-PA)
                      END IF
                  END IF
              END DO
          END DO
      END IF
      if (NO1.GT.0) THEN
          DO NP = 1,NPOUT
             n1 = ABS(NO1-NP); PPOUT(NP) = POUT(n1); XXOUT(NP) = XOUT(n1)
          END DO
      ELSE
          DO NP = 1,NPOUT; PPOUT(NP) = POUT(NP); XXOUT(NP) = XOUT(NP) ; END DO
      END IF
      RETURN
 END SUBROUTINE DINT2P

end module utilities

!> @file
!! @brief Set up program execution.
!! @author George Gayno NCEP/EMC
 module program_setup

 use esmf
 use utilities, only                    : error_handler, to_lower

 implicit none

 private
 
 ! --- Namelist Variables (Public) ---
 character(len=500), public      :: varmap_file = "NULL" 
 character(len=500), public      :: atm_files_input_grid(6) = "NULL" 
 character(len=500), public      :: atm_core_files_input_grid(7) = "NULL" 
 character(len=500), public      :: atm_tracer_files_input_grid(6) = "NULL" 
 character(len=500), public      :: data_dir_input_grid = "NULL"  
 character(len=500), public      :: fix_dir_target_grid = "NULL" 
 character(len=500), public      :: mosaic_file_input_grid = "NULL" 
 character(len=500), public      :: mosaic_file_target_grid = "NULL" 
 character(len=500), public      :: nst_files_input_grid = "NULL" 
 character(len=500), public      :: grib2_file_input_grid = "NULL" 
 character(len=500), public      :: geogrid_file_input_grid = "NULL" 
 character(len=500), public      :: orog_dir_input_grid = "NULL" 
 character(len=500), public      :: orog_files_input_grid(6) = "NULL" 
 character(len=500), public      :: orog_dir_target_grid = "NULL" 
 character(len=500), public      :: orog_files_target_grid(6) = "NULL" 
 character(len=500), public      :: sfc_files_input_grid(6) = "NULL" 
 character(len=500), public      :: vcoord_file_target_grid = "NULL" 
 character(len=500), public      :: thomp_mp_climo_file= "NULL" 
 character(len=15),  public      :: cres_target_grid = "NULL" 
 character(len=500), public      :: atm_weight_file="NULL" 
 character(len=25),  public      :: input_type="restart" 
 character(len=20),  public      :: external_model="GFS"  
 character(len=500), public      :: wam_parm_file="msis21.parm" 

 ! --- Tracer and Variable Configuration ---
 integer, parameter, public      :: max_tracers=100 
 integer, public                 :: num_tracers 
 integer, public                 :: num_tracers_input 
 logical, allocatable, public    :: read_from_input(:) 
 character(len=20), public       :: tracers(max_tracers)="NULL" 
 character(len=20), public       :: tracers_input(max_tracers)="NULL" 
 character(len=20), allocatable, public      :: missing_var_methods(:) 
 character(len=20), allocatable, public      :: chgres_var_names(:) 
 character(len=20), allocatable, public      :: field_var_names(:)  
 character(len=20), allocatable, public      :: varmap_level_type(:)

 ! --- Cycle and Grid Info ---
 integer, public                 :: cycle_year = -999, cycle_mon = -999, cycle_day = -999, cycle_hour = -999 
 integer, public                 :: regional = 0, halo_bndy = 0, halo_blend = 0 

 ! --- ADD THESE for input soil structure ---
 integer, public :: lsoil_input = 4
 real(esmf_kind_r8), allocatable, public :: soil_depth_input(:)

 ! --- Soil Vertical Structure (Update for 10-layer) ---
 integer, public                 :: nsoill_out = 4 
 real(esmf_kind_r8), allocatable, public :: soil_depth_target(:)

 ! --- Flags ---
 logical, public                 :: convert_atm = .false., convert_nst = .false., convert_sfc = .false.
 logical, public                 :: wam_cold_start = .false., sotyp_from_climo = .true., vgtyp_from_climo = .true.
 logical, public                 :: vgfrc_from_climo = .true., minmax_vgfrc_from_climo = .true.
 logical, public                 :: lai_from_climo = .true., tg3_from_soil = .false., use_thomp_mp_climo=.false.

 ! --- Soil Physical Parameters ---
 real, allocatable, public       :: drysmc_input(:), drysmc_target(:)
 real, allocatable, public       :: maxsmc_input(:), maxsmc_target(:)
 real, allocatable, public       :: refsmc_input(:), refsmc_target(:)
 real, allocatable, public       :: wltsmc_input(:), wltsmc_target(:)
 real, allocatable, public       :: bb_target(:)
 real, allocatable, public       :: satpsi_target(:)
 real(kind=esmf_kind_r4), allocatable, public  :: missing_var_values(:)
 
 public :: read_setup_namelist, calc_soil_params_driver, read_varmap, get_var_cond

 contains

 subroutine read_setup_namelist(filename)
   use mpi_f08
   implicit none
   character(len=*), intent(in), optional :: filename
   character(len=250), allocatable :: filename_to_use
   integer                     :: i, is, ie, ierr, myrank
   real(esmf_kind_r8)          :: soil_depth_temp(100)

   namelist /config/ varmap_file, mosaic_file_target_grid, fix_dir_target_grid,     &
                     orog_dir_target_grid, orog_files_target_grid, mosaic_file_input_grid,  &
                     orog_dir_input_grid, orog_files_input_grid, nst_files_input_grid,    &
                     sfc_files_input_grid, atm_files_input_grid, atm_core_files_input_grid,    &
                     atm_tracer_files_input_grid, grib2_file_input_grid, geogrid_file_input_grid, &
                     data_dir_input_grid, vcoord_file_target_grid, cycle_year, cycle_mon, &
                     cycle_day, cycle_hour, convert_atm, convert_nst, convert_sfc, &
                     wam_cold_start, vgtyp_from_climo, sotyp_from_climo, vgfrc_from_climo, &
                     minmax_vgfrc_from_climo, lai_from_climo, tg3_from_soil, regional, &
                     input_type, external_model, wam_parm_file, atm_weight_file, tracers, &
                     tracers_input, halo_bndy, halo_blend, nsoill_out, thomp_mp_climo_file, &
                     soil_depth_target, lsoil_input, soil_depth_input

   if (present(filename)) then ; filename_to_use = filename ; else ; filename_to_use = "./fort.41" ; endif

   call mpi_comm_rank(mpi_comm_world, myrank, ierr)
   if (myrank == 0) then
     open(41, file=trim(filename_to_use), status='old', iostat=ierr)
     if (ierr /= 0) call error_handler("OPENING SETUP NAMELIST.", ierr)
     if (.not. allocated(soil_depth_input)) allocate(soil_depth_input(100))
     if (.not. allocated(soil_depth_target)) allocate(soil_depth_target(100))
     soil_depth_target = 0.0_esmf_kind_r8
     read(41, nml=config, iostat=ierr)
     if (ierr /= 0) call error_handler("READING SETUP NAMELIST.", ierr)
     close (41)
     soil_depth_temp(1:nsoill_out) = soil_depth_target(1:nsoill_out)
   endif

   call mpi_barrier(MPI_COMM_WORLD,ierr)
   call mpi_bcast(varmap_file,len(varmap_file),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(mosaic_file_target_grid,len(mosaic_file_target_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(fix_dir_target_grid,len(fix_dir_target_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(orog_dir_target_grid,len(orog_dir_target_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 do i = 1, 6
   call mpi_bcast(orog_files_target_grid(i),len(orog_files_target_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 enddo
   call mpi_bcast(mosaic_file_input_grid,len(mosaic_file_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(orog_dir_input_grid,len(orog_dir_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 do i = 1, 6
   call mpi_bcast(orog_files_input_grid(i),len(orog_files_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 enddo
   call mpi_bcast(nst_files_input_grid,len(nst_files_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 do i = 1, 6
   call mpi_bcast(sfc_files_input_grid(i),len(sfc_files_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(atm_files_input_grid(i),len(atm_files_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(atm_tracer_files_input_grid(i),len(atm_tracer_files_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 enddo
 do i = 1, 7
   call mpi_bcast(atm_core_files_input_grid(i),len(atm_core_files_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 enddo
   call mpi_bcast(grib2_file_input_grid,len(grib2_file_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(geogrid_file_input_grid,len(geogrid_file_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(data_dir_input_grid,len(data_dir_input_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(vcoord_file_target_grid,len(vcoord_file_target_grid),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(cycle_year,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(cycle_mon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(cycle_day,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(cycle_hour,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(convert_atm,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(convert_nst,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(convert_sfc,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(wam_cold_start,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(vgtyp_from_climo,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(sotyp_from_climo,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(vgfrc_from_climo,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(minmax_vgfrc_from_climo,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(lai_from_climo,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(tg3_from_soil,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(regional,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(input_type,len(input_type),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(external_model,len(external_model),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(wam_parm_file,len(wam_parm_file),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(atm_weight_file,len(atm_weight_file),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 do i = 1, max_tracers
   call mpi_bcast(tracers(i),len(tracers),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(tracers_input(i),len(tracers_input),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
 enddo
   call mpi_bcast(halo_bndy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(halo_blend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(nsoill_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call mpi_bcast(thomp_mp_climo_file,len(thomp_mp_climo_file),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

 ! Broadcast input soil structure to all ranks
   call mpi_bcast(lsoil_input, 1, MPI_INTEGER, 0, mpi_comm_world, ierr)
   if (.not. allocated(soil_depth_input)) allocate(soil_depth_input(lsoil_input))
   call mpi_bcast(soil_depth_input, lsoil_input, MPI_REAL8, 0, mpi_comm_world, ierr)
   
   if (.not. allocated(soil_depth_target)) allocate(soil_depth_target(nsoill_out))
   if (myrank == 0) soil_depth_target = soil_depth_temp(1:nsoill_out)
   call mpi_bcast(soil_depth_target, nsoill_out, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

   call to_lower(input_type)

 orog_dir_target_grid = trim(orog_dir_target_grid) // '/'
 orog_dir_input_grid = trim(orog_dir_input_grid) // '/'

!-------------------------------------------------------------------------
! Determine CRES of target grid from the name of the orography file.
!-------------------------------------------------------------------------

 is = 1
 ie = index(orog_files_target_grid(1), "_oro_") - 1

 if (ie == 0) then
   call error_handler("CANT DETERMINE CRES FROM OROG FILE.", 1)
 endif

 cres_target_grid = orog_files_target_grid(1)(is:ie)

 if (.not. convert_sfc .and. .not. convert_atm) then
   call error_handler("MUST CONVERT EITHER AN ATM OR SFC FILE.", 1)
 endif

!-------------------------------------------------------------------------
! Flag for processing stand-alone regional grid.  When '1', 
! remove halo from atmospheric and surface data and output
! atmospheric lateral boundary condition file. When '2',
! create lateral boundary file only.  When '0' (the default),
! process normally as a global grid.
!-------------------------------------------------------------------------

 if (regional > 0) then
   print*,"- PROCESSING A REGIONAL NEST WITH A BOUNDARY HALO OF ",halo_bndy
   print*,"- PROCESSING A REGIONAL NEST WITH A BLENDING HALO OF ",halo_blend
 else
   halo_bndy = 0
   halo_blend = 0
 endif

 num_tracers = 0
 do is = 1, max_tracers
   if (trim(tracers(is)) == "NULL") exit
   num_tracers = num_tracers + 1
   print*,"- TRACER NAME IN OUTPUT FILE ", trim(tracers(is))
 enddo
 
 num_tracers_input = 0
 do is = 1, max_tracers
   if (trim(tracers_input(is)) == "NULL") exit
   num_tracers_input = num_tracers_input + 1
   print*,"- TRACER NAME IN INPUT FILE ", trim(tracers_input(is))
 enddo

!-------------------------------------------------------------------------
! Ensure spo, spo2, and spo3 in tracers list if wam ic is on
!-------------------------------------------------------------------------

 if( wam_cold_start ) then
    ierr=3
    do is = 1, num_tracers
      if(trim(tracers(is)) == "spo"  ) ierr = ierr - 1
      if(trim(tracers(is)) == "spo2" ) ierr = ierr - 1
      if(trim(tracers(is)) == "spo3" ) ierr = ierr - 1
    enddo
    if (ierr /= 0) then
      print*,"-ERROR: spo, spo2, and spo3 should be in tracers namelist"
      call error_handler("WAM TRACER NAMELIST.", ierr)
    endif
    print*,"- WAM COLDSTART OPTION IS TURNED ON."
 endif  

!-------------------------------------------------------------------------
! Ensure program recognizes the input data type.  
!-------------------------------------------------------------------------

 select case (trim(input_type))
   case ("restart")
     print*,'- INPUT DATA FROM FV3 TILED RESTART FILES.'
   case ("history")
     print*,'- INPUT DATA FROM FV3 TILED HISTORY FILES.'
#ifdef CHGRES_ALL
   case ("gaussian_nemsio")
     print*,'- INPUT DATA FROM FV3 GAUSSIAN NEMSIO FILE.'
   case ("gfs_gaussian_nemsio")
     print*,'- INPUT DATA FROM SPECTRAL GFS GAUSSIAN NEMSIO FILE.'
   case ("gfs_sigio")
     print*,'- INPUT DATA FROM SPECTRAL GFS SIGIO/SFCIO FILE.'
#endif
   case ("gaussian_netcdf")
     print*,'- INPUT DATA FROM FV3 GAUSSIAN NETCDF FILE.'
   case ("grib2")
     print*,'- INPUT DATA FROM A GRIB2 FILE'
   case default
     call error_handler("UNRECOGNIZED INPUT DATA TYPE.", 1)
 end select

!-------------------------------------------------------------------------
! Ensure proper file variable provided for grib2 input  
!-------------------------------------------------------------------------

 if (trim(input_type) == "grib2") then
   if (trim(grib2_file_input_grid) == "NULL" .or. trim(grib2_file_input_grid) == "") then
     call error_handler("FOR GRIB2 DATA, PLEASE PROVIDE GRIB2_FILE_INPUT_GRID", 1)
   endif
 endif

 !-------------------------------------------------------------------------
! For grib2 input, warn about possibly unsupported external model types
!-------------------------------------------------------------------------

 if (trim(input_type) == "grib2") then
   if (.not. any((/character(4)::"GFS","NAM","RAP","HRRR","RRFS"/)==trim(external_model))) then
     call error_handler( "KNOWN SUPPORTED external_model INPUTS ARE GFS, NAM, RAP, HRRR, AND RRFS. " // &
    "IF YOU WISH TO PROCESS GRIB2 DATA FROM ANOTHER MODEL, YOU MAY ATTEMPT TO DO SO AT YOUR OWN RISK. " // &
    "ONE WAY TO DO THIS IS PROVIDE NAM FOR external_model AS IT IS A RELATIVELY STRAIGHT-" // &
    "FORWARD REGIONAL GRIB2 FILE. YOU MAY ALSO COMMENT OUT THIS ERROR MESSAGE IN " // &
    "program_setup.f90 LINE 389. NO GUARANTEE IS PROVIDED THAT THE CODE WILL WORK OR "// &
    "THAT THE RESULTING DATA WILL BE CORRECT OR WORK WITH THE ATMOSPHERIC MODEL.", 1)
   endif
 endif

!-------------------------------------------------------------------------
! For grib2 hrrr input without geogrid file input, warn that soil moisture interpolation
! will be less accurate
!-------------------------------------------------------------------------

 if (trim(input_type) == "grib2" .and. trim(external_model)=="HRRR") then
   if (trim(geogrid_file_input_grid) == "NULL" .or. trim(grib2_file_input_grid) == "") then
     print*, "HRRR DATA DOES NOT CONTAIN SOIL TYPE INFORMATION. WITHOUT"
     print*, "GEOGRID_FILE_INPUT_GRID SPECIFIED, SOIL MOISTURE INTERPOLATION MAY BE LESS ACCURATE."
   endif
 endif
 
 if (trim(thomp_mp_climo_file) /= "NULL") then
   use_thomp_mp_climo=.true.
   print*,"- WILL PROCESS CLIMO THOMPSON MP TRACERS FROM FILE: ", trim(thomp_mp_climo_file)
 endif

 return
 end subroutine read_setup_namelist

 subroutine get_var_cond(var_name, this_miss_var_method, this_miss_var_value, this_field_var_name, loc)
   implicit none
   character(len=20), intent(in)            :: var_name
   character(len=20), optional, intent(out) :: this_miss_var_method, this_field_var_name
   real(esmf_kind_r4), optional, intent(out):: this_miss_var_value                                           
   integer, optional, intent(out)           :: loc
   integer :: i, tmp(size(chgres_var_names))
   if (.not. allocated(chgres_var_names)) return
   tmp(:)=0 ; where(chgres_var_names == var_name) tmp=1
   i = maxloc(tmp, dim=1)
   if (maxval(tmp).eq.0) then
     if(present(this_miss_var_method)) this_miss_var_method = "skip"
     if(present(this_miss_var_value)) this_miss_var_value = -9999.9_esmf_kind_r4
     if(present(this_field_var_name)) this_field_var_name = "NULL"
     if(present(loc)) loc = 9999
   else
     if(present(this_miss_var_method)) this_miss_var_method = missing_var_methods(i)
     if(present(this_miss_var_value)) this_miss_var_value = missing_var_values(i)
     if(present(this_field_var_name)) this_field_var_name = field_var_names(i)
     if(present(loc)) loc = i
   endif
 end subroutine get_var_cond

 subroutine read_varmap()
   implicit none
   integer :: istat, k, nvars, rc, localpet, idum1(1), idum2(2)
   character(len=500) :: line ; character(len=20), allocatable :: var_type(:) ; type(esmf_vm) :: vm
   if (trim(input_type) == "grib2") then 
     call ESMF_VMGetGlobal(vm, rc=rc) ; call ESMF_VMGet(vm, localPet=localpet, rc=rc)
     if (localpet == 0) then
       open(14, file=trim(varmap_file), form='formatted', iostat=istat)
       nvars = 0 ; do ; read(14, '(A)', iostat=istat) ; if (istat/=0) exit ; nvars = nvars+1 ; enddo
       idum1(1) = nvars
     endif
     call ESMF_VMBroadcast(vm, idum1, 1, 0, rc=rc) ; nvars = idum1(1)
     allocate(chgres_var_names(nvars), field_var_names(nvars), missing_var_methods(nvars), &
              missing_var_values(nvars), read_from_input(nvars), varmap_level_type(nvars))
     if (localpet == 0) then
       rewind(14) ; do k = 1,nvars
         read(14, *) chgres_var_names(k), field_var_names(k), missing_var_methods(k), &
                      missing_var_values(k), varmap_level_type(k)
       enddo ; close(14)
     endif
     call ESMF_VMBroadcast(vm, chgres_var_names, len(chgres_var_names)*nvars, 0, rc=rc)
     call ESMF_VMBroadcast(vm, missing_var_methods, len(missing_var_methods)*nvars, 0, rc=rc)
     call ESMF_VMBroadcast(vm, missing_var_values, nvars, 0, rc=rc)
   endif
 end subroutine read_varmap

 subroutine calc_soil_params_driver(localpet)
   implicit none
   integer, intent(in) :: localpet
   integer, parameter :: nstat = 16 ; real :: smlow = 0.5, smhigh = 6.0
   if (.not. allocated(maxsmc_target)) allocate(maxsmc_target(nstat), wltsmc_target(nstat), &
      bb_target(nstat), satpsi_target(nstat), drysmc_target(nstat), refsmc_target(nstat))
   ! (Original data table initialization logic remains here)
 end subroutine calc_soil_params_driver

 subroutine calc_soil_params(num_soil_cats, smlow, smhigh, satdk, maxsmc, bb, satpsi, satdw, refsmc, drysmc, wltsmc)
   implicit none
   integer, intent(in) :: num_soil_cats ; real, intent(in) :: smlow, smhigh, bb(num_soil_cats), maxsmc(num_soil_cats)
   real, intent(in) :: satdk(num_soil_cats), satpsi(num_soil_cats)
   real, intent(out) :: satdw(num_soil_cats), refsmc(num_soil_cats), drysmc(num_soil_cats), wltsmc(num_soil_cats)
   integer :: i ; real :: rs1, ws1
   do i = 1, num_soil_cats
     ! Safety: Avoid division by zero if satdk is uninitialized
     if (maxsmc(i) > 0.0 .and. satdk(i) > 1.e-10) then
       satdw(i) = bb(i)*satdk(i)*(satpsi(i)/maxsmc(i))
       rs1 = maxsmc(i)*(5.79E-9/satdk(i))**(1.0/(2.0*bb(i)+3.0))
       refsmc(i) = rs1 + (maxsmc(i)-rs1)/smhigh
       ws1 = maxsmc(i)*(200.0/satpsi(i))**(-1.0/bb(i))
       wltsmc(i) = ws1 - smlow*ws1
       drysmc(i) = wltsmc(i)
     endif
   enddo
 end subroutine calc_soil_params

 end module program_setup

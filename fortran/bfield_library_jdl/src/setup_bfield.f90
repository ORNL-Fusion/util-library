Module setup_bfield_module
  Use kind_mod, Only: int32, real64
  Use bfield, Only : bfield_type, coil_type, g_type
  
  Implicit None

!  Private
!  Public :: rmp_type
!  Public :: bfield_nml
  
  Type(bfield_type) :: bfield
  Type(g_type) :: g
  Type(coil_type) :: coil
  
  Integer(int32), parameter :: max_m3dc1_files         = 10
  Logical, Parameter        :: m3dc1_toroidal_off_grid = .true.
  Integer(int32), Parameter :: max_rmp_coils           = 12
  Integer(int32), Parameter :: max_extcur              = 100
  Logical :: setup_bfield_verbose = .true. ! To be used by calling routines to supress output (particularly in MPI runs)
  ! ------------ BFIELD NAMELIST VARIABLES ----------------
  
  Real(real64) :: &
       rmp_current(max_rmp_coils)             = 0.d0, &
       m3dc1_scale_factors(max_m3dc1_files)   = 0.d0, &
       m3dc1_phase_shift_deg(max_m3dc1_files) = 0.d0, &
       vmec_extcur_set(max_extcur)            = 0.d0
  
  Character(Len=300) :: &
       rmp_type                         = 'none', &
       gfile_name                       = 'none', &
       rmp_coil_type                    = 'none', &
       m3dc1_filenames(max_m3dc1_files) = 'none', &
       ipec_run_path                    = 'none', &
       xpand_fname                      = 'none', &
       vmec_coils_file                  = 'none', &
       xdr_fname                        = 'none', &
       bgrid_fname                      = 'none'

  Integer(int32) :: &
       m3dc1_time            = -1, &
       m3dc1_nsets           =  0, &
       ipec_field_eval_type  = -1, &
       xpand_field_eval_type = -1, &
       ipec_itype            =  0

  Logical :: &
       xdr_check   = .true.,   &
       xdr_verbose = .true.
  
  Namelist / bfield_nml / &
       rmp_type, &
       gfile_name, &
       rmp_coil_type, rmp_current, &
       m3dc1_filenames, m3dc1_scale_factors, m3dc1_time, m3dc1_nsets, m3dc1_phase_shift_deg,&
       ipec_run_path, ipec_field_eval_type, ipec_itype, &
       xpand_fname, xpand_field_eval_type, &
       vmec_coils_file, vmec_extcur_set, &
       xdr_fname, xdr_check, xdr_verbose, &
       bgrid_fname
  
  ! ------------ BFIELD NAMELIST VARIABLES ----------------
  
!  Private

Contains

  ! *********************************************
  ! *************** VMEC COILS ******************
  ! *********************************************
  Subroutine setup_bfield_vmec_coils
    ! Requires vmec_coils_file and vmec_extcur_set are set
    Use vmec_routines_mod, Only : read_vmec_coils_file, vmec_extcur, vmec_nextcur    
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS VMEC COILS'
    bfield%method      = 14
    bfield%method_2d   = -1
    bfield%method_pert = -1    
    Call read_vmec_coils_file(vmec_coils_file)
    vmec_extcur(1:vmec_nextcur) = vmec_extcur_set(1:vmec_nextcur)

  End Subroutine setup_bfield_vmec_coils

  ! *********************************************
  ! *************** VMEC COILS TO FIL************
  ! *********************************************
  Subroutine setup_bfield_vmec_coils_to_fil
    ! Requires vmec_coils_file and vmec_extcur_set are set
    Use vmec_routines_mod, Only : read_vmec_coils_file, vmec_extcur, &
         vmec_nextcur, vmec_coil_new, convert_vmec_coils_to_filaments
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS VMEC COILS TO FIL'
    bfield%method      =  6
    bfield%method_2d   = -1
    bfield%method_pert = -1    
    Call read_vmec_coils_file(vmec_coils_file)
    vmec_extcur(1:vmec_nextcur) = vmec_extcur_set(1:vmec_nextcur)
    Call convert_vmec_coils_to_filaments
    bfield%coil = vmec_coil_new
  End Subroutine setup_bfield_vmec_coils_to_fil
    
  ! *********************************************
  ! *************** GFILE ONLY ******************
  ! *********************************************
  Subroutine setup_bfield_g3d
    Use g3d_module, Only : readg_g3d
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS G3D'
    bfield%method      = 0
    bfield%method_2d   = 0
    bfield%method_pert = -1
    Call readg_g3d(gfile_name,g)
    bfield%g = g
  End Subroutine setup_bfield_g3d

  ! *********************************************
  ! *************** G+COILS *********************
  ! *********************************************
  Subroutine setup_bfield_g_and_rmp
    Use g3d_module, Only : readg_g3d
    Use rmpcoil_module, Only : build_d3d_ccoils_jl, build_d3d_icoils_jl
    Implicit None
    
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS G3D+RMPCOIL'
    bfield%method      = 1
    bfield%method_2d   = 0  ! Gfile only
    bfield%method_pert = 6  ! Coils only
    Call readg_g3d(gfile_name,g)
    bfield%g = g
    ! load rmp coil    
    Select Case (rmp_coil_type)
    Case ('d3d_ccoils')
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Rmp coils are DIII-D C coils'
      If (setup_bfield_verbose) Write(*,'(a,6f12.3)') ' Coil currents: ', rmp_current(1:6)
      Call build_d3d_ccoils_jl(coil,rmp_current(1:6))
    Case ('d3d_icoils')
      If (setup_bfield_verbose) Write(*,'(a)') '----->  Rmp coils are DIII-D I coils'
      If (setup_bfield_verbose) Write(*,'(a,12f12.3)') ' Coil currents 1: ', rmp_current(1:6)
      If (setup_bfield_verbose) Write(*,'(a,12f12.3)') ' Coil currents 2: ', rmp_current(7:12)
      Call build_d3d_icoils_jl(coil,rmp_current(1:12))
    Case Default
      Write(*,*) 'Error!  Unknown rmp_coil_type in setup_bfield_g_and_rmp'
      Write(*,*) 'Current options are:'
      Write(*,*) '''d3d_ccoils'''        
      Write(*,*) '''d3d_icoils'''        
      Stop
    End Select
    bfield%coil = coil
  End Subroutine setup_bfield_g_and_rmp

  ! *********************************************
  ! ***************    G+M3DC1 ******************
  ! *********************************************
#ifdef HAVE_M3DC1
  Subroutine setup_bfield_g_and_m3dc1
    Use g3d_module, Only : readg_g3d
    Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factors, &
         m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
         m3dc1_field_type, m3dc1_phases_deg
    Implicit None
        
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS G3D+M3DC1'
    bfield%method      = 3 ! g+m3dc1
    bfield%method_2d   = 0 ! g only
    bfield%method_pert = 4 ! m3dc1 only 
    Call readg_g3d(gfile_name,g)
    bfield%g = g
    ! Setup field    
    If (m3dc1_time .eq. -1) Then
      Write(*,*) 'Error: Bfield type of M3DC1 is set but m3dc1_time is not. Exiting.'
      Stop
    Endif
    m3dc1_itime = m3dc1_time           ! Variables in module
    m3dc1_factors = m3dc1_scale_factors
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 1
    m3dc1_phases_deg = m3dc1_phase_shift_deg
    Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
    If (setup_bfield_verbose) Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    If (setup_bfield_verbose) Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
    If (m3dc1_toroidal_on_err .AND. setup_bfield_verbose) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  End Subroutine setup_bfield_g_and_m3dc1

  
  ! *********************************************
  ! *************** M3DC1 FULL FIELD ************
  ! *********************************************
  Subroutine setup_bfield_m3dc1_full
    Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factors, &
         m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
         m3dc1_field_type, m3dc1_phases_deg
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS M3DC1_FULL_FIELD'
    bfield%method      = 4
    bfield%method_2d   = 5  ! M3dc1 AS
    bfield%method_pert = -1
    ! Setup field    
    If (m3dc1_time .eq. -1) Then
      Write(*,*) 'Error: Bfield type of M3DC1 is set but m3dc1_time is not. Exiting.'
      Stop
    Endif
    m3dc1_itime = m3dc1_time           ! Variables in module
    m3dc1_factors = m3dc1_scale_factors
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 0
    m3dc1_phases_deg = m3dc1_phase_shift_deg
    Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
    If (setup_bfield_verbose) Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    If (setup_bfield_verbose) Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
    If (m3dc1_toroidal_on_err .AND. setup_bfield_verbose) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  End Subroutine setup_bfield_m3dc1_full
  
  ! *********************************************
  ! ***************  M3DC1 AS  ******************
  ! *********************************************
  Subroutine setup_bfield_m3dc1_as
    Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factors, &
         m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
         m3dc1_field_type, m3dc1_phases_deg
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS m3dc1_as'
    bfield%method      = 5
    bfield%method_2d   = 5
    bfield%method_pert = -1
    ! Setup field    
    If (m3dc1_time .eq. -1) Then
      Write(*,*) 'Error: Bfield type of M3DC1 is set but m3dc1_time is not. Exiting.'
      Stop
    Endif
    m3dc1_itime = m3dc1_time           ! Variables in module
    m3dc1_factors = m3dc1_scale_factors
    m3dc1_toroidal_on_err = m3dc1_toroidal_off_grid
    m3dc1_field_type = 0
    m3dc1_phases_deg = m3dc1_phase_shift_deg
    Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
    If (setup_bfield_verbose) Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
    If (setup_bfield_verbose) Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
    If (m3dc1_toroidal_on_err .AND. setup_bfield_verbose) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
  End Subroutine setup_bfield_m3dc1_as
#endif
  ! *********************************************
  ! *************** IPEC ************************
  ! *********************************************    
  Subroutine setup_bfield_ipec
    Use g3d_module, Only : readg_g3d
    Use ipec_module, Only: open_ipec_fields
    Implicit None
    
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS IPEC'
    If (setup_bfield_verbose) Write(*,'(a,i0)') '-----> ipec_field_eval_type is: ',ipec_field_eval_type
    If (ipec_field_eval_type .eq. 0) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating IPEC fields as EQUILIBRIUM ONLY!'
      bfield%method = 7
      bfield%method_2d = 0
      bfield%method_pert = -1
    Elseif (ipec_field_eval_type .eq. 1) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating IPEC fields as EQ + VACUUM!'
      bfield%method = 8
      bfield%method_2d = 0
      bfield%method_pert = -1
    Elseif (ipec_field_eval_type .eq. 2) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating IPEC fields as EQ + PERT!'
      bfield%method = 9
      bfield%method_2d = 0
      bfield%method_pert = -1
    Elseif (ipec_field_eval_type .eq. 3) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating IPEC fields as VAC PERT ONLY!'
      bfield%method = 12
      bfield%method_2d = 0
      bfield%method_pert = 12
    Elseif (ipec_field_eval_type .eq. 4) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating IPEC fields as RESPONSE PERT ONLY!'
      bfield%method = 13
      bfield%method_2d = 0
      bfield%method_pert = 13
    Else
      Stop "Did not recognize ipec_field_eval_type"
    Endif
    Call readg_g3d(gfile_name,g)
    bfield%g = g
    Call open_ipec_fields(ipec_run_path,ipec_itype)
  End Subroutine setup_bfield_ipec

  ! *********************************************
  ! *************** XPAND ***********************
  ! *********************************************    
  Subroutine setup_bfield_xpand
    Use g3d_module, Only : readg_g3d
    Use xpand_module, Only: open_xpand_fields
    Implicit None
    
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS XPAND'
    If (setup_bfield_verbose) Write(*,'(a,i0)') '-----> xpand_field_eval_type is: ',xpand_field_eval_type
    If (xpand_field_eval_type .eq. 0) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating XPAND fields as PERTURBED!'
      bfield%method      = 10
      bfield%method_2d   = 0
      bfield%method_pert = -1
    Elseif (xpand_field_eval_type .eq. 1) Then
      If (setup_bfield_verbose) Write(*,'(a)') '-----> Evaluating XPAND fields as VACUUM!'
      bfield%method      = 11
      bfield%method_2d   = 0
      bfield%method_pert = -1
    Else
      Stop "Did not recognize xpand_field_eval_type"
    Endif
    Call readg_g3d(gfile_name,g)
    bfield%g = g
    Call open_xpand_fields(xpand_fname)    
  End Subroutine setup_bfield_xpand

#ifdef HAVE_FXDR    
  ! *********************************************
  ! ***************  XDR  ***********************
  ! *********************************************    
  Subroutine setup_bfield_xdr
    Use xdr_routines_mod, Only : readbgrid_xdr
    Implicit None
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS XDR'
    Call readbgrid_xdr(xdr_fname,xdr_check,xdr_verbose)
    bfield%method      = 15
    bfield%method_2d   = -1
    bfield%method_pert = -1
  End Subroutine setup_bfield_xdr
#endif
  
  ! *********************************************
  ! *************** BGRID ***********************
  ! *********************************************    
  Subroutine setup_bfield_bgrid
    Use bgrid_module, Only: open_bgrid_fields
    Implicit None
    
    If (setup_bfield_verbose) Write(*,'(a)') '-----> BFIELD METHOD IS BGRID'
      bfield%method      = 16
      bfield%method_2d   = -1
      bfield%method_pert = -1
    Call open_bgrid_fields(bgrid_fname)    
  End Subroutine setup_bfield_bgrid

  
End Module setup_bfield_module
  

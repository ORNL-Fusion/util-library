Module bnorm_bfield_mod
  Use kind_mod, Only: int32, real64
  Use bfield, Only : bfield_type, coil_type, g_type

  Implicit None

  Integer(int32), parameter :: max_m3dc1_files = 10
  Logical, Parameter :: m3dc1_toroidal_off_grid = .true.
  Integer(int32), parameter :: max_rmp_coils       = 12

  Real(real64) :: &
       rmp_current(max_rmp_coils) = 0.d0, &
       m3dc1_scale_factors(max_m3dc1_files) = 0.d0
  
  Character(Len=120) :: &
       rmp_type                         = 'none',  &
       gfile_name                       = 'none', &
       rmp_coil_type                    = 'none', &
       m3dc1_filenames(max_m3dc1_files) = 'none', &
       ipec_run_path                    = 'none', &
       xpand_fname                      = 'none'

  Integer(int32) :: &
       m3dc1_time = -1, &
       m3dc1_nsets = 0, &
       ipec_field_eval_type = -1, &
       xpand_field_eval_type = -1, &
       ipec_itype = 0

  Type(bfield_type) :: bfield
  Type(g_type) :: g
  Type(coil_type) :: coil

  
  Namelist / bfield_nml / &
       rmp_type, &
       gfile_name, &
       rmp_coil_type, rmp_current, &
       m3dc1_filenames, m3dc1_scale_factors, m3dc1_time, m3dc1_nsets, &
       ipec_run_path, ipec_field_eval_type, ipec_itype, &
       xpand_fname, xpand_field_eval_type
  
  
Contains

  Subroutine setup_bfield_bnorm
    Use ipec_module, Only: open_ipec_fields
    Use xpand_module, Only: open_xpand_fields
    Use rmpcoil_module, Only : build_d3d_ccoils_jl, build_d3d_icoils_jl
    Use M3DC1_routines_mod, Only : prepare_m3dc1_fields, m3dc1_factors, m3dc1_itime, m3dc1_toroidal_on_err,bfield_m3dc1, &
         m3dc1_field_type
    Use g3d_module, Only : readg_g3d
    
    Implicit None
    Integer(int32) :: iocheck
    
    Open(99,file="bnorm_settings.nml",status="old",form="formatted",iostat=iocheck)
    If ( iocheck /= 0 ) Then
      Write(*,*) 'Error opening namelist file'
      Stop 'Exiting: I/O Error in bnorm_driver.f90'
    Endif
    Read(99,nml=bfield_nml)
    Close(99)


    
    ! Setup rmp field
    Select Case (rmp_type)
    Case ('g3d+rmpcoil')
      Write(*,'(a)') '-----> BFIELD METHOD IS G3D+RMPCOIL'
      bfield%method = 1
      bfield%method_2d = 0
      bfield%method_pert = 6
      Call readg_g3d(gfile_name,g)
      bfield%g = g
      ! load rmp coil    
      Select Case (rmp_coil_type)
      Case ('d3d_ccoils')
        Write(*,'(a)') '-----> Rmp coils are DIII-D C coils'
        Write(*,'(a,6f12.3)') ' Coil currents: ', rmp_current(1:6)
        Call build_d3d_ccoils_jl(coil,rmp_current(1:6))
      Case ('d3d_icoils')
        Write(*,'(a)') '----->  Rmp coils are DIII-D I coils'
        Write(*,'(a,12f12.3)') ' Coil currents 1: ', rmp_current(1:6)
        Write(*,'(a,12f12.3)') ' Coil currents 2: ', rmp_current(7:12)
        Call build_d3d_icoils_jl(coil,rmp_current(1:12))
      Case Default
        Write(*,*) 'Error!  Unknown rmp_coil_type'
        Write(*,*) 'Current options are:'
        Write(*,*) '''d3d_ccoils'''        
        Write(*,*) '''d3d_icoils'''        
        Stop
      End Select
      bfield%coil = coil
    Case ('g3d+m3dc1')
      Write(*,'(a)') '-----> BFIELD METHOD IS G3D+M3DC1'
      bfield%method    = 3
      bfield%method_2d = 0
      bfield%method_pert = 4
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
      Call prepare_m3dc1_fields(m3dc1_filenames(1:m3dc1_nsets))
      Write(*,'(a,i0)') '---------> M3DC1 time (0 vacuum, 1 response): ',m3dc1_itime
      Write(*,*) '---------> M3DC1 scale factor: ',m3dc1_factors(1:m3dc1_nsets)
      If (m3dc1_toroidal_on_err) Write(*,'(a)') '---------> M3DC1 fields will be set to B=Bt=1 off grid!'
    Case ('ipec')
      Write(*,'(a)') '-----> BFIELD METHOD IS IPEC'
      Write(*,'(a,i0)') '-----> ipec_field_eval_type is: ',ipec_field_eval_type
      If (ipec_field_eval_type .eq. 3) Then
        Write(*,'(a)') '-----> Evaluating IPEC fields as VAC PERT ONLY!'
        bfield%method = 12
        bfield%method_2d = 0
        bfield%method_pert = 12
      Elseif (ipec_field_eval_type .eq. 4) Then
        Write(*,'(a)') '-----> Evaluating IPEC fields as RESPONSE PERT ONLY!'
        bfield%method = 13
        bfield%method_2d = 0
        bfield%method_pert = 13
      Else
        Stop "Only ipec_field_eval_type == (2 or 3) is valid for bnorm"
      Endif
      Call readg_g3d(gfile_name,g)
      bfield%g = g
      Call open_ipec_fields(ipec_run_path,ipec_itype)
    Case Default
      Write(*,*) 'Unknown rmp_type in bnorm_driver:',rmp_type
      Write(*,*) 'Current options are:'
      Write(*,*) '''g3d+rmpcoil'''
      Write(*,*) '''g3d+m3dc1'''
      Write(*,*) '''ipec'''
      Stop      
    End Select



  End Subroutine setup_bfield_bnorm

End Module bnorm_bfield_mod

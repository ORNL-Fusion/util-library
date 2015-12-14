Program test_int_two_lines
  Use kind_mod, Only: int32, real64
  Use math_geo_module, Only: int_two_lines
  Implicit None
  
  Real(real64), Dimension(2) :: p1, p2, p3, p4, pint
  Real(real64), Dimension(2) :: pans
  Real(real64) :: u1, u2, diff
  Integer(int32) :: itest, ipass, ierr_ans, ierr
  Real(real64), Parameter :: dtol = 2.e-16_real64
  Integer(int32), Parameter :: ntest = 6

  Logical, Parameter :: DEBUG = .false.
  
  ipass = 0
  Do itest = 1,ntest

    Select Case (itest)
    Case (1)
      p1 = [0.5_real64, 0.5_real64]      ! INTERSECTION WITHIN SEGMENTS
      p2 = [1.0_real64, 1.0_real64]
      p3 = [0.5_real64, 1.0_real64]
      p4 = [1.0_real64, 0.5_real64]

      pans = [0.75_real64, 0.75_real64]
      ierr_ans = 0
      
    Case (2)
      p1 = [0.5_real64, 0.5_real64]      ! INTERSECTION BEYOND SEGMENTS
      p2 = [0.55_real64, 0.55_real64]
      p3 = [0.5_real64, 1.0_real64]
      p4 = [0.55_real64, 0.95_real64]

      pans = [0.75_real64, 0.75_real64]
      ierr_ans = 0

    Case (3)
      p1 = [0.5_real64, 0.5_real64]      ! COLLINEAR LINES
      p2 = [1.0_real64, 1.0_real64]
      p3 = [0.25_real64, 0.25_real64]
      p4 = [1.25_real64, 1.25_real64]

      pans = [0._real64, 0._real64]
      ierr_ans = 1

    Case (4)
      p1 = [0.5_real64, 0.5_real64]      ! PARALLEL LINES
      p2 = [1.0_real64, 1.0_real64]
      p3 = [0.5_real64, 1.0_real64]
      p4 = [1.0_real64, 1.5_real64]

      pans = [0._real64, 0._real64]
      ierr_ans = 1

    Case (5)
      p1 = [0.5_real64, 0.5_real64]      ! INTERSECTION WITH LINE ENDPOINTS EQUAL
      p2 = [1.0_real64, 1.0_real64]
      p3 = [0.5_real64, 1.0_real64]
      p4 = [1.0_real64, 1.0_real64]

      pans = [1._real64, 1._real64]
      ierr_ans = 0

    Case (6)
      p1 = [0.5_real64, 0.5_real64]      ! INTERSECTION WITH INT ON LINE
      p2 = [1.0_real64, 1.0_real64]
      p3 = [0.5_real64, 1.0_real64]
      p4 = [0.75_real64, 0.75_real64]

      pans = [0.75_real64, 0.75_real64]
      ierr_ans = 0
            
    Case Default
      Stop "Should not be here"

    End Select
    Write(*,'(a,i0,a,i0,a)') 'Running test (',itest,'/',ntest,')'
    If (debug) Then
      Write(*,*) 'p1',p1
      Write(*,*) 'p2',p2
      Write(*,*) 'p3',p3
      Write(*,*) 'p4',p4
    Endif
    
    Call int_two_lines(p1,p2,p3,p4,u1,u2,ierr,pint)

    If (DEBUG) Then
      Write(*,*) 'P_int',pint
      Write(*,*) 'P_ans',pans
      Write(*,*) 'ierr    ',ierr
      Write(*,*) 'ierr_ans',ierr_ans
    Endif
    
    diff = Sqrt(Sum(pans-pint)**2)
    If ((ierr == ierr_ans) .AND. (diff .le. dtol)) Then
      Write(*,*) 'PASS'
      ipass = ipass + 1
    Else
      Write(*,*) 'FAIL'
    Endif
    
    
  Enddo

  Write(*,'(a)') '********************************************'
  Write(*,'(i0,a,i0,a)') ipass,' of ',ntest,' tests passed'
  Write(*,'(a)') '********************************************'
  
End Program test_int_two_lines

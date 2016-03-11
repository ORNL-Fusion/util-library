Program test_int_curve_curve
  Use kind_mod, Only: int32, real64
  Use math_geo_module, Only: int_curve_curve
  Implicit None
  
  Real(real64), Allocatable, Dimension(:) :: rline1, zline1, rline2, zline2
  Real(real64), Dimension(2) :: pans, pint
  Real(real64) :: diff
  Integer(int32) :: itest, ipass, ierr_ans, ierr, nline1, nline2, found_ind1_ans, found_ind2_ans
  Integer(int32) :: int_count_ans, found_ind1, found_ind2, int_count
  Logical :: first
  Logical, Parameter :: DEBUG = .true.

  Real(real64), Parameter :: dtol = 2.e-16_real64
  Integer(int32), Parameter :: ntest = 6

  
  ipass = 0
  Do itest = 1,ntest

    Select Case (itest)
    Case (1)
      nline1 = 4
      nline2 = 4
      Allocate(rline1(nline1), zline1(nline1))
      Allocate(rline2(nline2), zline2(nline2))
      
      rline1 = [0._real64, 1._real64, 2._real64, 3._real64]
      zline1 = [0._real64, 1._real64, 4._real64, 3._real64]
      rline2 = [0._real64, 1._real64, 2._real64, 3._real64]
      zline2 = [1._real64,-1._real64, 3._real64, 5._real64]

      first = .true.
      
      pans(1) = 1._real64/3._real64
      pans(2) = 1._real64/3._real64
      ierr_ans = 0
      found_ind1_ans = 1
      found_ind2_ans = 1
      int_count_ans = 1
      
    Case (2)
      nline1 = 2
      nline2 = 2
      Allocate(rline1(nline1), zline1(nline1))
      Allocate(rline2(nline2), zline2(nline2))
      
      rline1 = [0._real64, 1._real64]
      zline1 = [0._real64, 1._real64]
      rline2 = [0._real64, 1._real64]
      zline2 = [1._real64,-1._real64]

      first = .true.
      
      pans(1) = 1._real64/3._real64
      pans(2) = 1._real64/3._real64
      ierr_ans = 0
      found_ind1_ans = 1
      found_ind2_ans = 1
      int_count_ans = 1

    Case (3)
      nline1 = 4
      nline2 = 4
      Allocate(rline1(nline1), zline1(nline1))
      Allocate(rline2(nline2), zline2(nline2))

      rline1 = [ 0._real64, 1._real64, 2._real64, 3._real64]
      zline1 = [ 0._real64, 1._real64, 4._real64, 3._real64]
      rline2 = [-2._real64,-1._real64, 0._real64, 1._real64]
      zline2 = [ 1._real64,-1._real64, 1._real64,-1._real64]

      first = .true.
      
      pans(1) = 1._real64/3._real64
      pans(2) = 1._real64/3._real64
      ierr_ans = 0
      found_ind1_ans = 1
      found_ind2_ans = 3
      int_count_ans = 1

    Case (4)
      nline1 = 4
      nline2 = 4
      Allocate(rline1(nline1), zline1(nline1))
      Allocate(rline2(nline2), zline2(nline2))

      rline1 = [ 0._real64, 1._real64, 2._real64, 3._real64]
      zline1 = [ 1._real64, 1._real64, 4._real64, 3._real64]
      rline2 = [-2._real64,-1._real64, 0._real64, 1._real64]
      zline2 = [ 1._real64,-1._real64, 1._real64,-1._real64]

      first = .true.
      
      pans(1) = 0._real64
      pans(2) = 1._real64
      ierr_ans = 0
      found_ind1_ans = 1
      found_ind2_ans = 3
      int_count_ans = 1

    Case (5)
      nline1 = 4
      nline2 = 4
      Allocate(rline1(nline1), zline1(nline1))
      Allocate(rline2(nline2), zline2(nline2))

      rline1 = [-3._real64, 1._real64, 2._real64, 3._real64]
      zline1 = [ 1._real64, 1._real64, 4._real64, 3._real64]
      rline2 = [-2._real64,-1._real64, 0._real64, 1._real64]
      zline2 = [ 1._real64,-1._real64, 1._real64,-1._real64]

      first = .true.
      
      pans(1) =-2._real64
      pans(2) = 1._real64
      ierr_ans = 0
      found_ind1_ans = 1
      found_ind2_ans = 1
      int_count_ans = 1
      
    Case (6)
      nline1 = 4
      nline2 = 4
      Allocate(rline1(nline1), zline1(nline1))
      Allocate(rline2(nline2), zline2(nline2))

      rline1 = [-3._real64, 1._real64, 2._real64, 3._real64]
      zline1 = [ 2._real64, 2._real64, 4._real64, 3._real64]
      rline2 = [-2._real64,-1._real64, 0._real64, 1._real64]
      zline2 = [ 1._real64, 1._real64, 1._real64,-1._real64]

      first = .true.
      
      pans(1) = 0._real64
      pans(2) = 0._real64
      ierr_ans = 1
      found_ind1_ans = 0
      found_ind2_ans = 0
      int_count_ans = 0

    Case Default
      Stop "Should not be here"

    End Select
    Write(*,'(a,i0,a,i0,a)') 'Running test (',itest,'/',ntest,')'
    If (debug) Then
      Write(*,*) 'rline1',rline1
      Write(*,*) 'zline1',zline1
      Write(*,*) 'rline2',rline2
      Write(*,*) 'zline2',zline2
    Endif
    
    Call int_curve_curve(rline1,zline1,nline1,rline2,zline2,nline2,first,pint,ierr,found_ind1,found_ind2,int_count)

    If (DEBUG) Then
      Write(*,*) 'P_int',pint
      Write(*,*) 'P_ans',pans
      Write(*,*) 'ierr    ',ierr
      Write(*,*) 'ierr_ans',ierr_ans
      Write(*,*) 'found_ind1, ans',found_ind1, found_ind1_ans
      Write(*,*) 'found_ind2, ans',found_ind2, found_ind2_ans
    Endif
    
    diff = Sqrt(Sum(pans-pint)**2)
    If ((ierr == ierr_ans) .AND. (diff .le. dtol) .AND. (found_ind1 == found_ind1_ans) &
         .AND. (found_ind2 == found_ind2_ans) .AND. (int_count == int_count_ans)) Then
      Write(*,*) 'PASS'
      ipass = ipass + 1
    Else
      if (DEBUG) Then
        Write(*,*) 'ierr equal?',(ierr == ierr_ans)
        Write(*,*) 'found_ind1 equal?',(found_ind1 == found_ind1_ans)
        Write(*,*) 'found_ind2 equal?',(found_ind2 == found_ind2_ans)
        Write(*,*) 'int_count equal?',(int_count == int_count_ans)
        Write(*,*) 'diff ok?',diff,dtol,(diff .le. dtol)
      Endif
      Write(*,*) 'FAIL'
    Endif

    Deallocate(rline1,zline1,rline2,zline2)
    
  Enddo

  Write(*,'(a)') '********************************************'
  Write(*,'(i0,a,i0,a)') ipass,' of ',ntest,' tests passed'
  Write(*,'(a)') '********************************************'
  
End Program test_int_curve_curve

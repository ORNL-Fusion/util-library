!-----------------------------------------------------------------------------
!+ Kind specifications
!-----------------------------------------------------------------------------
Module kind_mod
!
! Description:
!   This module contains the kind specifications for all subroutines 
!   and modules.
!
! History:
! Version   Date      Comment
! -------   ----      -------
!  1.0     04/12/2011   Adapted from PENTA.  JL
! 
! Author(s): J. Lore 7/2009 - 4/12/2011
!
  Implicit none
  Integer, parameter :: rknd = selected_real_kind(15,307) 
  Integer, parameter :: iknd = selected_int_kind(8)       
End module kind_mod
!- End of header -------------------------------------------------------------


Module bfield_typedef
  Use kind_mod, Only : int32, real64
  Use g_typedef, Only : g_type
  Use coil_typedef, Only : coil_type
  Implicit None
  Private
  Type, Public :: bfield_type
    Integer(int32) :: method = 0
    Type(g_type) :: g
    Type(coil_type) :: coil    
  End Type bfield_type
End Module bfield_typedef

!   Set bfield_method to control the fieldline deriviative calls.
!   Appropriate loading must be done before calls to any fieldline following routine 
!   (e.g., gfile loading, rmp coil generation)
!
!    bfield%method == 
!                     0 -- gfile field only
!                     1 -- gfile + rmp coil field
!                     2 -- gfile + Pavel's screened fields with bspline interpolation (not fully implemented!)
!                     3 -- gfile + M3DC1 perturbed field
!                     4 -- M3DC1 total field
!                     5 -- M3DC1 2D field only
!                     6 -- Just coils
!                     7 -- ipec eq only
!                     8 -- ipec vacuum
!                     9 -- ipec perturbed
!                    10 -- xpand perturbed
!                    11 -- xpand vacuum
! 


Module bfield
  Use bfield_typedef, Only : bfield_type
  Implicit None
  Private
  Public :: bfield_type
  
End Module bfield


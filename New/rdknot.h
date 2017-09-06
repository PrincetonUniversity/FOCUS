
!title (axis) ! Reads axis parameters from file.

!latex \briefly{The parameters describing a closed curve are read from file.}

!latex \calledby{\link{focus}}
!latex \calls{\link{knotxx}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] A user-prescribed closed curve in three-dimensional space will serve as the proxy magnetic axis.
!latex \item[2.] The parameters describing the knot, namely \internal{axisxc}, \internal{axisyss} and \internal{axiszs}, are read from \verb+ext.axis+.
!latex \item[3.] A circular or elliptical cross section surface is constructed.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

subroutine rdknot
  
  use globals, only : zero, one, half, ten, pi2, myid, ncpu, ounit, runit, &
                      ext, &
                      axisNF, axisxc, axisxs, axisyc, axisys, axiszc, axiszs, &
                      IsSymmetric, Nfp, Npc, discretefactor, &
                      knotsurf, Nteta, Nzeta, surf
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  LOGICAL              :: exist
  INTEGER              :: iostat, astat, ierr, ii, jj
  REAL                 :: teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3), ds(1:3), dd
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  inquire( file=trim(ext)//".axis", exist=exist )
  if( exist ) then
   if( myid.eq.0 ) then
    open(runit, file=trim(ext)//".axis", status="old", action='read', iostat=iostat)  
    write(ounit,'("rdknot  : reading ext.axis ;")') 
   endif
  else
   inquire( file=".axis", exist=exist )
   if( exist ) then
    if( myid.eq.0 ) then
     open(runit, file=".axis", status="old", action='read', iostat=iostat)  
     write(ounit,'("rdcoil  : reading .axis ;")') 
    endif
   endif
  endif
  
  FATAL( rdknot, .not.exist, neither ext.axis nor .axis found )
  
  if( myid.eq.0 ) then
   read( runit, *, iostat=iostat )
   read( runit, *, iostat=iostat ) axisNF
  endif
  
  IlBCAST( axisNF, 1, 0 )
  
  SALLOCATE( axisxc, (0:axisNF), zero )
  SALLOCATE( axisxs, (0:axisNF), zero )
  SALLOCATE( axisyc, (0:axisNF), zero )
  SALLOCATE( axisys, (0:axisNF), zero )
  SALLOCATE( axiszc, (0:axisNF), zero )
  SALLOCATE( axiszs, (0:axisNF), zero )
  
  if( myid.eq.0 ) then
   
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisxc(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisxs(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisyc(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axisys(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axiszc(0:axisNF)
   read( runit, *, iostat=iostat ) 
   read( runit, *, iostat=iostat ) axiszs(0:axisNF)
   
   close(runit,iostat=iostat)
   
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
  RlBCAST( axisxc(0:axisNF), axisNF+1, 0 )
  RlBCAST( axisxs(0:axisNF), axisNF+1, 0 )
  RlBCAST( axisyc(0:axisNF), axisNF+1, 0 )
  RlBCAST( axisys(0:axisNF), axisNF+1, 0 )
  RlBCAST( axiszc(0:axisNF), axisNF+1, 0 )
  RlBCAST( axiszs(0:axisNF), axisNF+1, 0 )
  
  if( myid.eq.0 ) then
   write(ounit,'("rdknot  : axisxc =",    999es11.03)') axisxc(0:axisNF)
   write(ounit,'("rdknot  : axisxs =",11x,998es11.03)') axisxs(1:axisNF)
   write(ounit,'("rdknot  : axisyc =",    999es11.03)') axisyc(0:axisNF)
   write(ounit,'("rdknot  : axisys =",11x,999es11.03)') axisys(1:axisNF)
   write(ounit,'("rdknot  : axiszc =",    999es11.03)') axiszc(0:axisNF)
   write(ounit,'("rdknot  : axiszs =",11x,999es11.03)') axiszs(1:axisNF)
  endif
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  
! discretize the surface data; 2017/03/28; czhu;
  
! select case (IsSymmetric)
! case ( 0 )
!    Nfp = 1                          !reset Nfp to 1;
!    Npc = 1                          !number of coils periodicity 
! case ( 1 )                          !plasma periodicity enabled;
!    Npc = 1
! case ( 2 )                          !plasma and coil periodicity enabled;
!    Npc = Nfp
! end select
! discretefactor = discretefactor/Nfp

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!
  return
  
end subroutine rdknot

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!!


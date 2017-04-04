
!title (winding) ! The winding surface is read from file.

!latex \briefly{A Fourier representation for the winding surface is read from file.}

!latex \calledby{\link{xoptim}}
!l tex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] The Fourier harmonics of the winding surface are required in \verb+winding.surface+.
!latex The format of this file is as follows: \begin{verbatim}
!latex wmn wNfp   ! integers
!latex wim(1:wmn) ! integers: poloidal mode identification;
!latex win(1:wmn) ! integers: toroidal mode identification;
!latex Rwc(1:wmn) ! real    : cylindrical R cosine harmonics;
!latex Rws(1:wmn) ! real    : cylindrical R   sine harmonics;
!latex Zwc(1:wmn) ! real    : cylindrical Z cosine harmonics;
!latex Zws(1:wmn) ! real    : cylindrical Z   sine harmonics; \end{verbatim}
!latex \item[2.] Note that immediately after reading (and broadcasting) \verb+win+, the field periodicity factor is included, i.e. \verb+win = win * wNfp+.
!latex \item[3.] The winding surface, ${\bf x}_W(\t,\z)$, is constructed in \link{windsf}.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine rwsurf
  
  use knotopt, only : zero, myid, ncpu, ounit, lunit, ext, &
                    wmn, wNfp, wim, win, Rwc, Rws, Zwc, Zws
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  LOGICAL :: exist
  INTEGER :: iosta, iostb, iostc, iostd, ioste, iostf, astat, ierr
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  inquire( file="winding.surface", exist=exist)
  
  FATAL(rwsurf, .not.exist, winding.surface does not exist)
  
  if( myid.eq.0 ) open(lunit, file="winding.surface", status='old', action='read', iostat=iosta)
  
  IlBCAST(iosta,1,0)
  
  FATAL(rwsurf, iosta.ne.0, error opening winding.surface)
  
  if( myid.eq.0 ) read(lunit,*,iostat=iosta) wmn, wNfp
  
  IlBCAST(iosta,1,0)
  
  FATAL(rwsurf, iosta.ne.0, error reading wmn wNfp from winding.surface)
  
  IlBCAST(wmn ,1,0)
  IlBCAST(wNfp,1,0)
  
  FATAL(rwsurf, wmn.le.0, invalid)
  
  SALLOCATE(wim,(1:wmn),0)
  SALLOCATE(win,(1:wmn),0)
  
  SALLOCATE(Rwc,(1:wmn),zero)
  SALLOCATE(Rws,(1:wmn),zero)
  SALLOCATE(Zwc,(1:wmn),zero)
  SALLOCATE(Zws,(1:wmn),zero)
  
  if( myid.eq.0 ) read(lunit,*,iostat=iosta) wim(1:wmn)
  if( myid.eq.0 ) read(lunit,*,iostat=iostb) win(1:wmn)
  if( myid.eq.0 ) read(lunit,*,iostat=iostc) Rwc(1:wmn)
  if( myid.eq.0 ) read(lunit,*,iostat=iostd) Rws(1:wmn)
  if( myid.eq.0 ) read(lunit,*,iostat=ioste) Zwc(1:wmn)
  if( myid.eq.0 ) read(lunit,*,iostat=iostf) Zws(1:wmn)
 
  IlBCAST(iosta,1,0)
  IlBCAST(iostb,1,0)
  IlBCAST(iostc,1,0)
  IlBCAST(iostd,1,0)
  IlBCAST(ioste,1,0)
  IlBCAST(iostf,1,0)

  FATAL(rwsurf, iosta.ne.0, error reading wim)
  FATAL(rwsurf, iostb.ne.0, error reading win)
  FATAL(rwsurf, iostc.ne.0, error reading Rwc)
  FATAL(rwsurf, iostd.ne.0, error reading Rws)
  FATAL(rwsurf, ioste.ne.0, error reading Zwc)
  FATAL(rwsurf, iostf.ne.0, error reading Zws)

  IlBCAST(wim(1:wmn),wmn,0)
  IlBCAST(win(1:wmn),wmn,0)
  
  win(1:wmn) = win(1:wmn) * wNfp
  
  RlBCAST(Rwc(1:wmn),wmn,0)
  RlBCAST(Rws(1:wmn),wmn,0)
  RlBCAST(Zwc(1:wmn),wmn,0)
  RlBCAST(Zws(1:wmn),wmn,0)

  if( myid.eq.0 ) close(lunit,iostat=iosta)
  
  IlBCAST(iosta,1,0)
  
  FATAL(rwsurf, iosta.ne.0, error closing winding.surface)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( myid.eq.0 ) then
   write(ounit,'("rwsurf : " 10x " : wim ="99i13   )') wim(1:wmn)
   write(ounit,'("rwsurf : " 10x " : win ="99i13   )') win(1:wmn)
   write(ounit,'("rwsurf : " 10x " : Rwc ="99es13.5)') Rwc(1:wmn)
   write(ounit,'("rwsurf : " 10x " : Zws ="99es13.5)') Zws(1:wmn)
  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  return
  
end subroutine rwsurf

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

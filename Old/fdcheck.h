!title (fd check) ! Read input file, initialize global variables.

!latex \briefly{Check the derivatives using finite difference method}

!latex \calledby{\link{knotopt}}
!latex \calls{\link{}}

!latex \tableofcontents

!latex \subsection{overview}
!latex \bi
!latex \item[1.] 
!latex \ei
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine testderivs1

  use kmodule, only: zero, half, machprec, coil, Ncoils, Cdof, NFcoil, Ndof, t1E, totalenergy, ncpu, myid, ounit, mincc, &
                     weight_eqarc, n1E, Inorm, Gnorm
                     
  implicit none
  include "mpif.h"
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER         :: icoil, ifd, idof, iteta, jzeta, jcoil, astat, ierr, infou
  REAL            :: tmpE(-1:1), diff, rdiff, xdof(Ndof), bdof(1:Ndof), start, finish, Ax(0:Cdof), Ay(0:Cdof), Az(0:Cdof), lx, ly, lz
  REAL, parameter :: fd = 1.0E-6

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!  
  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : Checking the first derivatives using finite-difference method")')

!!$!----------------------------------------------
!!$
!!$  ALLOCATE(mincc(1:Ncoils, 1:Ncoils))
!!$  icoil = 1; iteta = 1; jzeta = 1; jcoil = 2
!!$  call pack( xdof )  !back up the current coils data
!!$  call epotent1(icoil, jcoil, Ax(0:Cdof))
!!$  
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : coil# : dof : dof# :   numerical value"6X" :   fd-method value"6X" :   difference"11X" :   relative diff")')
!!$
!!$        
!!$!-------X_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%xc(idof) = coil(icoil)%xc(idof) + half * fd * ifd; call discretecoil
!!$           call epotent0(icoil,jcoil,lx); tmpE(ifd)            = lx
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = Ax(idof+1); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'XC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$     enddo
!!$
!!$!-------X_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%xs(idof) = coil(icoil)%xs(idof) + half * fd * ifd; call discretecoil
!!$           call epotent0(icoil,jcoil,lx); tmpE(ifd)            = lx
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = Ax(idof+NFcoil+2); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'XS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Y_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%yc(idof) = coil(icoil)%yc(idof) + half * fd * ifd; call discretecoil
!!$           call epotent0(icoil,jcoil,lx); tmpE(ifd)            = lx
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = Ax(idof+2*NFcoil+3); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'YC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Y_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%ys(idof) = coil(icoil)%ys(idof) + half * fd * ifd; call discretecoil
!!$           call epotent0(icoil,jcoil,lx); tmpE(ifd)            = lx
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = Ax(idof+3*NFcoil+4); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'YS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Z_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%zc(idof) = coil(icoil)%zc(idof) + half * fd * ifd; call discretecoil
!!$           call epotent0(icoil,jcoil,lx); tmpE(ifd)            = lx
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = Ax(idof+4*NFcoil+5); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'ZC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$!-------Z_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%zs(idof) = coil(icoil)%zs(idof) + half * fd * ifd; call discretecoil
!!$           call epotent0(icoil,jcoil,lx); tmpE(ifd)            = lx
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = Ax(idof+5*NFcoil+6); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'ZS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$ 
!!$  call unpack( xdof )  !recover data
!!$ 
!!$  return






  call cpu_time(start)
  call costfun(1); 
  call cpu_time(finish)
  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : First order derivatives of energy function takes " ES23.15 " seconds.")') finish - start

  call pack( bdof )  !back up the current coils data
  
  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : coil# : Fourier# :      numerical value"3X" :   fd-method value"6X" :   difference"11X" :   relative diff")')

  call cpu_time(start)
  do idof = 1, Ndof

     do ifd = -1, 1, 2
        xdof = bdof
        xdof(idof) = xdof(idof) + half * fd * ifd
        call unpack ( xdof )
        call costfun(0)
        tmpE(ifd) = totalenergy
     enddo ! end do ifd

     call DoFconvert(idof, icoil, infou)
     tmpE(0) = t1E(icoil,infou) ; diff = tmpE(0)- (tmpE(1) - tmpE(-1)) / fd

     if( abs(tmpE(0)) .lt. machprec ) then
         rdiff = 0
     else
         rdiff = diff / tmpE(0)
     endif
     
     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "I8, 4(" : "ES23.15))') icoil, infou, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
  enddo
  call cpu_time(finish)
  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : First order derivatives of energy function takes " ES23.15 " seconds.")') finish - start
!!$ 
!!$  !-----------------test Lagrange multipiler-------------------
!!$  if( weight_eqarc .gt. zero ) then
!!$
!!$     call ntchdms(1)
!!$     do icoil = 1, Ncoils
!!$
!!$        do idof = 0, NFcoil
!!$           do ifd = -1, 1, 2
!!$              coil(icoil)%lmdc(idof) = 1.0; coil(icoil)%lmdc(idof) = coil(icoil)%lmdc(idof) + half * fd * ifd
!!$              call costfun(0)     ; tmpE(ifd)            = totalenergy
!!$           enddo ! end do ifd
!!$           
!!$           tmpE(0) = n1E(icoil, Cdof+1+idof); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$           if( abs(tmpE(0)) .lt. machprec ) then
!!$              rdiff = 0
!!$           else
!!$              rdiff = diff / tmpE(0)
!!$           endif
!!$
!!$           if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                icoil, 'LC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$        enddo
!!$
!!$        do idof = 0, NFcoil
!!$           do ifd = -1, 1, 2
!!$              coil(icoil)%lmds(idof) = 1.0; coil(icoil)%lmds(idof) = coil(icoil)%lmds(idof) + half * fd * ifd
!!$              call costfun(0)     ; tmpE(ifd)            = totalenergy
!!$           enddo ! end do ifd
!!$
!!$           tmpE(0) = n1E(icoil, Cdof+NFcoil+2+idof); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$           if( abs(tmpE(0)) .lt. machprec ) then
!!$              rdiff = 0
!!$           else
!!$              rdiff = diff / tmpE(0)
!!$           endif
!!$
!!$           if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                icoil, 'LS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$        enddo
!!$     enddo
!!$
!!$  endif

  call unpack( bdof )  !recover data
 
  return
  
end subroutine testderivs1

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine testderivs2

  use kmodule, only: half, machprec, coil, Ncoils, NFcoil, Ndof, t1E, t2E, totalenergy, ncpu, myid, ounit, ccsep, t1C, t2C, Cdof, &
                     mincc, n2E, n1E,  weight_eqarc, Loptimize
  implicit none
  include "mpif.h"
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER         :: icoil, ifd, idof, jcoil, jdof, iteta, jzeta, astat, ierr, infou
  REAL            :: tmpE(-1:4), diff, rdiff, xdof(Ndof), bdof(Ndof), start, finish, lbx(0:Cdof,0:Cdof),lby(0:Cdof,0:Cdof),lbz(0:Cdof,0:Cdof), &
                     lbnorm, l1B(1:Ncoils,0:Cdof), l2B(1:Ncoils,0:Cdof,1:Ncoils,0:Cdof), lx(0:Cdof), ly(0:Cdof), lz(0:Cdof), fdvalue
  REAL, parameter :: fd = 1.0E-6

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!!$  call cpu_time(start)
!!$  call coilsep(2)
!!$  call cpu_time(finish)
!!$  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : Second order energy function 1 takes " ES23.15 " seconds.")') finish - start
!!$  
!!$  lbnorm = ccsep; l1B = t1C; l2B = t2C;
!!$
!!$  call cpu_time(start)
!!$  call coilsep(1)
!!$  call cpu_time(finish)
!!$  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : Second order energy function 2 takes " ES23.15 " seconds.")') finish - start
!!$
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "A5" : "4X,4(" : "ES23.15))') 'bnorm', ccsep, lbnorm, ccsep - lbnorm, (ccsep - lbnorm)/ccsep
!!$
!!$  do icoil = 1, Ncoils
!!$   do idof = 0, Cdof
!!$     if( abs(t1C(icoil,idof)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = (t1C(icoil,idof)-l1B(icoil,idof)) / t1C(icoil,idof)
!!$     endif
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "i4,4(" : "ES23.15))') icoil, idof, t1C(icoil,idof), l1B(icoil,idof), t1C(icoil,idof)-l1B(icoil,idof), rdiff
!!$    enddo
!!$   enddo
!!$   
!!$  jcoil = 1; jdof = 1;  
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : Checking the second derivatives using finite-difference method on t1E(",i4,", "i4")")') jcoil, jdof
!!$
!!$  do icoil = 1, Ncoils
!!$   do idof = 0, Cdof
!!$     if( abs(t2B(jcoil,jdof,icoil,idof)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = (t2B(jcoil,jdof,icoil,idof)-l2B(jcoil,jdof,icoil,idof)) / t2B(jcoil,jdof,icoil,idof)
!!$     endif
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "i4,4(" : "ES23.15))') icoil, idof, t2B(jcoil,jdof,icoil,idof), l2B(jcoil,jdof,icoil,idof), &
!!$                                                                                                   t2B(jcoil,jdof,icoil,idof)-l2B(jcoil,jdof,icoil,idof), rdiff
!!$    enddo
!!$   enddo
!!$  return

!-----------------------------------------

!!$  call cpu_time(start)
!!$  call bnormal(2)
!!$  call cpu_time(finish)
!!$  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : Second order energy function 1 takes " ES23.15 " seconds.")') finish - start
!!$  
!!$  lbx = coil(2)%Bx; lby = coil(2)%By; lbz = coil(2)%Bz
!!$
!!$  call cpu_time(start)
!!$  call bnormal2(2)
!!$  call cpu_time(finish)
!!$  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : Second order energy function 2 takes " ES23.15 " seconds.")') finish - start
!!$
!!$  do jdof = 0, Cdof
!!$   do idof = 0, Cdof
!!$     if( abs(coil(2)%Bx(idof,jdof)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = (coil(2)%Bx(idof,jdof)-lbx(idof,jdof)) / coil(2)%Bx(idof,jdof)
!!$     endif
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "i4,4(" : "ES23.15))') idof, jdof, coil(2)%Bx(idof,jdof), lbx(idof,jdof), &
!!$                                                                                                   coil(2)%Bx(idof,jdof)-lbx(idof,jdof), rdiff
!!$     if( abs(coil(2)%By(idof,jdof)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = (coil(2)%By(idof,jdof)-lby(idof,jdof)) / coil(2)%By(idof,jdof)
!!$     endif
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "i4,4(" : "ES23.15))') idof, jdof, coil(2)%By(idof,jdof), lby(idof,jdof), &
!!$                                                                                                   coil(2)%By(idof,jdof)-lby(idof,jdof), rdiff
!!$
!!$     if( abs(coil(2)%Bz(idof,jdof)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = (coil(2)%Bz(idof,jdof)-lbz(idof,jdof)) / coil(2)%Bz(idof,jdof)
!!$     endif
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "i4,4(" : "ES23.15))') idof, jdof, coil(2)%Bz(idof,jdof), lbz(idof,jdof), &
!!$                                                                                                   coil(2)%Bz(idof,jdof)-lbz(idof,jdof), rdiff
!!$    enddo
!!$   enddo

!!$ALLOCATE(mincc(1:Ncoils, 1:Ncoils))
!!$!----------------------------------------------
!!$  jcoil = 2; jdof = 2;   !chang here for testing different terms
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : Checking the second derivatives using finite-difference method on t1E(",i4,", "i4")")') jcoil, jdof
!!$
!!$  icoil = 1; iteta = 1; jzeta = 1
!!$  call pack( xdof )  !back up the current coils data
!!$  call epotent2(icoil, jcoil, icoil, lbz(0:Cdof,0:Cdof))
!!$  
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : coil# : dof : dof# :   numerical value"6X" :   fd-method value"6X" :   difference"11X" :   relative diff")')
!!$
!!$        
!!$!-------X_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%xc(idof) = coil(icoil)%xc(idof) + half * fd * ifd; call discretecoil
!!$           call epotent1(icoil,jcoil,lz(0:Cdof)); tmpE(ifd)            = lz(3)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = lbz(3,idof+1); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'XC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$     enddo
!!$
!!$!-------X_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%xs(idof) = coil(icoil)%xs(idof) + half * fd * ifd; call discretecoil
!!$           call epotent1(icoil,jcoil,lz(0:Cdof)); tmpE(ifd)            = lz(3)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = lbz(3,idof+NFcoil+2); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'XS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Y_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%yc(idof) = coil(icoil)%yc(idof) + half * fd * ifd; call discretecoil
!!$           call epotent1(icoil,jcoil,lz(0:Cdof)); tmpE(ifd)            = lz(3)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = lbz(3,idof+2*NFcoil+3); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'YC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Y_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%ys(idof) = coil(icoil)%ys(idof) + half * fd * ifd; call discretecoil
!!$           call epotent1(icoil,jcoil,lz(0:Cdof)); tmpE(ifd)            = lz(3)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = lbz(3,idof+3*NFcoil+4); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'YS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Z_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%zc(idof) = coil(icoil)%zc(idof) + half * fd * ifd; call discretecoil
!!$           call epotent1(icoil,jcoil,lz(0:Cdof)); tmpE(ifd)            = lz(3)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = lbz(3,idof+4*NFcoil+5); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'ZC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$!-------Z_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%zs(idof) = coil(icoil)%zs(idof) + half * fd * ifd; call discretecoil
!!$           call epotent1(icoil,jcoil,lz(0:Cdof)); tmpE(ifd)            = lz(3)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = lbz(3,idof+5*NFcoil+6); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'ZS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$ 
!!$  call unpack( xdof )  !recover data
!!$ 
!!$  return

!!$  icoil = 1; jcoil = 2; 
!!$  call pack( xdof )  !back up the current coils data
!!$  call epotent2(icoil, jcoil, icoil, lbz(0:Cdof,0:Cdof))
!!$  
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : coil# : dof : dof# :   numerical value"6X" :   fd-method value"6X" :   difference"11X" :   relative diff")')
!!$
!!$!-------X_COS-----------------
!!$  do idof = 0, NFcoil
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%xc(idof) = coil(icoil)%xc(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(1) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%xc(idof) = coil(icoil)%xc(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(2) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%xc(idof) = coil(icoil)%xc(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(3) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%xc(idof) = coil(icoil)%xc(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(4) = lbnorm
!!$
!!$     fdvalue = ( tmpE(1) + tmpE(4) - tmpE(2) -tmpE(3) ) / fd**2
!!$
!!$     tmpE(0) = lbz(3,idof+1); diff = tmpE(0) - fdvalue; rdiff = diff / tmpE(0)
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'XC', idof, tmpE(0), fdvalue, diff, rdiff
!!$  enddo
!!$!-------X_sin-----------------
!!$  do idof = 0, NFcoil
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%xs(idof) = coil(icoil)%xs(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(1) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%xs(idof) = coil(icoil)%xs(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(2) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%xs(idof) = coil(icoil)%xs(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(3) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%xs(idof) = coil(icoil)%xs(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(4) = lbnorm
!!$
!!$     fdvalue = ( tmpE(1) + tmpE(4) - tmpE(2) -tmpE(3) ) / fd**2
!!$
!!$     tmpE(0) = lbz(3,idof+NFcoil+2); diff = tmpE(0) - fdvalue; rdiff = diff / tmpE(0)
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'XS', idof, tmpE(0), fdvalue, diff, rdiff
!!$  enddo
!!$!-------Y_COS-----------------
!!$  do idof = 0, NFcoil
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%yc(idof) = coil(icoil)%yc(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(1) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%yc(idof) = coil(icoil)%yc(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(2) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%yc(idof) = coil(icoil)%yc(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(3) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%yc(idof) = coil(icoil)%yc(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(4) = lbnorm
!!$
!!$     fdvalue = ( tmpE(1) + tmpE(4) - tmpE(2) -tmpE(3) ) / fd**2
!!$
!!$     tmpE(0) = lbz(3,idof+2*NFcoil+3); diff = tmpE(0) - fdvalue; rdiff = diff / tmpE(0)
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'YC', idof, tmpE(0), fdvalue, diff, rdiff
!!$  enddo
!!$!-------Y_SIN-----------------
!!$  do idof = 0, NFcoil
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%ys(idof) = coil(icoil)%ys(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(1) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%ys(idof) = coil(icoil)%ys(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(2) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%ys(idof) = coil(icoil)%ys(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(3) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%ys(idof) = coil(icoil)%ys(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(4) = lbnorm
!!$
!!$     fdvalue = ( tmpE(1) + tmpE(4) - tmpE(2) -tmpE(3) ) / fd**2
!!$
!!$     tmpE(0) = lbz(3,idof+3*NFcoil+4); diff = tmpE(0) - fdvalue; rdiff = diff / tmpE(0)
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'YS', idof, tmpE(0), fdvalue, diff, rdiff
!!$  enddo
!!$!-------Z_COS-----------------
!!$  do idof = 0, NFcoil
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%zc(idof) = coil(icoil)%zc(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(1) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%zc(idof) = coil(icoil)%zc(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(2) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%zc(idof) = coil(icoil)%zc(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(3) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%zc(idof) = coil(icoil)%zc(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(4) = lbnorm
!!$
!!$     fdvalue = ( tmpE(1) + tmpE(4) - tmpE(2) -tmpE(3) ) / fd**2
!!$
!!$     tmpE(0) = lbz(3,idof+4*NFcoil+5); diff = tmpE(0) - fdvalue; rdiff = diff / tmpE(0)
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'ZC', idof, tmpE(0), fdvalue, diff, rdiff
!!$  enddo
!!$!-------Z_SIN-----------------
!!$  do idof = 0, NFcoil
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%zs(idof) = coil(icoil)%zs(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(1) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) - half * fd
!!$     coil(icoil)%zs(idof) = coil(icoil)%zs(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(2) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%zs(idof) = coil(icoil)%zs(idof) - half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(3) = lbnorm
!!$
!!$     call unpack ( xdof ); 
!!$     coil(icoil)%xc(2   ) = coil(icoil)%xc(2   ) + half * fd
!!$     coil(icoil)%zs(idof) = coil(icoil)%zs(idof) + half * fd
!!$     call discretecoil; call epotent0(icoil, jcoil, lbnorm); tmpE(4) = lbnorm
!!$
!!$     fdvalue = ( tmpE(1) + tmpE(4) - tmpE(2) -tmpE(3) ) / fd**2
!!$
!!$     tmpE(0) = lbz(3,idof+5*NFcoil+6); diff = tmpE(0) - fdvalue; rdiff = diff / tmpE(0)
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'ZS', idof, tmpE(0), fdvalue, diff, rdiff
!!$  enddo
!!$ 
!!$  call unpack( xdof )  !recover data
!!$ 
!!$  return

  !check lambda terms in eqarc.

!!$  call langrange(2)
!!$  icoil = 1
!!$
!!$  do idof = 0, NFcoil
!!$     do ifd = -1, 1, 2
!!$        coil(icoil)%lmdc(idof) = weight_eqarc; coil(icoil)%lmdc(idof) = coil(icoil)%lmdc(idof) + half * fd * ifd
!!$        call costfun(1)     ; tmpE(ifd)            = t1E(1,2)
!!$     enddo ! end do ifd
!!$
!!$     coil(icoil)%lmdc(idof) = weight_eqarc; call langrange(2)
!!$     tmpE(0) = dlc(icoil,idof,2); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$        rdiff = 0
!!$     else
!!$        rdiff = diff / tmpE(0)
!!$     endif
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'LC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$  enddo
!!$
!!$  do idof = 0, NFcoil
!!$     do ifd = -1, 1, 2
!!$        coil(icoil)%lmds(idof) = weight_eqarc; coil(icoil)%lmds(idof) = coil(icoil)%lmds(idof) + half * fd * ifd
!!$        call costfun(1)     ; tmpE(ifd)            = t1E(1,2)
!!$     enddo ! end do ifd
!!$
!!$     coil(icoil)%lmdc(idof) = weight_eqarc; call langrange(2)
!!$     tmpE(0) = dls(icoil,idof,2); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$        rdiff = 0
!!$     else
!!$        rdiff = diff / tmpE(0)
!!$     endif
!!$
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'LS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$  enddo


!----------------------------------------------
  jcoil = 2; jdof = 4;   !chang here for testing different terms
  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : Checking the second derivatives using finite-difference method on t1E(",i4,", "i4")")') jcoil, jdof

!!$  call bnormal(0)
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : bnormal value using old method with bscontant = ",ES23.15," : "ES23.15)') 1.0, bnorm
!!$  call bnormal2(0)
!!$  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : bnormal value using new method without NDcoil" 26X        " : "ES23.15)')       bnorm


  call cpu_time(start)
  call costfun(2)
  call cpu_time(finish)
  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : Second order derivatives of energy function takes " ES23.15 " seconds.")') finish - start

  call pack( bdof )  !back up the current coils data
  if(myid .eq. 0) write(ounit,'("fdcheck : " 10X " : The 2nd derivatives is over t1E("I3" , "I3 "). ")') jcoil, jdof
  if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : coil# : Fourier# :      numerical value"3X" :   fd-method value"6X" :   difference"11X" :   relative diff")')

  do idof = 1, Ndof

     do ifd = -1, 1, 2
        xdof = bdof
        xdof(idof) = xdof(idof) + half * fd * ifd
        call unpack ( xdof )
        call costfun(1)
        tmpE(ifd) = t1E(jcoil,jdof)
     enddo ! end do ifd

     call DoFconvert(idof, icoil, infou)
     tmpE(0) = t2E(jcoil, jdof, icoil,infou) ; diff = tmpE(0)- (tmpE(1) - tmpE(-1)) / fd

     if( abs(tmpE(0)) .lt. machprec ) then
         rdiff = 0
     else
         rdiff = diff / tmpE(0)
     endif
     
     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "I8, 4(" : "ES23.15))') &
                           icoil, infou, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
  enddo
  
!!$  do icoil = 1, Ncoils
!!$
!!$!-------current-----------------
!!$     do ifd = -1, 1, 2
!!$        call unpack ( xdof ); coil(icoil)%I = coil(icoil)%I + half * fd * ifd
!!$        call costfun(1)     ; tmpE(ifd)     = t1E(jcoil,jdof)
!!$     enddo ! end do ifd
!!$      
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,0) ; diff = tmpE(0)- (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "4X,4(" : "ES23.15))') &
!!$                           icoil, 'I', tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$        
!!$!-------X_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%xc(idof) = coil(icoil)%xc(idof) + half * fd * ifd
!!$           call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,idof+1); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$          icoil, 'XC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$     enddo
!!$
!!$!-------X_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%xs(idof) = coil(icoil)%xs(idof) + half * fd * ifd
!!$           call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,idof+NFcoil+2); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'XS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Y_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%yc(idof) = coil(icoil)%yc(idof) + half * fd * ifd
!!$           call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,idof+2*NFcoil+3); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'YC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Y_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%ys(idof) = coil(icoil)%ys(idof) + half * fd * ifd
!!$           call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,idof+3*NFcoil+4); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'YS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$
!!$!-------Z_COS-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%zc(idof) = coil(icoil)%zc(idof) + half * fd * ifd
!!$           call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,idof+4*NFcoil+5); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'ZC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$!-------Z_SIN-----------------
!!$     do idof = 0, NFcoil
!!$        do ifd = -1, 1, 2
!!$           call unpack ( xdof ); coil(icoil)%zs(idof) = coil(icoil)%zs(idof) + half * fd * ifd
!!$           call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$        enddo ! end do ifd
!!$
!!$     tmpE(0) = t2E(jcoil,jdof,icoil,idof+5*NFcoil+6); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$     if( abs(tmpE(0)) .lt. machprec ) then
!!$         rdiff = 0
!!$     else
!!$         rdiff = diff / tmpE(0)
!!$     endif
!!$     
!!$     if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                           icoil, 'ZS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$     enddo
!!$  
!!$  enddo ! end do icoil
!!$
!!$  if ( weight_eqarc .gt. 0.0 .and. Loptimize .eq. 3) then
!!$
!!$     call unpack ( xdof ); call costfun(2); call ntchdms(2)
!!$
!!$     do icoil = 1, Ncoils
!!$
!!$        do idof = 0, NFcoil
!!$           do ifd = -1, 1, 2
!!$              coil(icoil)%lmdc(idof) = 1.0; coil(icoil)%lmdc(idof) = coil(icoil)%lmdc(idof) + half * fd * ifd
!!$              call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$           enddo ! end do ifd
!!$
!!$           tmpE(0) = n2E(jcoil,jdof,icoil,idof+6*NFcoil+7); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$           if( abs(tmpE(0)) .lt. machprec ) then
!!$              rdiff = 0
!!$           else
!!$              rdiff = diff / tmpE(0)
!!$           endif
!!$
!!$           if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                icoil, 'LC', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$        enddo
!!$
!!$        do idof = 0, NFcoil
!!$           do ifd = -1, 1, 2
!!$              coil(icoil)%lmds(idof) = 1.0; coil(icoil)%lmds(idof) = coil(icoil)%lmds(idof) + half * fd * ifd
!!$              call costfun(1)     ; tmpE(ifd)            = t1E(jcoil,jdof)
!!$           enddo ! end do ifd
!!$
!!$           tmpE(0) = n2E(jcoil,jdof,icoil,idof+7*NFcoil+8); diff = tmpE(0) - (tmpE(1) - tmpE(-1)) / fd
!!$
!!$           if( abs(tmpE(0)) .lt. machprec ) then
!!$              rdiff = 0
!!$           else
!!$              rdiff = diff / tmpE(0)
!!$           endif
!!$
!!$           if( myid.eq.0 ) write(ounit,'("fdcheck : " 10X " : "I5" : "A3" : "i4,4(" : "ES23.15))') &
!!$                icoil, 'LS', idof, tmpE(0), (tmpE(1) - tmpE(-1)) / fd, diff, rdiff
!!$
!!$        enddo
!!$
!!$     enddo ! end do icoil
!!$
!!$  endif

  call unpack( xdof )  !recover data
 
  return
  
end subroutine testderivs2

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine NAGderiv1

  use kmodule, only: Ndof, ounit, myid, Idisplay
  implicit none
  include "mpif.h"

  Real    :: f
  Integer :: ifail, liw, lw, n

  Real    :: g(Ndof), w(3*Ndof), x(Ndof)
  Integer :: iw(1)

  external :: funct

  if (myid .eq. 0) Write (ounit,'("E04HCF  : " 10x " : " A )') 'Using NAG_E04HCF to test the 1sd derivatives'

  n = Ndof; lw = 3*n; liw = 1

  call pack(x(1:n))
  ifail = Idisplay
  Call e04hcf(n,funct,x,f,g,iw,liw,w,lw,ifail)

  if (ifail==0) Then
     if (myid .eq. 0) Write (ounit,'("E04HCF  : " 10x " : " A )') "1st derivatives are consistent with function values"
  else
     if (myid .eq. 0) Write (ounit,'("E04HCF  : " 10x " : " A )') 'Probable error in calculation of 1st derivatives'
  endif

 return

end subroutine NAGderiv1
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
subroutine NAGderiv2

  use kmodule, only: Ndof, ounit, myid, Idisplay
  implicit none
  include "mpif.h"

  Real    :: f
  Integer :: ifail, liw, lw, n, lh

  Real    :: g(Ndof), hesd(Ndof), hesl(Ndof*(Ndof-1)/2), w(5*Ndof), x(Ndof)
  Integer :: iw(1)

  external :: funct, h

  if (myid .eq. 0) Write (ounit,'("E04HDF  : " 10x " : " A )') 'Using NAG_E04HDF to test the 2nd derivatives'

  n = Ndof; lw = 5*n; liw = 1; lh = n*(n-1)/2

  if(lh .eq. 0) lh = 1

  call pack(x(1:n))

  ifail = Idisplay
  Call e04hcf(n,funct,x,f,g,iw,liw,w,lw,ifail)

  ifail = Idisplay
  Call e04hdf(n,funct,h,x,g,hesl,lh,hesd,iw,liw,w,lw,ifail)
 
  if (ifail==0) Then
     if (myid .eq. 0) Write (ounit,'("E04HDF  : " 10x " : " A )') "2nd derivatives are consistent with function values"
  elseif ( ifail == 2) then
     if (myid .eq. 0) Write (ounit,'("E04HDF  : " 10x " : " A )') 'Probable error in calculation of 1st derivatives'
  endif

 return

end subroutine NAGderiv2
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

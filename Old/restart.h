
!title (restart) ! Writes output and restart files.

!latex \briefly{Writes output and restart files.}

!latex \calledby{\link{knotopt} and \link{descent} }
!l tex \calls{\link{iccoil},\link{knotxx}}

!latex \tableofcontents

!latex \subsection{h5 format}
!latex \bi
!latex \item[1.] All the restart information etc. is written to file.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine restart( irestart )
  
  use kmodule

  use hdf5
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  INTEGER, intent(in):: irestart

  INTEGER            :: astat, ierr, ii, jj, ifail, llmodnp, idof, mm

  REAL               :: tt, aa, teta, zeta, ax(1:3), at(1:3), az(1:3), xx(1:3), xt(1:3), xz(1:3)
  REAL               :: x(0:1), y(0:1), z(0:1)

! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume
  
  CHARACTER          :: suffix*3, srestart*6

  REAL, allocatable  :: perA(:)

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
!  do icoil = 1, Ncoils
!   
!   if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop; 14 Apr 16;
!   
!   do ii = 0, NDcoil 
!    ifail = 0 ; tt = ii * pi2 / NDcoil  ; call iccoil( tt, x(0:1), y(0:1), z(0:1), ifail )
!    coil(icoil)%xx(ii) = x(0)
!    coil(icoil)%yy(ii) = y(0)
!    coil(icoil)%zz(ii) = z(0)
!   enddo ! end of do ii; 14 Apr 16;
!   
!   coilspace(icoil,1) = coil(icoil)%I
!   coilspace(icoil,2) = coil(icoil)%L
!
!  enddo ! end of do icoil; 14 Apr 16;
!  
!  do icoil = 1, Ncoils ; llmodnp = modulo(icoil-1,ncpu)
!   RlBCAST( coil(icoil)%xx(0:NDcoil), NDcoil+1, llmodnp )
!   RlBCAST( coil(icoil)%yy(0:NDcoil), NDcoil+1, llmodnp )
!   RlBCAST( coil(icoil)%zz(0:NDcoil), NDcoil+1, llmodnp )
!  enddo
    
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if(itau .gt. Ntauout) itau = Ntauout

  if( myid.ne.0 ) return

  SALLOCATE(perA, (1:Tdof), zero)
  if(allocated(coilspace)) then                                           !for saving coil current and coil length
  idof = 0
  do icoil = 1, Ncoils
      
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%I     ; perA(idof) = zero
                                                                    
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%xc( 0); perA(idof) = zero
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%yc( 0); perA(idof) = zero
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%zc( 0); perA(idof) = zero
                                                                    
      do mm = 1, NFcoil                                             
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%xc(mm); perA(idof) =   mm*coil(icoil)%xs(mm)
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%yc(mm); perA(idof) =   mm*coil(icoil)%ys(mm)
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%zc(mm); perA(idof) =   mm*coil(icoil)%zs(mm)
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%xs(mm); perA(idof) = - mm*coil(icoil)%xc(mm)
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%ys(mm); perA(idof) = - mm*coil(icoil)%yc(mm) 
        idof = idof + 1 ; coilspace(itau, idof) = coil(icoil)%zs(mm); perA(idof) = - mm*coil(icoil)%zc(mm)
      enddo

   enddo
  
   FATAL( restart , idof .ne. Tdof, counting error in packing to coilspace )
      
  endif
   
  do icoil = 1, Ncoils
   if (.not. allocated(coilsI) ) allocate( coilsI(1:Ncoils) )             !scale back the current
   coilsI(icoil) = coil(icoil)%I / antibscont
  enddo
 
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  nrestart = nrestart + 1 ; write(srestart,'(i6.6)') nrestart
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5open_f )
  
  call h5fcreate_f( trim(ext)//".fo.h5", H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
  FATAL( restart, hdfier.ne.0, error calling h5fcreate_f )

  HWRITEIV( 1                           ,   Idisplay                                ,   Idisplay                                                  )
  HWRITEIV( 1                           ,   Isymmetric                              ,   Isymmetric                                                )
  HWRITEIV( 1                           ,   Itopology                               ,   Itopology                                                 )
  HWRITERV( 1                           ,   knotsurf                                ,   knotsurf                                                  )
  HWRITERV( 1                           ,   ellipticity                             ,   ellipticity                                               )
  HWRITEIV( 1                           ,   Linitialize                             ,   Linitialize                                               )
  HWRITERV( 1                           ,   Rmaj                                    ,   Rmaj                                                      )
  HWRITERV( 1                           ,   rmin                                    ,   rmin                                                      )
  HWRITEIV( 1                           ,   Ic                                      ,   Ic                                                        )
  HWRITERV( 1                           ,   Io                                      ,   Io                                                        )
  HWRITERV( 1                           ,   Iw                                      ,   Iw                                                        )
  HWRITEIV( 1                           ,   Lc                                      ,   Lc                                                        )
  HWRITERV( 1                           ,   Lo                                      ,   Lo                                                        )
  HWRITERV( 1                           ,   Lw                                      ,   Lw                                                        )
  HWRITEIV( 1                           ,   Loptimize                               ,   Loptimize                                                 )
  HWRITEIV( 1                           ,   Lnormalize                              ,   Lnormalize                                                )
  HWRITEIV( 1                           ,   NFcoil                                  ,   NFcoil                                                    )
  HWRITEIV( 1                           ,   NDcoil                                  ,   NDcoil                                                    )
  HWRITERV( 1                           ,   weight_bnorm                            ,   tmpw_bnorm                                                )
  HWRITERV( 1                           ,   weight_tflux                            ,   tmpw_tflux                                                )
  HWRITERV( 1                           ,   target_tflux                            ,   tmpt_tflux                                                )
  HWRITERV( 1                           ,   weight_ttlen                            ,   tmpw_ttlen                                                )
  HWRITERV( 1                           ,   weight_eqarc                            ,   tmpw_eqarc                                                )
  HWRITERV( 1                           ,   weight_ccsep                            ,   tmpw_ccsep                                                )
  HWRITERV( 1                           ,   tauend                                  ,   tauend                                                    )
  HWRITERV( 1                           ,   tautol                                  ,   tautol                                                    )
  HWRITEIV( 1                           ,   Ntauout                                 ,   Ntauout                                                   )
  HWRITEIV( 1                           ,   Savfreq                                 ,   Savfreq                                                   )
  HWRITEIV( 1                           ,   Nteta                                   ,   Nteta                                                     )
  HWRITEIV( 1                           ,   Nzeta                                   ,   Nzeta                                                     )
  HWRITERV( 1                           ,   absacc                                  ,   absacc                                                    )
  HWRITERV( 1                           ,   absreq                                  ,   absreq                                                    )
  HWRITERV( 1                           ,   relreq                                  ,   relreq                                                    )
  HWRITERV( 1                           ,   xtol                                    ,   xtol                                                      )
  HWRITERV( 1                           ,   eta                                     ,   eta                                                       )
  HWRITERV( 1                           ,   stepmx                                  ,   stepmx                                                    )
  HWRITEIV( 1                           ,   Mpol                                    ,   Mpol                                                      )
  HWRITEIV( 1                           ,   Ntor                                    ,   Ntor                                                      )
  HWRITEIV( 1                           ,   Lpoincare                               ,   Lpoincare                                                 )
  HWRITERV( 1                           ,   odetol                                  ,   odetol                                                    )
  HWRITEIV( 1                           ,   Ppts                                    ,   Ppts                                                      )
  HWRITEIV( 1                           ,   Ptrj                                    ,   Ptrj                                                      )
  HWRITEIV( 1                           ,   iphi                                    ,   iphi                                                      )
  HWRITERV( 1                           ,   bstol                                   ,   bstol                                                     )
  HWRITEIV( 1                           ,   bsnlimit                                ,   bsnlimit                                                  )

  HWRITEIV( 1                           ,   bmn                                     ,   bmn                                                       )
  if( bmn.gt.0 ) then 
  HWRITEIV( bmn                         ,   bim                                     ,   bim(1:bmn)                                                )
  HWRITEIV( bmn                         ,   bin                                     ,   bin(1:bmn)                                                )
  HWRITERV( bmn                         ,   Rbc                                     ,   Rbc(1:bmn)                                                )
  HWRITERV( bmn                         ,   Rbs                                     ,   Rbs(1:bmn)                                                )
  HWRITERV( bmn                         ,   Zbc                                     ,   Zbc(1:bmn)                                                )
  HWRITERV( bmn                         ,   Zbs                                     ,   Zbs(1:bmn)                                                )
  endif

  HWRITEIV( 1                           ,   Ncoils                                  ,   Ncoils                                                    )

  HWRITEIV( 1                           ,   itau                                    ,   itau                                                      )

  HWRITERA( itau+1, 10                  ,   evolution                               ,   evolution(0:itau,0:9)                                     ) ! itau-1

  HWRITERA( Ndof, 6                     ,   deriv                                   ,   deriv(1:Ndof, 0:5)                                        )

  HWRITERA( itau+1, Tdof                ,   coilspace                               ,   coilspace(0:itau, 1:Tdof)                                 )
  HWRITERV( Tdof                        ,   perA                                    ,   perA(1:Tdof)                                              )
  
  HWRITERA( 1+Nteta,1+Nzeta             ,   xsurf                                   ,   surf(1)%xx(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   ysurf                                   ,   surf(1)%yy(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   zsurf                                   ,   surf(1)%zz(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   nx                                      ,   surf(1)%nx(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   ny                                      ,   surf(1)%ny(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   nz                                      ,   surf(1)%nz(0:Nteta,0:Nzeta)                               )

  if (allocated(tbn)) then
  HWRITERA( 1+Nteta,1+Nzeta             ,   tgtBn                                   ,   surf(1)%bnt(0:Nteta,0:Nzeta)                              )
  HWRITERA( 1+Nteta,1+Nzeta             ,   curBn                                   ,           tbn(0:Nteta,0:Nzeta)                              )

 !HWRITERA( 1+Cdof ,1+Cdof              ,   Bdx                                     ,   coil(icoil)%Bx(0:Cdof,0:Cdof)                             )
  HWRITERA( 1+Nteta,1+Nzeta             ,   Bx                                      ,       SaveBx(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   By                                      ,       SaveBy(0:Nteta,0:Nzeta)                               )
  HWRITERA( 1+Nteta,1+Nzeta             ,   Bz                                      ,       SaveBz(0:Nteta,0:Nzeta)                               )
  endif

  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )
  
  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  do icoil = 1, Ncoils
   write(suffix,'(i3.3)') icoil
   open( lunit, file=".fo.coil."//suffix, status="unknown" )
   write(lunit,1000) coil(icoil)%N, coil(icoil)%D
   write(lunit,1001) coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%Io, coil(icoil)%Iw
   write(lunit,1001) coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo, coil(icoil)%Lw
   write(lunit,1002) coil(icoil)%xc(0:NFcoil)
   write(lunit,1002) coil(icoil)%xs(0:NFcoil)
   write(lunit,1002) coil(icoil)%yc(0:NFcoil)
   write(lunit,1002) coil(icoil)%ys(0:NFcoil)
   write(lunit,1002) coil(icoil)%zc(0:NFcoil)
   write(lunit,1002) coil(icoil)%zs(0:NFcoil)
   close(lunit)
  enddo

1000 format(            2i9          )
1001 format(    es23.15, i2, 2es23.15)
1002 format(9999es23.15              )

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
 !select case( Itopology )
 !case( 1 )
 ! open( lunit, file=trim(ext)//".fo.knot", status="unknown" )
 ! write(lunit,'(i9,es23.15)')       NFcoil, knotphase
 ! write(lunit,'(999es23.15)') xkc(0:NFcoil)
 ! write(lunit,'(999es23.15)') xks(0:NFcoil)
 ! write(lunit,'(999es23.15)') ykc(0:NFcoil)
 ! write(lunit,'(999es23.15)') yks(0:NFcoil)
 ! write(lunit,'(999es23.15)') zkc(0:NFcoil)
 ! write(lunit,'(999es23.15)') zks(0:NFcoil)
 ! close(lunit)
 !end select

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  if( irestart.eq. 1 ) then

  open( lunit, file="."//trim(ext)//".fo.filaments."//srestart, status="unknown", form="unformatted" )  
  write(lunit) Ncoils, NDcoil
  do icoil = 1, Ncoils
   write(lunit) coil(icoil)%xx(0:NDcoil)
   write(lunit) coil(icoil)%yy(0:NDcoil)
   write(lunit) coil(icoil)%zz(0:NDcoil)
  enddo ! end of do icoil; 14 Apr 16;
  close( lunit )

  endif

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  
  if( irestart.eq.1 ) then

  open(lunit,file="coils."//trim(ext), status="unknown", form="formatted" )
  write(lunit,'("periods "I3)') bNfp
  write(lunit,'("begin filament")')
  write(lunit,'("mirror NIL")')
  do icoil = 1, Ncoils
   do ii = 0, NDcoil-1 ; write(lunit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), coilsI(icoil)
   enddo
   ;  ii =    NDcoil   ; write(lunit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), zero         , icoil, coil(icoil)%name
  enddo ! end of do icoil; 14 Apr 16;
  write(lunit,'("end")')
  close(lunit)

  endif
  
1010 format(4es23.15,:,i9,"  ",a10)


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  return
  
end subroutine restart

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_plasma
!-------------------------------------------------------------------------------!
! write down the unpdated plasma boundary information;                          !
! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
!-------------------------------------------------------------------------------!
  use kmodule, only : zero, half, pi2, myid, ncpu, ounit, lunit, ext, &
                      bmn, bNfp, nbf, bim, bin, bnim, bnin, Rbc, Rbs, Zbc, Zbs, bnc, bns, &
                      Nteta, Nzeta, surf
  
  implicit none  
  include "mpif.h"
!-------------------------------------------------------------------------------
  INTEGER :: imn
!-------------------------------------------------------------------------------
  if(myid .ne. 0) return

  open(lunit, file=trim(ext)//".plasma", status='unknown', action='write')

  write(lunit,*      ) "#bmn bNfp nbf"
  write(lunit,'(3I)' ) bmn, bNfp, nbf

  write(lunit,*      ) "#------- plasma boundary------"
  write(lunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
  do imn = 1, bmn
     write(lunit,'(2I, 4ES15.6)') bin(imn), bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
  enddo

  write(lunit,*      ) "#-------Bn harmonics----------"
  write(lunit,*      ) "#  n  m  bnc   bns"
  if (nbf .gt. 0) then
  do imn = 1, nbf
     write(lunit,'(2I, 2ES15.6)') bnin(imn), bnim(imn), bnc(imn), bns(imn)
  enddo
  else
     write(lunit,'(2I, 2ES15.6)') 0, 0, 0.0, 0.0
  endif

  close(lunit)
END SUBROUTINE write_plasma

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

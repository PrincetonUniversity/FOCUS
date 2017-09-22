!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

!title (restart) ! Writes output and restart files.

!latex \briefly{Writes output and restart files.}

!latex \calledby{\link{solvers}}
!l tex \calls{}

!latex \tableofcontents

!latex \section{h5 format}
!latex \bi
!latex \item[1.] All the restart information etc. is written to file.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

subroutine archive
  
  use globals
  
  use hdf5
  
  implicit none
  
  include "mpif.h"
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  INTEGER            :: ii, icoil, ierr
  REAL               :: tnow
  CHARACTER(LEN=10)  :: version='h1.0.0'

  CHARACTER          :: srestart*6
    
! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank, imn
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  REAL               :: tvolume
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  if( myid.ne.0 ) return

  write(srestart,'(i6.6)') iarchive

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5open_f )
  
  call h5fcreate_f( trim(hdf5file), H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
  FATAL( restart, hdfier.ne.0, error calling h5fcreate_f )
  
  HWRITECH( 10               ,   version       ,   version                    )
  
  HWRITEIV( 1                ,   Icheck        ,   Icheck                     )

  HWRITEIV( 1                ,   Isurface      ,   Isurface                   )
  HWRITERV( 1                ,   minorrad      ,   minorrad                   )
  HWRITERV( 1                ,   ellipticity   ,   ellipticity                )

  HWRITEIV( 1                ,   Nteta         ,   Nteta                      )
  HWRITEIV( 1                ,   Nzeta         ,   Nzeta                      )
  HWRITEIV( 1                ,   Initialize    ,   Initialize                 )
  HWRITEIV( 1                ,   Ncoils        ,   Ncoils                     )
  HWRITERV( 1                ,   init_current  ,   init_current               )
  HWRITERV( 1                ,   init_radius   ,   init_radius                )
  HWRITEIV( 1                ,   IsVaryCurrent ,   IsVaryCurrent              )
  HWRITEIV( 1                ,   IsVaryGeometry,   IsvaryGeometry             )
  HWRITEIV( 1                ,   NFcoil        ,   NFcoil                     )
  HWRITEIV( 1                ,   Nsegments     ,   Nsegments                  )
  HWRITEIV( 1                ,   case_length   ,   case_length                )
  HWRITERV( 1                ,   weight_bnorm  ,   weight_bnorm               )
  HWRITERV( 1                ,   weight_tflux  ,   weight_tflux               )
  HWRITERV( 1                ,   target_tflux  ,   target_tflux               )
  HWRITERV( 1                ,   weight_ttlen  ,   weight_ttlen               )
  HWRITERV( 1                ,   target_length ,   target_length              )

  HWRITERV( 1                ,   converged     ,   converged                  )

  HWRITERV( 1                ,    tauend       ,    tauend                    )
  HWRITEIV( 1                ,   Ntauout       ,   Ntauout                    )
  HWRITERV( 1                ,    tautol       ,    tautol                    )

  HWRITERV( 1                ,   friction      ,   friction                   )
  HWRITEIV( 1                ,   Iminimize     ,   Iminimize                  )

! HWRITEIV( 1                ,   Nfp           ,   Nfp                          )
  HWRITERV( 1                ,   area          ,   surf%area                    )
  HWRITERA( Nteta,Nzeta      ,   xsurf         ,   surf%xx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ysurf         ,   surf%yy(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   zsurf         ,   surf%zz(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nx            ,   surf%nx(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   ny            ,   surf%ny(0:Nteta-1,0:Nzeta-1) )
  HWRITERA( Nteta,Nzeta      ,   nz            ,   surf%nz(0:Nteta-1,0:Nzeta-1) )

  HWRITERV( 1                ,   totlengt      ,   totlengt(0)                  )
  
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )
  
  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  open( wunit, file=trim(coilfile), status="unknown", form="formatted") ! FOCUS-format; coilfile = ext.focus; defined in initial;
  
  write(wunit, *) "# Total number of coils"
  write(wunit, '(8X,I6)') Ncoils
  
  do icoil = 1, Ncoils
   
   write(wunit, *) "#-----------------", icoil, "---------------------------" 
   
   write(wunit, *) "#coil_type     coil_name"
   write(wunit,'(3X,I3,4X, A10)') 1, coil(icoil)%name
   
   write(wunit, '(3(A6, A15, 8X))') " #Nseg", "current",  "Ifree", "Length", "Lfree", "target_length"
   write(wunit,'(2X,I4,ES23.15,3X,I3,ES23.15,3X,I3,ES23.15)') Ns, coil(icoil)%I, coil(icoil)%Ifree, coil(icoil)%dL(0), coil(icoil)%Lfree, coil(icoil)%Lo
   
   write(wunit, *) "#NFcoil"
   write(wunit, '(I3)') coil(icoil)%NF
   write(wunit, *) "#Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs) "
   write(wunit, 1000) coil(icoil)%xc(0:coil(icoil)%NF)
   write(wunit, 1000) coil(icoil)%xs(0:coil(icoil)%NF)
   write(wunit, 1000) coil(icoil)%yc(0:coil(icoil)%NF)
   write(wunit, 1000) coil(icoil)%ys(0:coil(icoil)%NF)
   write(wunit, 1000) coil(icoil)%zc(0:coil(icoil)%NF)
   write(wunit, 1000) coil(icoil)%zs(0:coil(icoil)%NF)
   
  enddo ! end of do icoil = 1, Ncoils ;
  
  close(wunit)
  
1000 format(999es23.15)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  open(funit,file=trim(outcoils), status="unknown", form="formatted" ) ! xgrid-format; outcoils = ext.coils; defined in initial;
  
  write(funit,'("periods "I3)') Nfp
  write(funit,'("begin filament")')
  write(funit,'("mirror NIL")')
  
  do icoil = 1, Ncoils
   
   do ii = 0, Ns-1
    write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), coil(icoil)%I
   enddo
   
   ii =  Ns
   
   write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), zero, icoil, coil(icoil)%name
   
  enddo
  
  write(funit,'("end")')
  
  close(funit)
  
1010 format(4es23.15,:,i9,"  ",a10)
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
  open( funit, file="."//trim(ext)//".fo.filaments."//srestart, status="unknown", form="unformatted" )  
  
  write(funit) Ncoils, Ns
  
  do icoil = 1, Ncoils
   write(funit) coil(icoil)%xx(0:Ns)
   write(funit) coil(icoil)%yy(0:Ns)
   write(funit) coil(icoil)%zz(0:Ns)
  enddo
  
  close( funit )
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
    
  iarchive = iarchive + 1

  return
  
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-
  
end subroutine archive

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

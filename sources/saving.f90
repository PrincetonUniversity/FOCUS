
!title (restart) ! Writes output and restart files.

!latex \briefly{Writes output and restart files.}

!latex \calledby{\link{solvers}}
!l tex \calls{}

!latex \tableofcontents

!latex \section{h5 format}
!latex \bi
!latex \item[1.] All the restart information etc. is written to file.
!latex \ei

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

subroutine saving

  use globals
  use mpi
  use hdf5

  implicit none

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!


  INTEGER            :: ii, jj, icoil, NF, icpu

  ! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume

  CHARACTER          :: srestart*6




  !--------------------------write focus coil file-----------------------------------------
  if( save_coils == 1 ) then
     if (myid==0) then
        open( wunit, file=trim(out_focus), status="unknown", form="formatted")
        write(wunit, '("Total number of coils")') 
        write(wunit, '(2X,I8)') Ncoils_total ! note the fixed coils are not written
#ifdef TOPO
        write(wunit, '(2(A4,", "), A13, ", ", 3(A15,", "), 2(A2,", ",A15,", ",A15,", "))') &
             "type", "symm.", "coilname", "ox", "oy", "oz", "Ic", "M_0", "pho", "Lc", "mp", "mt"
#endif
        close(wunit)
     endif

     call MPI_barrier( MPI_COMM_WORLD, ierr )

     do icpu = 0, ncpu-1
        if (myid/=0 .and. myid == icpu) then                 ! each cpu write the data independently
           open( wunit, file=trim(out_focus), status="old", position="append", action="write")

           do icoil = 1, Ncoils
#ifdef TOPO
              write(wunit, '(2(I4, ", ")A13,", ",3(ES15.8,", "),2(I2,", ",ES15.8,", ",ES15.8,", "))') &
                   coil(icoil)%itype, coil(icoil)%symmetry, coil(icoil)%name, &
                   coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                   coil(icoil)%Ic, coil(icoil)%moment, coil(icoil)%pho, &
                   coil(icoil)%Lc, coil(icoil)%mp, coil(icoil)%mt
#else
              write(wunit, *) "#-----------------", icoil, "---------------------------" 
              write(wunit, *) "#coil_type     coil_name"
              write(wunit,'(3X,I3,4X,I3,4X,A15)') coil(icoil)%itype, coil(icoil)%symmetry, coil(icoil)%name

              select case (coil(icoil)%itype)
              case (1)
                 write(wunit, '(3(A6, A15, 8X))') " #Nseg", "current",  "Ifree", "Length", "Lfree", "target_length"
                 write(wunit,'(2X, I4, ES23.15, 3X, I3, ES23.15, 3X, I3, ES23.15)') &
                      coil(icoil)%NS, coil(icoil)%I, coil(icoil)%Ic, coil(icoil)%L, coil(icoil)%Lc, coil(icoil)%Lo
                 NF = FouCoil(icoil)%NF ! shorthand;
                 write(wunit, *) "#NFcoil"
                 write(wunit, '(I3)') NF
                 write(wunit, *) "#Fourier harmonics for coils ( xc; xs; yc; ys; zc; zs) "
                 write(wunit, 1000) FouCoil(icoil)%xc(0:NF)
                 write(wunit, 1000) FouCoil(icoil)%xs(0:NF)
                 write(wunit, 1000) FouCoil(icoil)%yc(0:NF)
                 write(wunit, 1000) FouCoil(icoil)%ys(0:NF)
                 write(wunit, 1000) FouCoil(icoil)%zc(0:NF)
                 write(wunit, 1000) FouCoil(icoil)%zs(0:NF)
              case (2) 
                 write(wunit, *) "#  Lc  ox   oy   oz  Ic  I  mt  mp"
                 write(wunit,'(2(I3, 3ES23.15))') coil(icoil)%Lc, coil(icoil)%ox, coil(icoil)%oy, coil(icoil)%oz, &
                      coil(icoil)%Ic, coil(icoil)%I , coil(icoil)%mt, coil(icoil)%mp  
              case (3)
                 write(wunit, *) "# Ic     I    Lc  Bz  (Ic control I; Lc control Bz)"
                 write(wunit,'(I3, ES23.15, I3, ES23.15)') coil(icoil)%Ic, coil(icoil)%I, &
                      coil(icoil)%Lc, coil(icoil)%Bz
              case default
                 FATAL(restart, .true., not supported coil types)
              end select
#endif
           enddo
           close(wunit)
1000       format(9999ES23.15)
        endif
        call MPI_barrier( MPI_COMM_WORLD, ierr )
     enddo
  endif

  if( myid.ne.0 ) return

  !SALLOCATE(perA, (1:Tdof), zero)

  !store the current derivatives;
  if (allocated(deriv)) then
     deriv = zero
     if(allocated(t1E)) deriv(1:Ndof,0) = t1E(1:Ndof)
     if(allocated(t1B)) deriv(1:Ndof,1) = t1B(1:Ndof)
     if(allocated(t1F)) deriv(1:Ndof,2) = t1F(1:Ndof)
     if(allocated(t1L)) deriv(1:Ndof,3) = t1L(1:Ndof)
     if(allocated(t1S)) deriv(1:Ndof,4) = t1S(1:Ndof)
     if(allocated(t1C)) deriv(1:Ndof,5) = t1C(1:Ndof)
     if(allocated(t1H)) deriv(1:Ndof,6) = t1H(1:Ndof)
  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  call h5open_f( hdfier ) ! initialize Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5open_f )

  call h5fcreate_f( trim(hdf5file), H5F_ACC_TRUNC_F, file_id, hdfier ) ! create new file;
  FATAL( restart, hdfier.ne.0, error calling h5fcreate_f )

  HWRITECH( 10               ,   version       ,   version                       )

  !INPUT namelist;
  HWRITEIV( 1                ,   IsQuiet       ,   IsQuiet                       )
  HWRITEIV( 1                ,   IsSymmetric   ,   IsSymmetric                   )
  HWRITECH( 100              ,   input_surf    ,   input_surf                    )
  HWRITECH( 100              ,   input_coils   ,   input_coils                   )
  HWRITECH( 100              ,   input_harm    ,   input_harm                    )
  HWRITEIV( 1                ,   case_surface  ,   case_surface                  )
  HWRITERV( 1                ,   knotsurf      ,   knotsurf                      )
  HWRITEIV( 1                ,   Nteta         ,   Nteta                         )
  HWRITEIV( 1                ,   Nzeta         ,   Nzeta                         )
  HWRITEIV( 1                ,   case_init     ,   case_init                     )
  HWRITEIV( 1                ,   case_coils    ,   case_coils                    )
  HWRITEIV( 1                ,   Ncoils        ,   Ncoils                        )
  HWRITERV( 1                ,   init_current  ,   init_current                  )
  HWRITERV( 1                ,   init_radius   ,   init_radius                   )
  HWRITEIV( 1                ,   IsVaryCurrent ,   IsVaryCurrent                 )
  HWRITEIV( 1                ,   IsVaryGeometry,   IsvaryGeometry                )
  HWRITEIV( 1                ,   NFcoil        ,   NFcoil                        )
  HWRITEIV( 1                ,   Nseg          ,   Nseg                          )
  HWRITEIV( 1                ,   case_optimize ,   case_optimize                 )
  HWRITERV( 1                ,   exit_xtol     ,   exit_tol                      )
  HWRITEIV( 1                ,   IsNormalize   ,   IsNormalize                   )
  HWRITEIV( 1                ,   ISNormWeight  ,   IsNormWeight                  )
  HWRITEIV( 1                ,   case_bnormal  ,   case_bnormal                  )
  HWRITEIV( 1                ,   case_length   ,   case_length                   )
  HWRITERV( 1                ,   weight_bnorm  ,   weight_bnorm                  )
  HWRITERV( 1                ,   weight_bharm  ,   weight_bharm                  )
  HWRITERV( 1                ,   weight_tflux  ,   weight_tflux                  )
  HWRITERV( 1                ,   target_tflux  ,   target_tflux                  )
  HWRITERV( 1                ,   weight_ttlen  ,   weight_ttlen                  )
  HWRITERV( 1                ,   target_length ,   target_length                 )
  HWRITERV( 1                ,   weight_specw  ,   weight_specw                  )
  HWRITERV( 1                ,   weight_ccsep  ,   weight_ccsep                  )
  HWRITERV( 1                ,   weight_gnorm  ,   weight_gnorm                  )
  HWRITERV( 1                ,   weight_inorm  ,   weight_inorm                  )
  HWRITERV( 1                ,   weight_mnorm  ,   weight_mnorm                  )
  HWRITERV( 1                ,   DF_tausta     ,   DF_tausta                     )
  HWRITERV( 1                ,   DF_tauend     ,   DF_tauend                     )
  HWRITERV( 1                ,   DF_xtol       ,   DF_xtol                       )
  HWRITEIV( 1                ,   DF_maxiter    ,   DF_maxiter                    )
  HWRITEIV( 1                ,   CG_maxiter    ,   CG_maxiter                    )
  HWRITERV( 1                ,   CG_xtol       ,   CG_xtol                       )
  HWRITERV( 1                ,   CG_wolfe_c1   ,   CG_wolfe_c1                   )
  HWRITERV( 1                ,   CG_wolfe_c2   ,   CG_wolfe_c2                   )
  HWRITEIV( 1                ,   HN_maxiter    ,   HN_maxiter                    )
  HWRITERV( 1                ,   HN_xtol       ,   HN_xtol                       )
  HWRITERV( 1                ,   HN_factor     ,   HN_factor                     )
  HWRITEIV( 1                ,   TN_maxiter    ,   TN_maxiter                    )
  HWRITEIV( 1                ,   TN_reorder    ,   TN_reorder                    )
  HWRITERV( 1                ,   TN_xtol       ,   TN_xtol                       )
  HWRITERV( 1                ,   TN_cr         ,   TN_cr                         )
  HWRITEIV( 1                ,   case_postproc ,   case_postproc                 )
  HWRITEIV( 1                ,   save_freq     ,   save_freq                     )
  HWRITEIV( 1                ,   save_coils    ,   save_coils                    )
  HWRITEIV( 1                ,   save_harmonics,   save_harmonics                )
  HWRITEIV( 1                ,   save_filaments,   save_filaments                )
  HWRITEIV( 1                ,   update_plasma ,   update_plasma                 )
  HWRITERV( 1                ,   pp_phi        ,   pp_phi                        )
  HWRITERV( 1                ,   pp_raxis      ,   pp_raxis                      )
  HWRITERV( 1                ,   pp_zaxis      ,   pp_zaxis                      )
  HWRITERV( 1                ,   pp_rmax       ,   pp_rmax                       )
  HWRITERV( 1                ,   pp_zmax       ,   pp_zmax                       )
  HWRITEIV( 1                ,   pp_ns         ,   pp_ns                         )
  HWRITEIV( 1                ,   pp_maxiter    ,   pp_maxiter                    )
  HWRITERV( 1                ,   pp_xtol       ,   pp_xtol                       )

  HWRITEIV( 1                ,   Nfp           ,   Nfp_raw                         )

  if (allocated(bn)) then
     HWRITERA( Nteta,Nzeta      ,   plas_Bn       ,   surf(1)%pb(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,        Bn       ,   surf(1)%bn(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bx            ,   surf(1)%Bx(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   By            ,   surf(1)%By(0:Nteta-1,0:Nzeta-1) )
     HWRITERA( Nteta,Nzeta      ,   Bz            ,   surf(1)%Bz(0:Nteta-1,0:Nzeta-1) )
  endif

  if (allocated(ppr)) then
     HWRITERA( pp_ns, pp_maxiter+1,   ppr         ,  ppr(1:pp_ns, 0:pp_maxiter) )
     HWRITERA( pp_ns, pp_maxiter+1,   ppz         ,  ppz(1:pp_ns, 0:pp_maxiter) )
     HWRITERV( pp_ns              ,   iota        ,  iota(1:pp_ns)              )
  endif

  if (allocated(XYZB)) then
     HWRITERC( total_num,4, pp_ns ,   XYZB        ,   XYZB(1:total_num, 1:4, 1:pp_ns) )
     HWRITERA( booz_mn,      pp_ns,  booz_mnc     ,   booz_mnc(1:booz_mn, 1:pp_ns)    )    
     HWRITERA( booz_mn,      pp_ns,  booz_mns     ,   booz_mns(1:booz_mn, 1:pp_ns)    )    
     HWRITEIV( booz_mn,              bmim         ,   bmim(1:booz_mn)                 )  
     HWRITEIV( booz_mn,              bmin         ,   bmin(1:booz_mn)                 )  
  endif

  HWRITERV( 1                ,  time_initialize,   time_initialize               )
  HWRITERV( 1                ,  time_optimize  ,   time_optimize                 )
  HWRITERV( 1                ,  time_postproc  ,   time_postproc                 )

  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )


  !--------------------------write coils.ext file-----------------------------------------------  

  if( save_coils == 1 ) then

     open(funit,file=trim(out_coils), status="unknown", form="formatted" )
     write(funit,'("periods "I3)') Nfp_raw
     write(funit,'("begin filament")')
     write(funit,'("mirror NIL")')
     do icoil = 1, Ncoils
        do ii = 0, coil(icoil)%NS-1
           write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), coil(icoil)%I
        enddo
        ii =  0
        write(funit,1010) coil(icoil)%xx(ii), coil(icoil)%yy(ii), coil(icoil)%zz(ii), &
             zero, icoil, coil(icoil)%name
     enddo
     write(funit,'("end")')
     close(funit)

  endif

1010 format(4es23.15,:,i9,"  ",a10)

  !--------------------------write .ext.fo.filaments.xxx file-----------------------------------

  write(srestart,'(i6.6)') iout

  if( save_filaments == 1 ) then

     open( funit, file="."//trim(ext)//".filaments."//srestart, status="unknown", form="unformatted" )  
     write(funit) Ncoils, Nseg
     do icoil = 1, Ncoils
        write(funit) coil(icoil)%xx(0:coil(icoil)%NS)
        write(funit) coil(icoil)%yy(0:coil(icoil)%NS)
        write(funit) coil(icoil)%zz(0:coil(icoil)%NS)
     enddo
     close( funit )

  endif


  !--------------------------write ext.harmonics file-----------------------------------

  if (save_harmonics == 1 .and. allocated(Bmnc)) then

     open(wunit, file=trim(out_harm), status='unknown', action='write')
     write(wunit,'("#NBmn")')                     ! comment line;
     write(wunit,'(I6)') NBmn                     ! write dimensions
     write(wunit,'("# n  m   Bmnc  Bmns  wBmn")') ! comment line;
     do imn = 1, NBmn
        write(wunit,'(2(I3, 4X), 3(ES23.15,4X))') Bmnin(imn)/Nfp_raw, Bmnim(imn), Bmnc(imn), Bmns(imn), wBmn(imn)
     enddo
     close(wunit)

  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if (update_plasma == 1 ) call write_plasma

  return

end subroutine saving


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_plasma
!-------------------------------------------------------------------------------!
! write down the unpdated plasma boundary information;                          !
! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
!-------------------------------------------------------------------------------!
  use globals, only : dp, zero, half, two, pi2, myid, ncpu, ounit, wunit, ext, &
                      Nfou, Nfp, NBnf, bim, bin, Bnim, Bnin, Rbc, Rbs, Zbc, Zbs, Bnc, Bns,  &
                      Nteta, Nzeta, surf, Nfp_raw, bnorm, sqrtmachprec, out_plasma
  
  implicit none  
  include "mpif.h"

  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! predefined Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, maxN, maxM
  REAL                :: teta, zeta, arg, tol, tmpc, tmps
  !------------------------------------------------------------------------------- 

  mf = 24 ;  nf = 24
  FATAL(bnftran, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

  tmpc = zero ; tmps = zero

  if (bnorm .gt. sqrtmachprec ) then
     tol = 1.0E-8 * bnorm
  else
     tol = 1.0E-8
  endif

  if(myid .ne. 0) return

  if(Nbnf .gt. 0) then  ! if there is input Bn target
     DALLOCATE(bnim)
     DALLOCATE(bnin)
     DALLOCATE(bnc )
     DALLOCATE(bns )
  endif

  Nbnf = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

  SALLOCATE( bnim, (1:Nbnf), 0    )
  SALLOCATE( bnin, (1:Nbnf), 0    )
  SALLOCATE( bnc , (1:Nbnf), zero )
  SALLOCATE( bns , (1:Nbnf), zero )  
  
  imn = 0
  do in = -nf, nf
     do im = 0, mf

        tmpc = zero ; tmps = zero
        do ii = 0, Nteta-1 
           teta = ( ii + half ) * pi2 / Nteta
           do jj = 0, Nzeta-1
              zeta = ( jj + half ) * pi2 / Nzeta

              arg = im*teta - in*(Nfp_raw/Nfp)*zeta
              tmpc = tmpc + surf(1)%bn(ii,jj)*cos(arg)
              tmps = tmps + surf(1)%bn(ii,jj)*sin(arg)

           enddo ! end jj
        enddo ! end ii

        if ( (abs(tmpc) + abs(tmps)) .lt. tol ) cycle

        imn = imn + 1
        bnin(imn) = in * (Nfp_raw/Nfp) ; bnim(imn) = im

        if (im .eq. 0  ) then
           tmpc = tmpc*half
           tmps = tmps*half
        endif
        bnc(imn) = tmpc
        bns(imn) = tmps

     enddo ! end im
  enddo ! end in

  Nbnf = imn

  bnc = bnc * two / (Nteta*Nzeta)
  bns = bns * two / (Nteta*Nzeta)
  !----------------------------------------------

  open(wunit, file=trim(out_plasma), status='unknown', action='write')

  write(wunit,*      ) "#Nfou Nfp  Nbnf"
  write(wunit,'(3I6)' ) Nfou, Nfp_raw, Nbnf

  write(wunit,*      ) "#------- plasma boundary------"
  write(wunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
  do imn = 1, Nfou
     write(wunit,'(2I6, 4ES15.6)') bin(imn)/Nfp_raw, bim(imn), Rbc(imn), Rbs(imn), Zbc(imn), Zbs(imn)
  enddo

  write(wunit,*      ) "#-------Bn harmonics----------"
  write(wunit,*      ) "#  n  m  bnc   bns"
  if (Nbnf .gt. 0) then
     do imn = 1, Nbnf
        write(wunit,'(2I6, 2ES15.6)') bnin(imn)/(Nfp_raw/nfp), bnim(imn), bnc(imn), bns(imn)
     enddo
  else
     write(wunit,'(2I6, 2ES15.6)') 0, 0, 0.0, 0.0
  endif

  close(wunit)
END SUBROUTINE write_plasma

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-


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


  INTEGER            :: ii, jj, icoil, NF, ip, is, cs, Npc

  ! the following are used by the macros HWRITEXX below; do not alter/remove;
  INTEGER            :: hdfier, rank
  integer(hid_t)     :: file_id, space_id, dset_id                 ! warning: the string "INTEGER" is a macro;
  integer(hsize_t)   :: onedims(1:1), twodims(1:2), threedims(1:3) ! warning: the string "INTEGER" is a macro;
  INTEGER            :: imn
  REAL               :: tvolume

  CHARACTER          :: srestart*6


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
     if(allocated(t1CU)) deriv(1:Ndof,7)=t1CU(1:Ndof)
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
  HWRITEIV( 1                ,   case_curv     ,   case_curv                     )
  HWRITEIV( 1                ,   curv_alpha    ,   curv_alpha                    )
  HWRITERV( 1                ,   weight_bnorm  ,   weight_bnorm                  )
  HWRITERV( 1                ,   weight_bharm  ,   weight_bharm                  )
  HWRITERV( 1                ,   weight_tflux  ,   weight_tflux                  )
  HWRITERV( 1                ,   target_tflux  ,   target_tflux                  )
  HWRITERV( 1                ,   weight_ttlen  ,   weight_ttlen                  )
  HWRITERV( 1                ,   target_length ,   target_length                 )
  HWRITERV( 1                ,   k0            ,   k0                            ) !Might need to fix this 
  HWRITERV( 1                ,   weight_specw  ,   weight_specw                  )
  HWRITERV( 1                ,   weight_cssep  ,   weight_cssep                  )
  HWRITERV( 1                ,   weight_gnorm  ,   weight_gnorm                  )
  HWRITERV( 1                ,   weight_inorm  ,   weight_inorm                  )
  HWRITERV( 1                ,   weight_mnorm  ,   weight_mnorm                  )
  HWRITERV( 1                ,   weight_curv   ,   weight_curv                   )
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

  call h5fclose_f( file_id, hdfier ) ! terminate access;
  FATAL( restart, hdfier.ne.0, error calling h5fclose_f )

  call h5close_f( hdfier ) ! close Fortran interface;
  FATAL( restart, hdfier.ne.0, error calling h5close_f )


  return

end subroutine saving


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

SUBROUTINE write_plasma
!-------------------------------------------------------------------------------!
! write down the unpdated plasma boundary information;                          !
! CZHU; first version: 2017/01/11; last revised: 2017/01/11                     !
!-------------------------------------------------------------------------------!
  use globals, only : dp, zero, half, two, pi2, myid, ncpu, ounit, wunit, ext, &
                      plasma, MPI_COMM_FOCUS, IsSymmetric, &
                      Nteta, Nzeta, surf, bnorm, sqrtmachprec, out_plasma
  use mpi
  implicit none  

  !-------------------------------------------------------------------------------
  INTEGER             :: mf, nf  ! predefined Fourier modes size
  INTEGER             :: imn=0, ii, jj, im, in, astat, ierr, maxN, maxM, isurf
  REAL                :: teta, zeta, arg, tol, tmpc, tmps
  !------------------------------------------------------------------------------- 

  ! use plasma as default
  isurf = plasma
  mf = 24 ;  nf = 24
  FATAL(bnftran, mf .le. 0 .and. nf .le. 0, INVALID size for Fourier harmonics)

  tmpc = zero ; tmps = zero

  if (bnorm .gt. sqrtmachprec ) then
     tol = 1.0E-8 * bnorm
  else
     tol = 1.0E-8
  endif

  if(myid .ne. 0) return
  FATAL( write_plasma, IsSymmetric==2, option not supported for now)

  if(surf(isurf)%Nbnf .gt. 0) then  ! if there is input Bn target
     DALLOCATE(surf(isurf)%bnim)
     DALLOCATE(surf(isurf)%bnin)
     DALLOCATE(surf(isurf)%bnc )
     DALLOCATE(surf(isurf)%bns )
  endif

  surf(isurf)%Nbnf = (mf+1)*(2*nf+1) ! (0:mf)*(-nf:nf)

  SALLOCATE( surf(isurf)%bnim, (1:surf(isurf)%Nbnf), 0    )
  SALLOCATE( surf(isurf)%bnin, (1:surf(isurf)%Nbnf), 0    )
  SALLOCATE( surf(isurf)%bnc , (1:surf(isurf)%Nbnf), zero )
  SALLOCATE( surf(isurf)%bns , (1:surf(isurf)%Nbnf), zero )  
  
  imn = 0
  do in = -nf, nf
     do im = 0, mf
        tmpc = zero ; tmps = zero
        do jj = 0, Nzeta-1
           zeta = ( jj + half ) * pi2 / surf(isurf)%Nzeta
           do ii = 0, Nteta-1
              teta = ( ii + half ) * pi2 / surf(isurf)%Nteta
              arg = im*teta - in*surf(isurf)%Nfp*zeta
              tmpc = tmpc + surf(isurf)%bn(ii,jj)*cos(arg)
              tmps = tmps + surf(isurf)%bn(ii,jj)*sin(arg)

           enddo ! end jj
        enddo ! end ii

        if ( (abs(tmpc) + abs(tmps)) .lt. tol ) cycle

        imn = imn + 1
        surf(isurf)%bnin(imn) = in * surf(isurf)%Nfp
        surf(isurf)%bnim(imn) = im

        if (im .eq. 0  ) then
           tmpc = tmpc*half
           tmps = tmps*half
        endif
        surf(isurf)%bnc(imn) = tmpc
        surf(isurf)%bns(imn) = tmps

     enddo ! end im
  enddo ! end in

  surf(isurf)%Nbnf = imn

  surf(isurf)%bnc = surf(isurf)%bnc * two / (Nteta*Nzeta)
  surf(isurf)%bns = surf(isurf)%bns * two / (Nteta*Nzeta)
  !----------------------------------------------


  open(wunit, file=trim(out_plasma), status='unknown', action='write')

  write(wunit,*      ) "#Nfou Nfp  Nbnf"
  write(wunit,'(3I6)' ) surf(isurf)%Nfou, surf(isurf)%Nfp, surf(isurf)%Nbnf

  write(wunit,*      ) "#------- plasma boundary------"
  write(wunit,*      ) "#  n   m   Rbc   Rbs    Zbc   Zbs"
  do imn = 1, surf(isurf)%Nfou
     write(wunit,'(2I6, 4ES15.6)') surf(isurf)%bin(imn)/surf(isurf)%Nfp, surf(isurf)%bim(imn), & 
          surf(isurf)%Rbc(imn), surf(isurf)%Rbs(imn), surf(isurf)%Zbc(imn), surf(isurf)%Zbs(imn)
  enddo

  write(wunit,*      ) "#-------Bn harmonics----------"
  write(wunit,*      ) "#  n  m  bnc   bns"
  if (surf(isurf)%Nbnf .gt. 0) then
     do imn = 1, surf(isurf)%Nbnf
        write(wunit,'(2I6, 2ES15.6)') surf(isurf)%bnin(imn)/surf(isurf)%Nfp, surf(isurf)%bnim(imn), &
             surf(isurf)%bnc(imn), surf(isurf)%bns(imn)
     enddo
  else
     write(wunit,'(2I6, 2ES15.6)') 0, 0, 0.0, 0.0
  endif

  close(wunit)
END SUBROUTINE write_plasma

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-

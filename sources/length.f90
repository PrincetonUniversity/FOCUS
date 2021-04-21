!title (tlength) ! Calculate total length cost functon and its derivatives. (czhu)

!latex \briefly{The constraint on coil length can prevent coils to be too long. There are two different
!latex         forms of length constraint, exponential and quadratic, controlled by \inputvar{case\_length}.
!latex         Exponential form is used to shorten coil length as much as possible, while quadratic form 
!latex         is forcing coils to have a length of the \emph{targt\_length}.}

!latex \calledby{\link{solvers}}

!latex  \section{General}
!latex  Without a constraint on length, the coils can become arbitrarily long to lower the ripple, 
!latex  and more ``wiggles" can potentially be formed to better match the plasma shape.
!latex  Besides, the total coil length is directly related to the usage of materials, 
!latex  i.e. cost, which is without any relation to the physical requirements.
!latex  We include an objective function of in exponential form
!latex  \be
!latex  \ds f_L(\vect{X}) = \frac{1}{N_C} \; \sum_{i=1}^{N_C} \frac{e^{L_i}}{e^{L_{i,o}}},
!latex  \ee
!latex  or in quadratic form
!latex  \be
!latex  \ds f_L(\vect{X}) = \frac{1}{N_C} \; \sum_{i=1}^{N_C} \half \frac{(L_i - L{i,o})^2}{L_{i,o}^2},
!latex  \ee
!latex  where $L_i({\bf X})$ is the length of $i$-th coil,
!latex  \be
!latex  \ds L_i(\vect{X}) = \int_0^{2\pi} \!\! |\vect{x}_i'| \; \dd{t},
!latex  \ee 
!latex  and $L_{i,0}$ is a user-specified normalization.
!latex  The variations in $L$ resulting from a variation $\delta \vect{x}_i$ is
!latex  \be\label{eq:var_L}
!latex  \ds \delta L(\vect{X})  = \int_0^{2\pi} \frac{ ({\vect{x}'_i} \cdot {\vect{x}''_i}){\vect{x}'_i}
!latex  - ({\vect{x}'_i} \cdot {\vect{x}'_i}){\vect{x}''_i} }
!latex  {({\vect{x}_i'} \cdot {\vect{x}_i'})^{3/2}} \cdot \delta \vect{x}_i \ \dd{t}. 
!latex  \ee
!latex  
!latex  \section{First derivatives}
!latex  From \Eqn{var_L}, we can calculated the functional derivatives of coil length with respect to the 
!latex  coil geometries,
!latex  \begin{align}
!latex  \ds \frac{\delta L}{\delta x} & =  \frac{y'y''x' + z'z''x' - y'y'x'' - z'z'x''}{(x'x' + y'y' + z'z')^{3/2}} \ ; \\
!latex  \ds \frac{\delta L}{\delta y} & =  \frac{x'x''y' + z'z''y' - x'x'y'' - z'z'y''}
!latex                                          {(x'x' + y'y' + z'z')^{3/2}} \ ; \\
!latex  \ds \frac{\delta L}{\delta z} & =  \frac{x'x''z' + y'y''z' - x'x'z'' - y'y'z''}
!latex                                          {(x'x' + y'y' + z'z')^{3/2}} \ .
!latex  \end{align}
!latex  Then the derivatives of coil length with respect to coil parameters are
!latex  \be
!latex  \ds \pdv{L}{X_i} = \int_0^{2\pi} \frac{\delta L}{\delta x} \pdv{x}{X_i} +
!Latex                                   \frac{\delta L}{\delta y} \pdv{y}{X_i} +
!Latex                                   \frac{\delta L}{\delta z} \pdv{z}{X_i}  \dd{t} \ .
!latex  \ee
!latex  Here, $X_i$ is an arbitrary variable $X_i \in \vect{X} = \{I^1, X_{c,0}^1, \dots \}$.
!latex  So the first derivatives of coil length objective function are
!latex  \begin{align}
!latex  \ds \pdv{f_L}{X_i} & = \frac{1}{N_C} \frac{e^{L_i}}{e^{L_{i,o}}} \pdv{L}{X_i} \  ; \\
!latex  \ds \pdv{f_L}{X_i} & = \frac{1}{N_C} \frac{(Li - L{i,o})}{L_{i,o}^2} \  \pdv{L}{X_i} \ .
!latex  \end{align}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! not parallelized; communications may take more time;
subroutine length(ideriv)
  use globals, only: dp, zero, half, pi2, machprec, ncpu, myid, ounit, &
       coil, DoF, Ncoils, Nfixgeo, Ndof, ttlen, t1L, t2L, case_length, &
       ittlen, mttlen, LM_fvec, LM_fjac, weight_ttlen, length_delta, &
       MPI_COMM_FOCUS

  use mpi
  implicit none
  INTEGER, INTENT(in) :: ideriv

  INTEGER             :: astat, ierr, icoil, idof, ND, ivec
  REAL                :: d1L(1:Ndof), norm(1:Ncoils)

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  ttlen = zero
  ivec = 1

  if( ideriv >= 0 ) then
     do icoil = 1, Ncoils     !only care about unique coils;
        if(coil(icoil)%type == 1) then  ! only for Fourier
           !if( myid.ne.modulo(icoil-1,ncpu) ) cycle ! parallelization loop;
           call LenDeriv0(icoil, coil(icoil)%L)
           !RlBCAST( coil(icoil)%L, 1, modulo(icoil-1,ncpu) ) !broadcast each coil's length
           if ( coil(icoil)%Lc /= 0 ) then
              if (case_length == 1) then ! quadratic;
                 ttlen = ttlen +  half * (coil(icoil)%L - coil(icoil)%Lo)**2 / coil(icoil)%Lo**2
                 if (mttlen > 0) then ! L-M format of targets
                    LM_fvec(ittlen+ivec) = weight_ttlen * (coil(icoil)%L - coil(icoil)%Lo)
                    ivec = ivec + 1
                 endif
              elseif (case_length == 2) then ! exponential;
                 ttlen = ttlen + exp(coil(icoil)%L) / exp(coil(icoil)%Lo)
                 if (mttlen > 0) then ! L-M format of targets
                    LM_fvec(ittlen+ivec) = weight_ttlen * exp(coil(icoil)%L) / exp(coil(icoil)%Lo)
                    ivec = ivec + 1
                 endif
              elseif (case_length == 3) then ! Delta quadratic 
                 if (length_delta < 0.0) then
                    FATAL( length, .true. , length_delta cannot be negative )
                 endif
                 if (abs(coil(icoil)%L - coil(icoil)%Lo) .ge.  length_delta) then
                    ttlen = ttlen + half * (abs(coil(icoil)%L - coil(icoil)%Lo) - length_delta)**2 / coil(icoil)%Lo**2
                 endif
              else
                 FATAL( length, .true. , invalid case_length option )
              endif  
           endif
        endif
     enddo
     if (mttlen > 0) then ! L-M format of targets
        FATAL( length, ivec == mttlen, Errors in counting ivec for L-M )
     endif
     ttlen = ttlen / (Ncoils - Nfixgeo + machprec)
  endif

  !-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

  if ( ideriv >= 1 ) then
     t1L = zero ; d1L = zero ; norm = zero
     idof = 0 ; ivec = 1
     do icoil = 1, Ncoils
        ND = DoF(icoil)%ND
        if ( coil(icoil)%Ic /= 0 ) then !if current is free;
           idof = idof +1
        endif
        if ( coil(icoil)%Lc /= 0 ) then !if geometry is free;
           if(coil(icoil)%type .eq. 1) then ! only for Fourier
              ! calculate normalization
              if (case_length == 1) then
                 norm(icoil) = (coil(icoil)%L - coil(icoil)%Lo) / coil(icoil)%Lo**2  ! quadratic;
              elseif (case_length == 2) then
                 norm(icoil) = exp(coil(icoil)%L) / exp(coil(icoil)%Lo)       ! exponential;
              elseif (case_length == 3) then
                 if (coil(icoil)%L < coil(icoil)%Lo - length_delta) then      ! Delta quadratic
                    norm(icoil) = (coil(icoil)%L - coil(icoil)%Lo + length_delta) / coil(icoil)%Lo**2
                 elseif (coil(icoil)%L > coil(icoil)%Lo + length_delta) then 
                    norm(icoil) = (coil(icoil)%L - coil(icoil)%Lo - length_delta) / coil(icoil)%Lo**2
                 endif
              else
                 FATAL( length, .true. , invalid case_length option )
              endif  
              ! call lederiv1 to calculate the 1st derivatives
              call lenDeriv1( icoil, d1L(idof+1:idof+ND), ND )
              t1L(idof+1:idof+ND) = d1L(idof+1:idof+ND) * norm(icoil)
              if (mttlen > 0) then ! L-M format of targets
                 LM_fjac(ittlen+ivec, idof+1:idof+ND) = weight_ttlen * d1L(idof+1:idof+ND)
                 if (case_length == 2) &
                      & LM_fjac(ivec, idof+1:idof+ND) = LM_fjac(ivec, idof+1:idof+ND) &
                      & * exp(coil(icoil)%L) / exp(coil(icoil)%Lo)
                 ! Need to add in case_length == 3 
                 ivec = ivec + 1
              endif
           endif 
           idof = idof + ND
        endif

     enddo !end icoil;
     FATAL( length , idof .ne. Ndof, counting error in packing )

     if (mttlen > 0) then ! L-M format of targets
        FATAL( length, ivec == mttlen, Errors in counting ivec for L-M )
     endif

     t1L = t1L / (Ncoils - Nfixgeo + machprec)
  endif

  return
end subroutine length

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

! calculate coil length

subroutine LenDeriv0(icoil, length)

  use globals, only: dp, zero, coil, myid, ounit, Ncoils, MPI_COMM_FOCUS 
  use mpi
  implicit none

  INTEGER, intent(in)  :: icoil
  REAL   , intent(out) :: length
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: dlength

  FATAL( LenDeriv0, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  length = zero
  
  do kseg = 0, coil(icoil)%NS-1

     dlength = sqrt(coil(icoil)%xt(kseg)**2 + coil(icoil)%yt(kseg)**2 + coil(icoil)%zt(kseg)**2)
     length  = length + dlength * coil(icoil)%dd(kseg)

  enddo ! end kseg
 
  return

end subroutine LenDeriv0

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!calculate coil length in derivs(0) and first derivatives in derivs(1:Cdof)

subroutine LenDeriv1(icoil, derivs, ND)

  use globals, only: dp, zero, pi2, coil, DoF, myid, ounit, Ncoils, MPI_COMM_FOCUS
  use mpi
  implicit none

  INTEGER, intent(in)  :: icoil, ND
  REAL   , intent(out) :: derivs(1:1, 1:ND)
!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!
  INTEGER              :: kseg, astat, ierr
  REAL                 :: dl3, xt, yt, zt, xa, ya, za
  REAL, dimension(1:1, 0:coil(icoil)%NS-1) :: dLx, dLy, dLz

  FATAL( LenDeriv1, icoil .lt. 1 .or. icoil .gt. Ncoils, icoil not in right range )
  
  derivs = zero
  
  do kseg = 0, coil(icoil)%NS-1
     
     xt = coil(icoil)%xt(kseg) ; yt = coil(icoil)%yt(kseg) ; zt = coil(icoil)%zt(kseg)
     xa = coil(icoil)%xa(kseg) ; ya = coil(icoil)%ya(kseg) ; za = coil(icoil)%za(kseg)

     dl3 = sqrt(xt*xt + yt*yt + zt*zt)**3

     dLx(1,kseg) = ( yt*ya*xt + zt*za*xt - yt*yt*xa - zt*zt*xa ) / dl3 * coil(icoil)%dd(kseg)
     dLy(1,kseg) = ( xt*xa*yt + zt*za*yt - xt*xt*ya - zt*zt*ya ) / dl3 * coil(icoil)%dd(kseg)
     dLz(1,kseg) = ( xt*xa*zt + yt*ya*zt - xt*xt*za - yt*yt*za ) / dl3 * coil(icoil)%dd(kseg)

  enddo ! end kseg

  derivs(1:1, 1:ND) = matmul(dLx, DoF(icoil)%xof) + matmul(dLy, DoF(icoil)%yof) + matmul(dLz, DoF(icoil)%zof)

  return

end subroutine LenDeriv1


!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

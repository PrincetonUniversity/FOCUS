!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!title (numerics) ! Some miscellaneous numerical routines.

!latex \briefly{briefly}
!latex \calledby{\link{}}
!latex \calls{\link{}}
!latex \tableofcontents

!latex \begin{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item \verb+tfft+

subroutine tfft( Nt, Nz, ijreal, ijimag, isr, trigm, trign, trigwk, mn, im, in, efmn, ofmn, cfmn, sfmn, ifail )

 !use constants
 !use inputlist, only : Nfp
 !use allglobal, only : pi2nfp

  use kmodule, only : zero, one, pi2

  implicit none
  
  INTEGER   :: Nt, Nz, mn, im(mn), in(mn), Ntz, imn, ifail, mm, nn, NFP
  REAL      :: ijreal(1:Nt*Nz), ijimag(1:Nt*Nz), trigm(2*Nt), trign(2*Nz), trigwk(2*Nt*Nz), efmn(1:mn), ofmn(1:mn), cfmn(1:mn), sfmn(1:mn)
  CHARACTER :: isr
  
  LOGICAL   :: check=.false.
  INTEGER   :: jj, kk
  REAL      :: jireal(1:Nt*Nz), jiimag(1:Nt*Nz), arg, ca, sa

  NFP = 1

  if( check ) then ; jireal=ijreal ; jiimag=ijimag
  endif

  Ntz=Nt*Nz

  call C06FUF( Nt , Nz , ijreal(1:Ntz) , ijimag(1:Ntz) , isr , trigm , trign , trigwk , ifail )
  
  ijreal(:) = ijreal(:)/sqrt(one*Ntz) ; ijimag(:) = ijimag(:)/sqrt(one*Ntz)
  
  cfmn=zero ; sfmn=zero ; efmn=zero ; ofmn=zero
  
  do imn = 1,mn ; mm = im(imn) ; nn = in(imn) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)+(   nn)*Nt) + ijreal(1+(   mm)+(Nz-nn)*Nt)
    ofmn(imn) =   ijimag(1+(Nt-mm)+(   nn)*Nt) - ijimag(1+(   mm)+(Nz-nn)*Nt)
    cfmn(imn) =   ijimag(1+(Nt-mm)+(   nn)*Nt) + ijimag(1+(   mm)+(Nz-nn)*Nt)
    sfmn(imn) = - ijreal(1+(Nt-mm)+(   nn)*Nt) + ijreal(1+(   mm)+(Nz-nn)*Nt) 
    
   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)           ) + ijreal(1+(   mm)           )
    ofmn(imn) =   ijimag(1+(Nt-mm)           ) - ijimag(1+(   mm)           )
    cfmn(imn) =   ijimag(1+(Nt-mm)           ) + ijimag(1+(   mm)           )
    sfmn(imn) = - ijreal(1+(Nt-mm)           ) + ijreal(1+(   mm)           )
    
   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
    efmn(imn) =   ijreal(1+(Nt-mm)+(Nz+nn)*Nt) + ijreal(1+(   mm)+(  -nn)*Nt)
    ofmn(imn) =   ijimag(1+(Nt-mm)+(Nz+nn)*Nt) - ijimag(1+(   mm)+(  -nn)*Nt)
    cfmn(imn) =   ijimag(1+(Nt-mm)+(Nz+nn)*Nt) + ijimag(1+(   mm)+(  -nn)*Nt)
    sfmn(imn) = - ijreal(1+(Nt-mm)+(Nz+nn)*Nt) + ijreal(1+(   mm)+(  -nn)*Nt) 
    
   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
    efmn(imn) =   ijreal(1+        (Nz-nn)*Nt) + ijreal(1+        (   nn)*Nt)
    ofmn(imn) = - ijimag(1+        (Nz-nn)*Nt) + ijimag(1+        (   nn)*Nt)
    cfmn(imn) =   ijimag(1+        (Nz-nn)*Nt) + ijimag(1+        (   nn)*Nt)
    sfmn(imn) =   ijreal(1+        (Nz-nn)*Nt) - ijreal(1+        (   nn)*Nt)
    
   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    efmn(imn) =   ijreal(1)
    ofmn(imn) =     zero
    cfmn(imn) =   ijimag(1)
    sfmn(imn) =     zero
    
   endif
   
  enddo

  if( .not.check ) return
  
  ijreal(1:Ntz)=zero ; ijimag(1:Ntz)=zero
  
  do jj=0,Nt-1
   do kk=0,Nz-1
    do imn=1,mn ; arg=im(imn)*jj*pi2/Nt-in(imn)*kk*(pi2/nfp)/Nz ; ca=cos(arg) ; sa=sin(arg)
     
     ijreal(1+jj+kk*Nt) = ijreal(1+jj+kk*Nt) + efmn(imn)*ca + ofmn(imn)*sa
     ijimag(1+jj+kk*Nt) = ijimag(1+jj+kk*Nt) + cfmn(imn)*ca + sfmn(imn)*sa
     
    enddo
   enddo
  enddo
  
  write(*,'("tfft : Fourier reconstruction error = "2es15.5)')sqrt(sum((ijreal-jireal)**2)/Ntz),sqrt(sum((ijimag-jiimag)**2)/Ntz)

  return

end subroutine tfft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \item \verb+invfft+

subroutine invfft( mn , im , in , efmn , ofmn , cfmn , sfmn , Nt , Nz , ijreal , ijimag , isr , trigm , trign , trigwk )

 !use constants, only : zero, half, one
 !use inputlist, only : Nfp

  use kmodule, only : zero, half, one

  implicit none
  INTEGER  , intent(in)    :: mn, im(mn), in(mn)
  REAL     , intent(in)    :: efmn(mn), ofmn(mn), cfmn(mn), sfmn(mn)
  INTEGER  , intent(in)    :: Nt, Nz
  REAL     , intent(out)   :: ijreal(Nt*Nz), ijimag(Nt*Nz) ! output real space;
  CHARACTER, intent(inout) :: isr
  REAL     , intent(inout) :: trigm(2*Nt), trign(2*Nz), trigwk(2*Nt*Nz)
  
  INTEGER                   :: Ntz, imn, jj, kk, c06fuffail, c06gcffail, mm, nn, NFP
  
  NFP = 1

  Ntz = Nt*Nz ; ijreal(1:Ntz) = zero ; ijimag(1:Ntz) = zero

  do imn = 1,mn ; mm = im(imn) ; nn = in(imn) / Nfp
   
   if    ( mm.gt.0 .and. nn.gt.0 ) then
    
    ijreal(1+(Nt-mm)+(   nn)*Nt) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)+(Nz-nn)*Nt) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)+(   nn)*Nt) = (   ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)+(Nz-nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.gt.0 .and. nn.eq.0 ) then
    
    ijreal(1+(Nt-mm)           ) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)           ) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)           ) = ( + ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)           ) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.gt.0 .and. nn.lt.0 ) then
    
    ijreal(1+(Nt-mm)+(Nz+nn)*Nt) = (   efmn(imn) - sfmn(imn) ) * half
    ijreal(1+(   mm)+(  -nn)*Nt) = (   efmn(imn) + sfmn(imn) ) * half

    ijimag(1+(Nt-mm)+(Nz+nn)*Nt) = ( + ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+(   mm)+(  -nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.eq.0 .and. nn.gt.0 ) then
    
    ijreal(1+        (Nz-nn)*Nt) = ( + efmn(imn) + sfmn(imn) ) * half
    ijreal(1+        (   nn)*Nt) = ( + efmn(imn) - sfmn(imn) ) * half
    
    ijimag(1+        (Nz-nn)*Nt) = ( - ofmn(imn) + cfmn(imn) ) * half
    ijimag(1+        (   nn)*Nt) = ( + ofmn(imn) + cfmn(imn) ) * half

   elseif( mm.eq.0 .and. nn.eq.0 ) then
    
    ijreal(1) = efmn(imn)
    ijimag(1) = cfmn(imn)
    
   endif
   
  enddo

  ijreal(1:Ntz) = ijreal(1:Ntz) * sqrt(one*Ntz)
  ijimag(1:Ntz) = ijimag(1:Ntz) * sqrt(one*Ntz)
  
  c06gcffail=0 ; call C06GCF( ijimag, Ntz , c06gcffail )
  c06fuffail=0 ; call C06FUF( Nt , Nz , ijreal , ijimag , isr , trigm , trign , trigwk , c06fuffail )
  c06gcffail=0 ; call C06GCF( ijimag , Ntz , c06gcffail )
  
  return

end subroutine invfft

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

!latex \end{itemize}

!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!

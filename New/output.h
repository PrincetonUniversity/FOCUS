
subroutine output (mark)

  use globals, only : zero, ounit, myid, ierr, astat, iout, Nouts, Ncoils, save_freq, Tdof, &
       coil, coilspace, FouCoil, chi, t1E, Bdotnsquared, bharm, toroidalfluxerror, ttlen, specw, ccsep, evolution, xdof, DoF, &
       exit_tol, exit_signal

  implicit none  
  include "mpif.h"

  REAL, INTENT( IN ) :: mark

  INTEGER            :: idof, NF, icoil
  REAL               :: sumdE


  iout = iout + 1
  
  FATAL( output , iout > Nouts+2, maximum iteration reached )

  sumdE = sqrt(sum(t1E**2)) ! Eucliean norm 2; 

  if (myid == 0) write(ounit, '("output  : "I6" : "9(ES12.5," ; "))') iout, mark, chi, sumdE, Bdotnsquared, bharm, &
       toroidalfluxerror, ttlen, specw, ccsep

  ! save evolution data;
  if (allocated(evolution)) then
     evolution(iout,0) = mark
     evolution(iout,1) = chi
     evolution(iout,2) = sumdE
     evolution(iout,3) = Bdotnsquared
     evolution(iout,4) = bharm
     evolution(iout,5) = toroidalfluxerror
     evolution(iout,6) = ttlen
     evolution(iout,7) = specw
     evolution(iout,8) = ccsep
  endif

  ! exit the optimization if no obvious changes in past 5 outputs; 07/20/2017
  if (iout>5) then
     if ( abs(evolution(iout,1) - evolution(iout-5, 1)) / evolution(iout,1) < exit_tol ) exit_signal = .True.
  end if
  
  !save all the coil parameters;
  if (allocated(coilspace)) then
     idof = 0
     do icoil = 1, Ncoils
        coilspace(iout, idof+1 ) = coil(icoil)%I ;  idof = idof + 1

        select case (coil(icoil)%itype)
        case (1)
           NF = FouCoil(icoil)%NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%xc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%xs(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%yc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%ys(1:NF) ; idof = idof + NF
           coilspace(iout, idof+1:idof+NF+1) = FouCoil(icoil)%zc(0:NF) ; idof = idof + NF +1
           coilspace(iout, idof+1:idof+NF  ) = FouCoil(icoil)%zs(1:NF) ; idof = idof + NF
        case default
           FATAL(solvers, .true., not supported coil types)
        end select
     enddo
     FATAL( output , idof .ne. Tdof, counting error in restart )
  endif

  if(mod(iout,save_freq) .eq. 0) call archive

  call MPI_BARRIER( MPI_COMM_WORLD, ierr ) ! wait all cpus;

  return  

end subroutine output


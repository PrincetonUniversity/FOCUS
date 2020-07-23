SUBROUTINE eval_spline_basis(icoil)
! Compute the basis functions of order 0,1,2 and 3 for the values of t contained in eval_points and store the value in Splines
    use globals, only : dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet,coil,Splines,astat  
    implicit none
    integer :: N,i,j,icoil
    REAL,allocatable :: eval_points(:)
    integer :: icp,ipos

    SALLOCATE(eval_points, (0:coil(icoil)%NS-1),zero)

    N = coil(icoil)%NS
    eval_points = Splines(icoil)%eval_points

    Splines(icoil)%basis_0 = 0
    Splines(icoil)%basis_1 = 0 
    Splines(icoil)%basis_2 = 0
    Splines(icoil)%basis_3 = 0

    do ipos=0,N-1 
    	do icp=0,Splines(icoil)%NCP+2
            if(eval_points(ipos)>=Splines(icoil)%vect(icp).AND.eval_points(ipos)<Splines(icoil)%vect(icp+1)) Splines(icoil)%basis_0(ipos,icp) = 1
        enddo
    enddo


    do icp=0,Splines(icoil)%NCP+1   
        do ipos=0,N -1      
            Splines(icoil)%basis_1(ipos,icp) = 1.0*(eval_points(ipos) - Splines(icoil)%vect(icp)) /        &
                (Splines(icoil)%vect(icp+1) - Splines(icoil)%vect(icp))  *Splines(icoil)%basis_0(ipos,icp) &
                                             + 1.0*( - eval_points(ipos) + Splines(icoil)%vect(icp+2))/    &
                (Splines(icoil)%vect(icp+2) - Splines(icoil)%vect(icp+1))*Splines(icoil)%basis_0(ipos,icp+1)
        enddo    
    enddo

    do icp=0,Splines(icoil)%NCP
        do ipos=0,N -1               
            Splines(icoil)%basis_2(ipos,icp) = 1.0*(   eval_points(ipos) - Splines(icoil)%vect(icp))/      & 
                (Splines(icoil)%vect(icp+2) - Splines(icoil)%vect(icp))*Splines(icoil)%basis_1(ipos,icp)   &
                                             + 1.0*( - eval_points(ipos) + Splines(icoil)%vect(icp+3))/    &
                (Splines(icoil)%vect(icp+3) - Splines(icoil)%vect(icp+1))*Splines(icoil)%basis_1(ipos,icp+1)
        enddo    
    enddo
    
    do icp=0,Splines(icoil)%NCP-1
        do ipos=0,N-1             
            Splines(icoil)%basis_3(ipos,icp) = 1.0*(   eval_points(ipos) - Splines(icoil)%vect(icp))/    &
                (Splines(icoil)%vect(icp+3) - Splines(icoil)%vect(icp))*Splines(icoil)%basis_2(ipos,icp) &
                                             + 1.0*( - eval_points(ipos) + Splines(icoil)%vect(icp+4))/  &
                (Splines(icoil)%vect(icp+4) - Splines(icoil)%vect(icp+1))*Splines(icoil)%basis_2(ipos,icp+1)
        enddo
    enddo

    return 
END SUBROUTINE eval_spline_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eval_spline_basis1(icoil)
! Compute the first derivatives of the third order basis functions in each points of eval_points
    use globals, only : dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet,coil,Splines,astat   
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: icp,ipos
    REAL :: vect(0:Splines(icoil)%NT-1)
       
    vect = Splines(icoil)%vect
    N = coil(icoil)%NS
    eval_points = Splines(icoil)%eval_points
    Splines(icoil)%db_dt=0

    do icp=0,Splines(icoil)%NCP-1
        do ipos=0,N-1 
            if(eval_points(ipos)>=vect(icp).AND.eval_points(ipos)<vect(icp+1)) then
                Splines(icoil)%db_dt(ipos,icp) = &
                     Splines(icoil)%basis_2(ipos,icp  )/(vect(icp+3)-vect(icp  )) + (eval_points(ipos)-vect(icp  ))/( vect(icp+3)-vect(icp  ))* &
                    (Splines(icoil)%basis_1(ipos,icp  )/(vect(icp+2)-vect(icp  )) + (eval_points(ipos)-vect(icp  ))/((vect(icp+2)-vect(icp  ))* &
		    (vect(icp+1)-vect(icp))))
						
            elseif(eval_points(ipos)>=vect(icp+1).AND.eval_points(ipos)<vect(icp+2)) then 
                Splines(icoil)%db_dt(ipos,icp) = &
                     Splines(icoil)%basis_2(ipos,icp  )/(vect(icp+3)-vect(icp  )) + (eval_points(ipos)-vect(icp  ))/( vect(icp+3)-vect(icp  ))* &
                    (Splines(icoil)%basis_1(ipos,icp  )/(vect(icp+2)-vect(icp  )) - (eval_points(ipos)-vect(icp  ))/((vect(icp+2)-vect(icp  ))* &
                    (vect(icp+2)-vect(icp+1)))                                                                                                  & 
                   - Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1)) + (vect(icp+3)-eval_points(ipos))/((vect(icp+3)-vect(icp+1))* &
                    (vect(icp+2)-vect(icp+1))))                                                                                                 &
                   - Splines(icoil)%basis_2(ipos,icp+1)/(vect(icp+4)-vect(icp+1)) + (vect(icp+4)-eval_points(ipos))/( vect(icp+4)-vect(icp+1))* &
                    (Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1)) + (eval_points(ipos)-vect(icp+1))/((vect(icp+3)-vect(icp+1))* &
                    (vect(icp+2)-vect(icp+1))))
						
            elseif(eval_points(ipos)>=vect(icp+2).AND.eval_points(ipos)<vect(icp+3)) then 
                Splines(icoil)%db_dt(ipos,icp) =  &
                     Splines(icoil)%basis_2(ipos,icp  )/(vect(icp+3)-vect(icp  )) + (eval_points(ipos)-vect(icp  ))/( vect(icp+3)-vect(icp  ))* &
                   (-Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1)) - (vect(icp+3)-eval_points(ipos))/((vect(icp+3)-vect(icp+1))* & 
                   (vect(icp+3)-vect(icp+2))))                                                                                                  &
                   - Splines(icoil)%basis_2(ipos,icp+1)/(vect(icp+4)-vect(icp+1)) + (vect(icp+4)-eval_points(ipos))/( vect(icp+4)-vect(icp+1))* &
                    (Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1)) - (eval_points(ipos)-vect(icp+1))/((vect(icp+3)-vect(icp+1))* &
                   (vect(icp+3)-vect(icp+2)))                                                                                                   &
                   - Splines(icoil)%basis_1(ipos,icp+2)/(vect(icp+4)-vect(icp+2)) + (vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+2))* &
                   (vect(icp+3)-vect(icp+2))))
						
            elseif(eval_points(ipos)>=vect(icp+3).AND.eval_points(ipos)<vect(icp+4)) then 
                Splines(icoil)%db_dt(ipos,icp) = &
                    -Splines(icoil)%basis_2(ipos,icp+1)/(vect(icp+4)-vect(icp+1)) + (vect(icp+4)-eval_points(ipos))/( vect(icp+4)-vect(icp+1))* &
                   (-Splines(icoil)%basis_1(ipos,icp+2)/(vect(icp+4)-vect(icp+2)) - (vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+2))* &
                   (vect(icp+4)-vect(icp+3))))
            
            endif
        enddo
    enddo


    return 
END SUBROUTINE eval_spline_basis1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE eval_spline_basis2(icoil)
! Compute the second derivatives of the third order basis functions in each points of eval_points
    use globals, only : dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet,coil,Splines,astat    
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: icp,ipos
    REAL :: vect(0:Splines(icoil)%NT-1)

    Splines(icoil)%db_dt_2 = 0
    vect = Splines(icoil)%vect
    N = coil(icoil)%NS
    eval_points = Splines(icoil)%eval_points
    
    do icp=0,Splines(icoil)%NCP-1
        do ipos=0,N-1 
            if(eval_points(ipos)>=vect(icp).AND.eval_points(ipos)<vect(icp+1)) then
                Splines(icoil)%db_dt_2(ipos,icp) = &
                    2.0/(vect(icp+3)-vect(icp)) * (Splines(icoil)%basis_1(ipos,icp)/(vect(icp+2)-vect(icp))             + &
	                (eval_points(ipos)-vect(icp))/((vect(icp+2)-vect(icp)) * (vect(icp+1)-vect(icp))))              + &
                    2.0*(eval_points(ipos)-vect(icp))/((vect(icp+3)-vect(icp)) * (vect(icp+1)-vect(icp))*(vect(icp+2)-vect(icp)))
						
            else if(eval_points(ipos)>=vect(icp+1).AND.eval_points(ipos)<vect(icp+2)) then 
                Splines(icoil)%db_dt_2(ipos,icp) = &
                    2.0/(vect(icp+3)-vect(icp  ))      *( Splines(icoil)%basis_1(ipos,icp  )/(vect(icp+2)-vect(icp  ))  - &
                        (eval_points(ipos)-vect(icp  ))/((vect(icp+2)-vect(icp  )) * (vect(icp+2)-vect(icp+1)))         - &
		                                          Splines(icoil)%basis_1(ipos,icp+1) /( vect(icp+3)-vect(icp+1))+ & 
		        (vect(icp+3)-eval_points(ipos))/((vect(icp+3)-vect(icp+1)) * (vect(icp+2)-vect(icp+1))))        - &
                        (eval_points(ipos)-vect(icp  ))/((vect(icp+3)-vect(icp  )) * (vect(icp+2)-vect(icp+1)))         * &
                   (2.0/(vect(icp+2)-vect(icp))        + 2.0/( vect(icp+3)-vect(icp+1)))                                - &
                    2.0/(vect(icp+4)-vect(icp+1))      *( Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1))  + &
                        (eval_points(ipos)-vect(icp+1))/((vect(icp+3)-vect(icp+1)) * (vect(icp+2)-vect(icp+1))))        + &
                    2.0*(vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+1)) * (vect(icp+3)-vect(icp+1)) * (vect(icp+2)-vect(icp+1))) 

            else if(eval_points(ipos)>=vect(icp+2).AND.eval_points(ipos)<vect(icp+3)) then
                Splines(icoil)%db_dt_2(ipos,icp) = &
                    2.0/(vect(icp+3)-vect(icp))        *(-Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1))                  - &
                        (vect(icp+3)-eval_points(ipos))/((vect(icp+3)-vect(icp+1))*(vect(icp+3)-vect(icp+2))))                          + &
                    2.0*(eval_points(ipos)-vect(icp  ))/((vect(icp+3)-vect(icp))  *(vect(icp+3)-vect(icp+1))*(vect(icp+3)-vect(icp+2))) - &
                    2.0/(vect(icp+4)-vect(icp+1))      *( Splines(icoil)%basis_1(ipos,icp+1)/(vect(icp+3)-vect(icp+1))                  - &
                        (eval_points(ipos)-vect(icp+1))/((vect(icp+3)-vect(icp+1))*(vect(icp+3)-vect(icp+2)))                           - &
                                                          Splines(icoil)%basis_1(ipos,icp+2)/( vect(icp+4)-vect(icp+2))                 + &
	                (vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+2))*(vect(icp+3)-vect(icp+2))))                          - &
                        (vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+1))*(vect(icp+3)-vect(icp+2)))                           * &
                   (2.0/(vect(icp+3)-vect(icp+1)) + 2.0/(vect(icp+4)-vect(icp+2)))
						
                else if(eval_points(ipos)>=vect(icp+3).AND.eval_points(ipos)<vect(icp+4)) then 
                    Splines(icoil)%db_dt_2(ipos,icp) =  &
                   -2.0/(vect(icp+4)-vect(icp+1))      *(-Splines(icoil)%basis_1(ipos,icp+2)/(vect(icp+4)-vect(icp+2))  - &
                        (vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+2))*(vect(icp+4)-vect(icp+3))))          + &
                    2.0*(vect(icp+4)-eval_points(ipos))/((vect(icp+4)-vect(icp+1))*(vect(icp+4)-vect(icp+2))*(vect(icp+4)-vect(icp+3)))
						
                endif
        enddo
    enddo


    return 
END SUBROUTINE eval_spline_basis2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE enforce_spline_periodicity(icoil)
!Copy derivatives in respect to last three control points in place of the derivates in respect to first three (they are the same points)
    use globals, only : dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet,coil,Splines,astat  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0: coil(icoil)%NS-1)
    integer :: icp,ipos


    N = coil(icoil)%NS
    eval_points = Splines(icoil)%eval_points

    do ipos = 0,N-1
        Splines(icoil)%basis_3(ipos,0:2)  = Splines(icoil)%basis_3(ipos,0:2) + &
			                    Splines(icoil)%basis_3(ipos,Splines(icoil)%NCP-3:Splines(icoil)%NCP-1)
	Splines(icoil)%db_dt  (ipos,0:2)  = Splines(icoil)%db_dt  (ipos,0:2) + &
				            Splines(icoil)%db_dt  (ipos,Splines(icoil)%NCP-3:Splines(icoil)%NCP-1)
	Splines(icoil)%db_dt_2(ipos,0:2)  = Splines(icoil)%db_dt_2(ipos,0:2) + & 
				            Splines(icoil)%db_dt_2(ipos,Splines(icoil)%NCP-3:Splines(icoil)%NCP-1)

	Splines(icoil)%basis_3(ipos,Splines(icoil)%NCP-3:Splines(icoil)%NCP-1)  =0
	Splines(icoil)%db_dt  (ipos,Splines(icoil)%NCP-3:Splines(icoil)%NCP-1)  =0
	Splines(icoil)%db_dt_2(ipos,Splines(icoil)%NCP-3:Splines(icoil)%NCP-1)  =0	
    enddo		
END SUBROUTINE enforce_spline_periodicity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_eval_spline_basis(icoil)
  !Check numerically the first and second derivatives of the basis functions and print difference on console
    use globals, only : dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet,coil,Splines,astat   
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: icp,ipos
    REAL :: vect(0:Splines(icoil)%NT-1)
    REAL :: sigma = 1E-06
    REAL :: bas_temp_3(0:coil(icoil)%NS-1,0:Splines(icoil)%NCP-1)
    REAL :: db_dt_temp(0:coil(icoil)%NS-1,0:Splines(icoil)%NCP-1)	
       
    vect = Splines(icoil)%vect
    N = coil(icoil)%NS 
    eval_points = Splines(icoil)%eval_points
    Splines(icoil)%eval_points= Splines(icoil)%eval_points + sigma

    call eval_spline_basis(icoil)	
    bas_temp_3 = Splines(icoil)%basis_3

    Splines(icoil)%eval_points = Splines(icoil)%eval_points - 2.0*sigma
    call eval_spline_basis(icoil)

    do ipos = 0,coil(icoil)%NS-1
        do icp = 0,Splines(icoil)%NCP-1
	    if (myid == 0) then
	        if (     Splines(icoil)%db_dt(ipos,icp) > 0 .AND. &
		    ABS((Splines(icoil)%db_dt(ipos,icp) - (bas_temp_3(ipos,icp) - Splines(icoil)%basis_3(ipos,icp))/2.0/sigma ))> 0 )     &
		        write (ounit,'("first pos " I7.3 " CP " I7.3 " basis 3 der" E15.7 ",rel " E15.7 " " E15.7 " " E15.7)') ipos,icp,  &
		        ( Splines(icoil)%db_dt(ipos,icp) - (bas_temp_3(ipos,icp) - Splines(icoil)%basis_3(ipos,icp))/2.0/sigma ),         &
                        ((Splines(icoil)%db_dt(ipos,icp) - (bas_temp_3(ipos,icp) - Splines(icoil)%basis_3(ipos,icp))/2.0/sigma ))/        &
			  Splines(icoil)%db_dt(ipos,icp),                                                                                 &
		          Splines(icoil)%db_dt(ipos,icp),                                                                                 &
		        (bas_temp_3(ipos,icp) - Splines(icoil)%basis_3(ipos,icp))/2.0/sigma
	    endif 
	enddo
    enddo


    Splines(icoil)%eval_points = eval_points 
    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    db_dt_temp = Splines(icoil)%db_dt

    Splines(icoil)%eval_points = Splines(icoil)%eval_points - eval_points*2.0*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)


    do ipos = 0,coil(icoil)%NS-1
	do icp = 0,Splines(icoil)%NCP-1
            if (myid == 0) then
	        if (eval_points(ipos) >0 .AND. Splines(icoil)%db_dt_2(ipos,icp) > 0 .AND. &
		    ABS((Splines(icoil)%db_dt_2(ipos,icp) - (db_dt_temp(ipos,icp) - Splines(icoil)%db_dt(ipos,icp))/2.0/eval_points(ipos)/sigma ))> 0 ) &
		    write (ounit,'("second pos " I7.3 " CP " I7.3 " basis 3 der" E15.7 ",rel " E15.7 " " E15.7 " " E15.7)') ipos,icp,            &
                    ( Splines(icoil)%db_dt_2(ipos,icp) - (db_dt_temp(ipos,icp) - Splines(icoil)%db_dt(ipos,icp))/2.0/ eval_points(ipos)/sigma ), &
		    ((Splines(icoil)%db_dt_2(ipos,icp) - (db_dt_temp(ipos,icp) - Splines(icoil)%db_dt(ipos,icp))/2.0/eval_points(ipos)/sigma ))/ &
		      Splines(icoil)%db_dt_2(ipos,icp),                                                                                          &
	              Splines(icoil)%db_dt_2(ipos,icp), &
		    (db_dt_temp(ipos,icp) - Splines(icoil)%db_dt(ipos,icp))/2.0/eval_points(ipos)/sigma
	    endif
	enddo
    enddo

    Splines(icoil)%eval_points = eval_points 
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
END SUBROUTINE check_eval_spline_basis

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_xt_xa(icoil)
   !check numerically tangent vectors  
    use globals, only : dp, zero, one, myid, ounit, sqrtmachprec, IsQuiet,coil,Splines,astat   
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: icp,ipos
    REAL :: vect(0:Splines(icoil)%NT-1)
    REAL :: sigma = 1E-6
    REAL :: x_temp(0:coil(icoil)%NS)
    REAL :: xt_temp(0:coil(icoil)%NS)
    REAL :: y_temp(0:coil(icoil)%NS)
    REAL :: yt_temp(0:coil(icoil)%NS)
    REAL :: z_temp(0:coil(icoil)%NS)
    REAL :: zt_temp(0:coil(icoil)%NS)	
       
    vect = Splines(icoil)%vect
    N = coil(icoil)%NS
   
    eval_points = Splines(icoil)%eval_points

    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call discoil(1) 
    x_temp = coil(icoil)%xx
 
    Splines(icoil)%eval_points= Splines(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_spline_basis(icoil)
    call discoil(1) 

    do ipos = 0,coil(icoil)%NS-1
	if (.true.) then
	    if (eval_points(ipos) >0 .AND. coil(icoil)%xt(ipos) > 0 ) &
	        write (ounit,'("xt pos " I7.3 " xt der" E15.7 ",rel " E15.7 " " E15.7 " " E15.7)') ipos,                            &
		( coil(icoil)%xt(ipos) - (x_temp(ipos) - coil(icoil)%xx(ipos))/2.0/ eval_points(ipos)/sigma ),                      &
		((coil(icoil)%xt(ipos) - (x_temp(ipos) - coil(icoil)%xx(ipos))/2.0/eval_points(ipos)/sigma ))/coil(icoil)%xt(ipos), &
	          coil(icoil)%xt(ipos),                                                                                             &
		(x_temp(ipos) - coil(icoil)%xx(ipos))/2.0/eval_points(ipos)/sigma
	endif
    enddo

    Splines(icoil)%eval_points = eval_points
    call eval_spline_basis(icoil)
    call discoil(1) 

    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 
    xt_temp = coil(icoil)%xt

    Splines(icoil)%eval_points= Splines(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

    do ipos = 0,coil(icoil)%NS-1
	if (myid==0) then
	    if (eval_points(ipos) >0 .AND. eval_points(ipos) < 2.AND. coil(icoil)%xa(ipos) > 0 ) &
	        write (ounit,'("xa pos " I7.3 " xa der" E15.7 ",rel " E15.7 "an " E15.7 "num " E15.7)') ipos,   &
	        ( coil(icoil)%xa(ipos) - (xt_temp(ipos) - coil(icoil)%xt(ipos))/2.0/eval_points(ipos)/sigma ),  &
		((coil(icoil)%xa(ipos) - (xt_temp(ipos) - coil(icoil)%xt(ipos))/2.0/eval_points(ipos)/sigma ))/ &
	          coil(icoil)%xa(ipos),                                                                         &
	 	  coil(icoil)%xa(ipos),                                                                         &
		(xt_temp(ipos) - coil(icoil)%xt(ipos))/2.0/eval_points(ipos)/sigma
	endif

    enddo

    Splines(icoil)%eval_points = eval_points
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call discoil(1) 
    y_temp = coil(icoil)%yy
 
    Splines(icoil)%eval_points= Splines(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_spline_basis(icoil)
    call discoil(1) 

    do ipos = 0,coil(icoil)%NS-1
	if (.true.) then
	    if (eval_points(ipos) >0 .AND. eval_points(ipos) < 2.AND. coil(icoil)%yt(ipos) > 0 ) &
		write (ounit,'("yt pos " I7.3 " yt der" E15.7 ",rel " E15.7 " " E15.7 " " E15.7)') ipos,         &
	        ( coil(icoil)%yt(ipos) - (y_temp(ipos) - coil(icoil)%yy(ipos))/2.0/eval_points(ipos)/sigma ),    &
		((coil(icoil)%yt(ipos) - (y_temp(ipos) - coil(icoil)%yy(ipos))/2.0/eval_points(ipos)/sigma ))/   &
		  coil(icoil)%yt(ipos),                                                                          &
	          coil(icoil)%yt(ipos),                                                                          &
		(y_temp(ipos) - coil(icoil)%yy(ipos))/2.0/eval_points(ipos)/sigma
	endif

    enddo

    Splines(icoil)%eval_points = eval_points
    call eval_spline_basis(icoil)
    call discoil(1) 

    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 
    yt_temp = coil(icoil)%yt
 
    Splines(icoil)%eval_points= Splines(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

    do ipos = 0,coil(icoil)%NS-1
	if (myid==0) then
	    if (eval_points(ipos) >0 .AND. eval_points(ipos) < 2.AND. coil(icoil)%ya(ipos) > 0 ) &
	        write (ounit,'("ya pos " I7.3 " ya der" E15.7 ",rel " E15.7 "an " E15.7 "num " E15.7)') ipos,    &
		( coil(icoil)%ya(ipos) - (yt_temp(ipos) - coil(icoil)%yt(ipos))/2.0/eval_points(ipos)/sigma ),   &
		((coil(icoil)%ya(ipos) - (yt_temp(ipos) - coil(icoil)%yt(ipos))/2.0/eval_points(ipos)/sigma ))/  &
		  coil(icoil)%ya(ipos),                                                                          &
		  coil(icoil)%ya(ipos),                                                                          &
		(yt_temp(ipos) - coil(icoil)%yt(ipos))/2.0/eval_points(ipos)/sigma
	endif

    enddo

    Splines(icoil)%eval_points = eval_points
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 
	
    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

    zt_temp = coil(icoil)%zt
    Splines(icoil)%eval_points= Splines(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

    do ipos = 0,coil(icoil)%NS-1
	if (myid==0) then
	    if (eval_points(ipos) >0 .AND. eval_points(ipos) < 2.AND. coil(icoil)%za(ipos) > 0 ) &
	        write (ounit,'("za pos " I7.3 " za der" E15.7 ",rel " E15.7 "an " E15.7 "num " E15.7)') ipos,     &
	 	( coil(icoil)%za(ipos) - (zt_temp(ipos) - coil(icoil)%zt(ipos))/2.0/eval_points(ipos)/sigma ),    &
		((coil(icoil)%za(ipos) - (zt_temp(ipos) - coil(icoil)%zt(ipos))/2.0/eval_points(ipos)/sigma ))/   &
		  coil(icoil)%za(ipos),                                                                           &
	          coil(icoil)%za(ipos),                                                                           &
		(zt_temp(ipos) - coil(icoil)%zt(ipos))/2.0/eval_points(ipos)/sigma
	endif

    enddo

    Splines(icoil)%eval_points = eval_points
    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

    Splines(icoil)%eval_points= Splines(icoil)%eval_points + eval_points*sigma
    call eval_spline_basis(icoil)
    call discoil(1) 
    z_temp = coil(icoil)%zz
 
    Splines(icoil)%eval_points= Splines(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_spline_basis(icoil)
    call discoil(1) 
    do ipos = 0,coil(icoil)%NS-1
	if (.true.) then
	    if (eval_points(ipos) >0 .AND. eval_points(ipos) < 2 .AND. coil(icoil)%zt(ipos) > 0 )                &
	        write (ounit,'("zt pos " I7.3 " zt der" E15.7 ",rel " E15.7 " " E15.7 " " E15.7)') ipos,         &
	 	( coil(icoil)%zt(ipos) - (z_temp(ipos) - coil(icoil)%zz(ipos))/2.0/eval_points(ipos)/sigma ),    &
		((coil(icoil)%zt(ipos) - (z_temp(ipos) - coil(icoil)%zz(ipos))/2.0/eval_points(ipos)/sigma ))/   &
		  coil(icoil)%zt(ipos),                                                                          &
		  coil(icoil)%zt(ipos), &
		(z_temp(ipos) - coil(icoil)%zz(ipos))/2.0/eval_points(ipos)/sigma
	endif

    enddo

    Splines(icoil)%eval_points = eval_points

    call eval_spline_basis(icoil)
    call eval_spline_basis1(icoil)
    call discoil(1) 

END SUBROUTINE check_xt_xa	

SUBROUTINE eval_basis(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0: coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos

    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points


	CPCoil(icoil)%basis_0 = 0
	CPCoil(icoil)%basis_1 = 0
	CPCoil(icoil)%basis_2 = 0
	CPCoil(icoil)%basis_3 = 0


    do iter_pos=0,N-1 
    	do iter_cp=0,CPCoil(icoil)%NCP+2
                if(eval_points(iter_pos)>=CPCoil(icoil)%vect(iter_cp).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(iter_cp+1))then
                 CPCoil(icoil)%basis_0(iter_pos,iter_cp) = 1
		!if (iter_cp <CPCoil(icoil)%NCP .AND. myid==0 .AND. iter_pos > 7 .AND. iter_pos < 9) &
		!							write (ounit,'(" 0 pos " I7.3"cp  " I7.3 "x " F15.7 "  "F15.7)')iter_pos,iter_cp,&
		!									 eval_points(iter_pos),CPCoil(icoil)%basis_0(iter_pos,iter_cp) 
		endif
        enddo
    enddo
	!Loop around basis functions to enforce periodicity
  	  !do iter_pos=0,N-1 
!	if(eval_points(iter_pos)>=CPCoil(icoil)%vect(CPCoil(icoil)%NCP-3).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(CPCoil(icoil)%NCP-2)) then
!                CPCoil(icoil)%basis_0(iter_pos,3) = 1.0/2
!	 	 CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-3) = 1.0/2
!       elseif(eval_points(iter_pos)>=CPCoil(icoil)%vect(CPCoil(icoil)%NCP-2).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(CPCoil(icoil)%NCP-1))then
!               CPCoil(icoil)%basis_0(iter_pos,4) = 1.0/2
!	 	 CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-2) = 1.0/2
!       elseif(eval_points(iter_pos)>=CPCoil(icoil)%vect(CPCoil(icoil)%NCP-1).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(CPCoil(icoil)%NCP))then
!               CPCoil(icoil)%basis_0(iter_pos,5) = 1.0/2
!	 	 CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-1) = 1.0/2
!	elseif(eval_points(iter_pos)>=CPCoil(icoil)%vect(3).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(4)) then
!                CPCoil(icoil)%basis_0(iter_pos,3) = 1.0/2
!	 	 CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-3) = 1.0/2
!       elseif(eval_points(iter_pos)>=CPCoil(icoil)%vect(4).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(5))then
!               CPCoil(icoil)%basis_0(iter_pos,4) = 1.0/2
!	 	 CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-2) = 1.0/2
!       elseif(eval_points(iter_pos)>=CPCoil(icoil)%vect(5).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(6))then
!               CPCoil(icoil)%basis_0(iter_pos,5) = 1.0/2
!	 	 CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-1) = 1.0/2
!	endif
!   enddo

    do iter_cp=0,CPCoil(icoil)%NCP+1   
        do iter_pos=0,N -1
                if( (CPCoil(icoil)%basis_0(iter_pos,iter_cp) /= 0 &
									.OR. CPCoil(icoil)%basis_0(iter_pos,iter_cp+1) /= 0))then
                         CPCoil(icoil)%basis_1(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp)) / &
                        (CPCoil(icoil)%vect(iter_cp+1) - CPCoil(icoil)%vect(iter_cp))*CPCoil(icoil)%basis_0(iter_pos,iter_cp) &
                         + 1.0*( - eval_points(iter_pos) + CPCoil(icoil)%vect(iter_cp+2))/ &
                         (CPCoil(icoil)%vect(iter_cp+2) - CPCoil(icoil)%vect(iter_cp+1))*CPCoil(icoil)%basis_0(iter_pos,iter_cp+1)
		!if((eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp+1)) < 0.001) CPCoil(icoil)%basis_1(iter_pos,iter_cp) = 1.0  !algorithm gives problems on knots
		!if (iter_cp <CPCoil(icoil)%NCP .AND. myid==0 .AND. iter_pos > 7 .AND. iter_pos < 9) &
		!							write (ounit,'(" 1 pos " I7.3"cp  " I7.3 "x " F15.7 "  "F15.7)')iter_pos,iter_cp,&
		!									 eval_points(iter_pos),CPCoil(icoil)%basis_1(iter_pos,iter_cp) 
		endif
                enddo    
        enddo

    do iter_cp=0,CPCoil(icoil)%NCP
        do iter_pos=0,N -1
                if (  ((CPCoil(icoil)%basis_1(iter_pos,iter_cp) .NE. 0) &
									.OR. (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1) .NE. 0))) &
                         CPCoil(icoil)%basis_2(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp))/ & 
                        (CPCoil(icoil)%vect(iter_cp+2) - CPCoil(icoil)%vect(iter_cp))*CPCoil(icoil)%basis_1(iter_pos,iter_cp) &
                        + 1.0*( - eval_points(iter_pos) + CPCoil(icoil)%vect(iter_cp+3))/ &
                        (CPCoil(icoil)%vect(iter_cp+3) - CPCoil(icoil)%vect(iter_cp+1))*CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)
		!if (iter_cp <CPCoil(icoil)%NCP .AND. myid==0 .AND. iter_pos > 7 .AND. iter_pos < 9) &
	!					write (ounit,'(" 2 pos " I7.3"cp  " I7.3 "x " F15.7 "  "F15.7)')iter_pos,iter_cp,&
	!										 eval_points(iter_pos),CPCoil(icoil)%basis_2(iter_pos,iter_cp) 
                enddo    
        enddo
    
    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if (((CPCoil(icoil)%basis_2(iter_pos,iter_cp) .NE. 0) &
							.OR. (CPCoil(icoil)%basis_2(iter_pos,iter_cp+1) .NE. 0))) then
                         CPCoil(icoil)%basis_3(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp))/ &
                        (CPCoil(icoil)%vect(iter_cp+3) - CPCoil(icoil)%vect(iter_cp))*CPCoil(icoil)%basis_2(iter_pos,iter_cp) &
                         + 1.0*( - eval_points(iter_pos) + CPCoil(icoil)%vect(iter_cp+4))/ &
                        (CPCoil(icoil)%vect(iter_cp+4) - CPCoil(icoil)%vect(iter_cp+1))*CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)
		!if (myid==0 .AND. iter_pos > 7 .AND. iter_pos < 9)write (ounit,'(" 3 pos " I7.3"cp  " I7.3 "x " F15.7 "  "F15.7)')iter_pos,iter_cp,&
		!									 eval_points(iter_pos),CPCoil(icoil)%basis_3(iter_pos,iter_cp) 
		endif
                enddo
        enddo




    return 
END SUBROUTINE eval_basis

SUBROUTINE eval_basis1(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT-1)
       
    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points
	CPCoil(icoil)%db_dt=0

    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=vect(iter_cp).AND.eval_points(iter_pos)<vect(iter_cp+1)) then
                                                 CPCoil(icoil)%db_dt(iter_pos,iter_cp) = &
                                                 CPCoil(icoil)%basis_2(iter_pos,iter_cp)/(vect(iter_cp+3)-vect(iter_cp)) &
                                                + (eval_points(iter_pos)-vect(iter_cp))/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) + &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+1)-vect(iter_cp))))
						!write(ounit,'("4 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                elseif(eval_points(iter_pos)>=vect(iter_cp+1).AND.eval_points(iter_pos)<vect(iter_cp+2)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) =  &
                                                CPCoil(icoil)%basis_2(iter_pos,iter_cp)/(vect(iter_cp+3)-vect(iter_cp)) &
                                                + (eval_points(iter_pos)-vect(iter_cp))/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) - &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1))) & 
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                 (vect(iter_cp+2)-vect(iter_cp+1)))) &
                                                - CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) + &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1))))
						!write(ounit,'("3 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                elseif(eval_points(iter_pos)>=vect(iter_cp+2).AND.eval_points(iter_pos)<vect(iter_cp+3)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) =  &
                                                 CPCoil(icoil)%basis_2(iter_pos,iter_cp)/(vect(iter_cp+3)-vect(iter_cp)) &
                                                + (eval_points(iter_pos)-vect(iter_cp))/(vect(iter_cp+3)-vect(iter_cp))* &
                                                 (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                - (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* & 
                                                (vect(iter_cp+3)-vect(iter_cp+2)))) &
                                                - CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) - &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2))) &
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+2)/(vect(iter_cp+4)-vect(iter_cp+2)) + &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2))))
						!write(ounit,'("2 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                elseif(eval_points(iter_pos)>=vect(iter_cp+3).AND.eval_points(iter_pos)<vect(iter_cp+4)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) = &
                                                 -CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+2)/(vect(iter_cp+4)-vect(iter_cp+2)) - &
                                                (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+4)-vect(iter_cp+3))))
						!write(ounit,'("1 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                endif
        enddo
    enddo


    return 
END SUBROUTINE eval_basis1

SUBROUTINE eval_basis2(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT-1)

    CPCoil(icoil)%db_dt_2 = 0
    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points
    
    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=vect(iter_cp).AND.eval_points(iter_pos)<vect(iter_cp+1)) then
                                                 CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 2.0/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) + &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+1)-vect(iter_cp)))) &
                                                + 2.0*(eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+3)-vect(iter_cp))* &
                                                (vect(iter_cp+1)-vect(iter_cp))*(vect(iter_cp+2)-vect(iter_cp)))
						!write(ounit,'("4 is " 7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                else if(eval_points(iter_pos)>=vect(iter_cp+1).AND.eval_points(iter_pos)<vect(iter_cp+2)) then 
                                               CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 2.0/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) - &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1))) & 
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1)))) &
                                                - (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+3)-vect(iter_cp))* &
                                                 (vect(iter_cp+2)-vect(iter_cp+1)))* &
                                                (2.0/(vect(iter_cp+2)-vect(iter_cp)) + 2.0/(vect(iter_cp+3)-vect(iter_cp+1))) &
                                                - 2.0/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) + &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1)))) + &
                                                2.0*(vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+1)) * (vect(iter_cp+2)-vect(iter_cp+1))) 
						!write(ounit,'("3 is "7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                else if(eval_points(iter_pos)>=vect(iter_cp+2).AND.eval_points(iter_pos)<vect(iter_cp+3)) then
                                                 CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 2.0/(vect(iter_cp+3)-vect(iter_cp))* &   
                                                 (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                - (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                 (vect(iter_cp+3)-vect(iter_cp+2)))) &
                                                + 2.0*(eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+3)-vect(iter_cp))* &
                                                (vect(iter_cp+3)-vect(iter_cp+1))*(vect(iter_cp+3)-vect(iter_cp+2))) &
                                                - 2.0/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) - &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2))) &
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+2)/(vect(iter_cp+4)-vect(iter_cp+2)) + &
                                                 (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2)))) &
                                                - (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2)))* &
                                                (2.0/(vect(iter_cp+3)-vect(iter_cp+1)) + 2.0/(vect(iter_cp+4)-vect(iter_cp+2)))
						!write(ounit,'("2 is " 7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                else if(eval_points(iter_pos)>=vect(iter_cp+3).AND.eval_points(iter_pos)<vect(iter_cp+4)) then 
                                                CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) =  &
                                                 -2.0/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+2)/(vect(iter_cp+4)-vect(iter_cp+2)) - &
                                                (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+4)-vect(iter_cp+3)))) &
                                                + 2.0*(vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (vect(iter_cp+4)-vect(iter_cp+2))*(vect(iter_cp+4)-vect(iter_cp+3)))
						!write(ounit,'("1 is " 7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                endif
        enddo
    enddo


    return 
END SUBROUTINE eval_basis2

SUBROUTINE enforce_periodicity(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0: coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos


    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points

    do iter_pos = 0,N-1
	
		CPCoil(icoil)%basis_3(iter_pos,0:2)=CPCoil(icoil)%basis_3(iter_pos,0:2) + &
						CPCoil(icoil)%basis_3(iter_pos,CPCoil(icoil)%NCP-3:CPCoil(icoil)%NCP-1)
		!CPCoil(icoil)%basis_2(iter_pos,0:3)=CPCoil(icoil)%basis_2(iter_pos,CPCoil(icoil)%NCP-4:CPCoil(icoil)%NCP-1)
		!CPCoil(icoil)%basis_1(iter_pos,0:4)=CPCoil(icoil)%basis_1(iter_pos,CPCoil(icoil)%NCP-5:CPCoil(icoil)%NCP-1)
		!CPCoil(icoil)%basis_0(iter_pos,0:5)=CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-6:CPCoil(icoil)%NCP-1)
		CPCoil(icoil)%db_dt(iter_pos,0:2)  = CPCoil(icoil)%db_dt(iter_pos,0:2) + &
							CPCoil(icoil)%db_dt(iter_pos,CPCoil(icoil)%NCP-3:CPCoil(icoil)%NCP-1)
		CPCoil(icoil)%db_dt_2(iter_pos,0:2)=CPCoil(icoil)%db_dt_2(iter_pos,0:2) + & 
							CPCoil(icoil)%db_dt_2(iter_pos,CPCoil(icoil)%NCP-3:CPCoil(icoil)%NCP-1)
		CPCoil(icoil)%basis_3(iter_pos,CPCoil(icoil)%NCP-3:CPCoil(icoil)%NCP-1)=0
		!CPCoil(icoil)%basis_2(iter_pos,CPCoil(icoil)%NCP-4:CPCoil(icoil)%NCP-1)=0
		!CPCoil(icoil)%basis_1(iter_pos,CPCoil(icoil)%NCP-5:CPCoil(icoil)%NCP-1)=0
		!CPCoil(icoil)%basis_0(iter_pos,CPCoil(icoil)%NCP-6:CPCoil(icoil)%NCP-1)=0
		CPCoil(icoil)%db_dt(iter_pos,CPCoil(icoil)%NCP-3:CPCoil(icoil)%NCP-1)  =0
		CPCoil(icoil)%db_dt_2(iter_pos,CPCoil(icoil)%NCP-3:CPCoil(icoil)%NCP-1)=0

	
    enddo		
END SUBROUTINE enforce_periodicity

SUBROUTINE eval_deriv1(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT-1)
       
	CPCoil(icoil)%db_dt=0
    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points

    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=vect(iter_cp).AND.eval_points(iter_pos)<vect(iter_cp+1)) then
                                                 CPCoil(icoil)%db_dt(iter_pos,iter_cp) = &
                                                 3.0 * (eval_points(iter_pos)-vect(iter_cp))**2/((vect(iter_cp+1)-vect(iter_cp)) * &
						(vect(iter_cp+2)-vect(iter_cp))*(vect(iter_cp+3)-vect(iter_cp)))
						!write(ounit,'("4 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                elseif(eval_points(iter_pos)>=vect(iter_cp+1).AND.eval_points(iter_pos)<vect(iter_cp+2)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) =  &
                                                -(eval_points(iter_pos) - vect(iter_cp+1))**2/((vect(iter_cp+2) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+4) - vect(iter_cp+1))) + &
						(2.0*(vect(iter_cp+4) - eval_points(iter_pos))*(eval_points(iter_pos) - vect(iter_cp+1)))/ &
						((vect(iter_cp+2) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+1))* &
						(vect(iter_cp+4) - vect(iter_cp+1))) + ((eval_points(iter_pos) - vect(iter_cp))* &
						(-(eval_points(iter_pos) - vect(iter_cp))/((vect(iter_cp+2) - vect(iter_cp))* &
						(vect(iter_cp+2) - vect(iter_cp+1))) + (vect(iter_cp+2) - eval_points(iter_pos))/ & 
						((vect(iter_cp+2) - vect(iter_cp))*(vect(iter_cp+2) - vect(iter_cp+1))) + &
						(vect(iter_cp+3) - eval_points(iter_pos))/((vect(iter_cp+2) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+1))) - (eval_points(iter_pos) - vect(iter_cp+1))/ &
						((vect(iter_cp+2) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+1)))))/ &
						(vect(iter_cp+3) - vect(iter_cp)) + (((eval_points(iter_pos) - vect(iter_cp))* &
						(vect(iter_cp+2) - eval_points(iter_pos)))/((vect(iter_cp+2) - vect(iter_cp))* &
						(vect(iter_cp+2) - vect(iter_cp+1))) + ((eval_points(iter_pos) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - eval_points(iter_pos)))/((vect(iter_cp+2) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+1))))/(vect(iter_cp+3) - vect(iter_cp))
						!write(ounit,'("3 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                elseif(eval_points(iter_pos)>=vect(iter_cp+2).AND.eval_points(iter_pos)<vect(iter_cp+3)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) =  &
                                                (vect(iter_cp+3) - eval_points(iter_pos))**2/((vect(iter_cp+3) - vect(iter_cp))* &
						(vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+2))) - &
						(2.0*(eval_points(iter_pos) - vect(iter_cp))*(vect(iter_cp+3) - eval_points(iter_pos)))/ &
						((vect(iter_cp+3) - vect(iter_cp))*(vect(iter_cp+3) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+2))) + ((vect(iter_cp+4) - eval_points(iter_pos))* &
						(-(eval_points(iter_pos) - vect(iter_cp+1))/((vect(iter_cp+3) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+2))) + (vect(iter_cp+3) - eval_points(iter_pos))/ &
						((vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+2))) + &
						(vect(iter_cp+4) - eval_points(iter_pos))/((vect(iter_cp+3) - vect(iter_cp+2))* &
						(vect(iter_cp+4) - vect(iter_cp+2))) - (eval_points(iter_pos) - vect(iter_cp+2))/ &
						((vect(iter_cp+3) - vect(iter_cp+2))*(vect(iter_cp+4) - vect(iter_cp+2)))))/ &
						(vect(iter_cp+4) - vect(iter_cp+1)) - (((eval_points(iter_pos) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - eval_points(iter_pos)))/((vect(iter_cp+3) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+2))) + ((eval_points(iter_pos) - vect(iter_cp+2))* &
						(vect(iter_cp+4) - eval_points(iter_pos)))/((vect(iter_cp+3) - vect(iter_cp+2))* &
						(vect(iter_cp+4) - vect(iter_cp+2))))/(vect(iter_cp+4) - vect(iter_cp+1))
						!write(ounit,'("2 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                elseif(eval_points(iter_pos)>=vect(iter_cp+3).AND.eval_points(iter_pos)<vect(iter_cp+4)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) = &
                                                (3.0*(vect(iter_cp+4) - eval_points(iter_pos))**2)/((vect(iter_cp+1) - vect(iter_cp+4))* &
						 (vect(iter_cp+4) - vect(iter_cp+2))*(vect(iter_cp+4) - vect(iter_cp+3)))
						!write(ounit,'("1 is " 7F20.10)')CPCoil(icoil)%db_dt(iter_pos,iter_cp)
                endif
        enddo
    enddo


    return 
END SUBROUTINE eval_deriv1

SUBROUTINE eval_deriv2(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT-1)

	CPCoil(icoil)%db_dt_2=0
    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points
    
    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=vect(iter_cp).AND.eval_points(iter_pos)<vect(iter_cp+1)) then
						CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 6.0 * (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+1)-vect(iter_cp)) * &
						(vect(iter_cp+2)-vect(iter_cp))*(vect(iter_cp+3)-vect(iter_cp)))			
						!write(ounit,'("4 is " 7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                else if(eval_points(iter_pos)>=vect(iter_cp+1).AND.eval_points(iter_pos)<vect(iter_cp+2)) then 
                                               CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                -(4.0*(eval_points(iter_pos) - vect(iter_cp+1)))/((vect(iter_cp+2) - vect(iter_cp+1))* &
						 (vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+4) - vect(iter_cp+1))) + & 
						((-2.0/((vect(iter_cp+2) - vect(iter_cp))*(vect(iter_cp+2) - vect(iter_cp+1))) - 2.0/ &
						((vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+2) - vect(iter_cp+1))))* & 
						(eval_points(iter_pos) - vect(iter_cp)))/(vect(iter_cp+3) - vect(iter_cp)) + (2.0* &
						(-(eval_points(iter_pos) - vect(iter_cp))/((vect(iter_cp+2) - vect(iter_cp))* &
						(vect(iter_cp+2) - vect(iter_cp+1))) + (vect(iter_cp+2) - eval_points(iter_pos))/ & 
						((vect(iter_cp+2) - vect(iter_cp))*(vect(iter_cp+2) - vect(iter_cp+1))) + &
						 (vect(iter_cp+3) - eval_points(iter_pos))/((vect(iter_cp+2) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+1))) - (eval_points(iter_pos) - vect(iter_cp+1))/ &
						((vect(iter_cp+2) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+1)))))/ &
						(vect(iter_cp+3) - vect(iter_cp)) + (2.0*(vect(iter_cp+4) - eval_points(iter_pos)))/ &
						((vect(iter_cp+2) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+1))* &
						(vect(iter_cp+4) - vect(iter_cp+1)))
						!write(ounit,'("3 is "7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                else if(eval_points(iter_pos)>=vect(iter_cp+2).AND.eval_points(iter_pos)<vect(iter_cp+3)) then
                                                 CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 (2.0*(eval_points(iter_pos) - vect(iter_cp)) - 4.0*(vect(iter_cp+3) - eval_points(iter_pos)))/ &
						((vect(iter_cp+3) - vect(iter_cp))*(vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+2)))&
						 + ((-2.0/((vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+2))) - 2.0/ &
						((vect(iter_cp+4) - vect(iter_cp+2))*(vect(iter_cp+3) - vect(iter_cp+2))))* &
						(vect(iter_cp+4) - eval_points(iter_pos)))/(vect(iter_cp+4) - vect(iter_cp+1)) - &
						(2.0*(-(eval_points(iter_pos) - vect(iter_cp+1))/((vect(iter_cp+3) - vect(iter_cp+1))* &
						(vect(iter_cp+3) - vect(iter_cp+2))) + (vect(iter_cp+3) - eval_points(iter_pos))/ &
						((vect(iter_cp+3) - vect(iter_cp+1))*(vect(iter_cp+3) - vect(iter_cp+2))) + &
						(vect(iter_cp+4) - eval_points(iter_pos))/((vect(iter_cp+3) - vect(iter_cp+2))* &
						(vect(iter_cp+4) - vect(iter_cp+2))) - (eval_points(iter_pos) - vect(iter_cp+2))/ &
						((vect(iter_cp+3) - vect(iter_cp+2))*(vect(iter_cp+4) - vect(iter_cp+2)))))/ &
						(vect(iter_cp+4) - vect(iter_cp+1))
						!write(ounit,'("2 is " 7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                else if(eval_points(iter_pos)>=vect(iter_cp+3).AND.eval_points(iter_pos)<vect(iter_cp+4)) then 
                                                CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) =  &
                                                (6.0*(vect(iter_cp+4) - eval_points(iter_pos)))/ &
						((vect(iter_cp+4) - vect(iter_cp+1))*(vect(iter_cp+4) - vect(iter_cp+2))* &
						(vect(iter_cp+4) - vect(iter_cp+3)))
						!write(ounit,'("1 is " 7F20.10)')CPCoil(icoil)%db_dt_2(iter_pos,iter_cp)
                endif
        enddo
    enddo


    return 
END SUBROUTINE eval_deriv2


SUBROUTINE check_eval_basis(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT-1)
    REAL :: sigma = 0.0001
    !REAL :: bas_temp_0(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    !REAL :: bas_temp_1(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    !REAL :: bas_temp_2(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    REAL :: bas_temp_3(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    REAL :: db_dt_temp(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)	
       
    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
   
    eval_points = CPCoil(icoil)%eval_points
    CPCoil(icoil)%eval_points= CPCoil(icoil)%eval_points + eval_points*sigma

    call eval_basis(icoil)	
    !bas_temp_0 = CPCoil(icoil)%basis_0
    !bas_temp_1 = CPCoil(icoil)%basis_1
    !bas_temp_2 = CPCoil(icoil)%basis_2
    bas_temp_3 = CPCoil(icoil)%basis_3

    CPCoil(icoil)%eval_points = CPCoil(icoil)%eval_points - eval_points*2.0*sigma
    call eval_basis(icoil)


		do iter_pos = 0,coil(icoil)%NS-1
	do iter_cp = 0,CPCoil(icoil)%NCP-1
		if (myid == 0) then
			if (eval_points(iter_pos) >0 .AND. CPCoil(icoil)%db_dt(iter_pos,iter_cp) > 0 .AND. ABS((CPCoil(icoil)%db_dt(iter_pos,iter_cp) - & 
								(bas_temp_3(iter_pos,iter_cp) - CPCoil(icoil)%basis_3(iter_pos,iter_cp))/2.0/ &
								eval_points(iter_pos)/sigma ))> 0 ) &
						 write (ounit,'("first pos " I7.3 " CP " I7.3 " basis 3 der" F15.7 ",rel " F15.7 " " F15.7 " " F15.7)') &
								iter_pos,iter_cp, &
								(CPCoil(icoil)%db_dt(iter_pos,iter_cp) - (bas_temp_3(iter_pos,iter_cp) - &
								CPCoil(icoil)%basis_3(iter_pos,iter_cp))/2.0/ &
								eval_points(iter_pos)/sigma ), ((CPCoil(icoil)%db_dt &
								(iter_pos,iter_cp) - (bas_temp_3(iter_pos,iter_cp) - &
								CPCoil(icoil)%basis_3(iter_pos,iter_cp))/2.0/eval_points(iter_pos)/sigma ))/ &
								CPCoil(icoil)%db_dt(iter_pos,iter_cp), &
								CPCoil(icoil)%db_dt(iter_pos,iter_cp), &
								(bas_temp_3(iter_pos,iter_cp) - &
								CPCoil(icoil)%basis_3(iter_pos,iter_cp))/2.0/ &
								eval_points(iter_pos)/sigma
								!bas_temp_3(iter_pos,iter_cp),CPCoil(icoil)%basis_3(iter_pos,iter_cp)
			endif 
		enddo
	enddo


	CPCoil(icoil)%eval_points = eval_points 
    CPCoil(icoil)%eval_points= CPCoil(icoil)%eval_points + eval_points*sigma
	call eval_basis(icoil)
	call eval_basis1(icoil)
    db_dt_temp = CPCoil(icoil)%db_dt
   CPCoil(icoil)%eval_points = CPCoil(icoil)%eval_points - eval_points*2.0*sigma
	call eval_basis(icoil)
        call eval_basis1(icoil)


		do iter_pos = 0,coil(icoil)%NS-1
	do iter_cp = 0,CPCoil(icoil)%NCP-1
		if (myid == 0) then
		if (eval_points(iter_pos) >0 .AND. CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) > 0 .AND. ABS((CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) - & 
								(db_dt_temp(iter_pos,iter_cp) - CPCoil(icoil)%db_dt(iter_pos,iter_cp))/2.0/ &
								eval_points(iter_pos)/sigma ))> 0 ) &
					write (ounit,'("second pos " I7.3 " CP " I7.3 " basis 3 der" F15.13 ",rel " F15.13 " " F15.13 " " F15.13)') &
								iter_pos,iter_cp, &
								(CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) - (db_dt_temp(iter_pos,iter_cp) - &
								CPCoil(icoil)%db_dt(iter_pos,iter_cp))/2.0/ &
								eval_points(iter_pos)/sigma ), ((CPCoil(icoil)%db_dt_2 &
								(iter_pos,iter_cp) - (db_dt_temp(iter_pos,iter_cp) - &
								CPCoil(icoil)%db_dt(iter_pos,iter_cp))/2.0/eval_points(iter_pos)/sigma ))/ &
								CPCoil(icoil)%db_dt_2(iter_pos,iter_cp), &
								CPCoil(icoil)%db_dt_2(iter_pos,iter_cp), &
								(db_dt_temp(iter_pos,iter_cp) - &
								CPCoil(icoil)%db_dt(iter_pos,iter_cp))/2.0/ &
								eval_points(iter_pos)/sigma
								!bas_temp_3(iter_pos,iter_cp),CPCoil(icoil)%basis_3(iter_pos,iter_cp)
	endif
		enddo
	enddo

	CPCoil(icoil)%eval_points = eval_points 
		call eval_basis(icoil)
	call eval_basis1(icoil)
END SUBROUTINE check_eval_basis

SUBROUTINE check_xt_xa(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT-1)
    REAL :: sigma = 1E-6
    !REAL :: bas_temp_0(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    !REAL :: bas_temp_1(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    !REAL :: bas_temp_2(0:coil(icoil)%NS-1,0:CPCoil(icoil)%NCP-1)
    REAL :: x_temp(0:coil(icoil)%NS-1)
    REAL :: xt_temp(0:coil(icoil)%NS-1)	
       
    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
   
    eval_points = CPCoil(icoil)%eval_points
    CPCoil(icoil)%eval_points= CPCoil(icoil)%eval_points + eval_points*sigma

    call eval_basis(icoil)

    call discoil(1) 

    x_temp = coil(icoil)%xx
 
    CPCoil(icoil)%eval_points= CPCoil(icoil)%eval_points - 2.0*eval_points*sigma

    call eval_basis(icoil)

    call discoil(1) 

	do iter_pos = 0,coil(icoil)%NS-1

	if (.true.) then
		if (eval_points(iter_pos) >0 .AND. coil(icoil)%xx(iter_pos) > 0 ) &
					write (ounit,'("xt pos " I7.3 " xt der" F15.13 ",rel " F15.13 " " F15.13 " " F15.13)') &
								iter_pos, &
								(coil(icoil)%xt(iter_pos) - (x_temp(iter_pos) - &
								coil(icoil)%xx(iter_pos))/2.0/ &
								eval_points(iter_pos)/sigma ), ((coil(icoil)%xt &
								(iter_pos) - (x_temp(iter_pos) - &
								coil(icoil)%xx(iter_pos))/2.0/eval_points(iter_pos)/sigma ))/ &
								coil(icoil)%xt(iter_pos), &
								coil(icoil)%xt(iter_pos), &
								(x_temp(iter_pos) - &
								coil(icoil)%xx(iter_pos))/2.0/ &
								eval_points(iter_pos)/sigma
								!bas_temp_3(iter_pos,iter_cp),CPCoil(icoil)%basis_3(iter_pos,iter_cp)
	endif

	enddo

      CPCoil(icoil)%eval_points = eval_points

    call eval_basis(icoil)
    call discoil(1) 

	CPCoil(icoil)%eval_points= CPCoil(icoil)%eval_points + eval_points*sigma
    call eval_basis(icoil)
    call eval_basis1(icoil)

    call discoil(1) 

    xt_temp = coil(icoil)%xt
 
    CPCoil(icoil)%eval_points= CPCoil(icoil)%eval_points - 2.0*eval_points*sigma
    call eval_basis(icoil)
    call eval_basis1(icoil)

    call discoil(1) 

	do iter_pos = 0,coil(icoil)%NS-1

	if (myid==0) then
		if (eval_points(iter_pos) >0 .AND. coil(icoil)%xt(iter_pos) > 0 ) &
					write (ounit,'("xa pos " I7.3 " xa der" F15.13 ",rel " F15.13 "an " F15.13 "num " F15.13)') &
								iter_pos, &
								(coil(icoil)%xa(iter_pos) - (xt_temp(iter_pos) - &
								coil(icoil)%xt(iter_pos))/2.0/ &
								eval_points(iter_pos)/sigma ), ((coil(icoil)%xa &
								(iter_pos) - (xt_temp(iter_pos) - &
								coil(icoil)%xt(iter_pos))/2.0/eval_points(iter_pos)/sigma ))/ &
								coil(icoil)%xa(iter_pos), &
								coil(icoil)%xa(iter_pos), &
								(xt_temp(iter_pos) - &
								coil(icoil)%xt(iter_pos))/2.0/ &
								eval_points(iter_pos)/sigma
								!bas_temp_3(iter_pos,iter_cp),CPCoil(icoil)%basis_3(iter_pos,iter_cp)
	endif

	enddo

      CPCoil(icoil)%eval_points = eval_points
    call eval_basis(icoil)
    call eval_basis1(icoil)
    call discoil(1) 
	
END SUBROUTINE check_xt_xa	

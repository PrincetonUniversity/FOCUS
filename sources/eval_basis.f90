SUBROUTINE eval_basis(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0: coil(icoil)%NS-1)
    integer :: iter_cp,iter_pos

    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points

    do iter_cp=0,CPCoil(icoil)%NCP+2
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=CPCoil(icoil)%vect(iter_cp).AND.eval_points(iter_pos)<CPCoil(icoil)%vect(iter_cp+1))&
                 CPCoil(icoil)%basis_0(iter_pos,iter_cp) = 1
        enddo
    enddo

    do iter_cp=0,CPCoil(icoil)%NCP+1   
        do iter_pos=0,N -1
                if( CPCoil(icoil)%basis_0(iter_pos,iter_cp) /= 0 .OR. CPCoil(icoil)%basis_0(iter_pos,iter_cp+1) /= 0)&
                         CPCoil(icoil)%basis_1(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp)) / &
                        (CPCoil(icoil)%vect(iter_cp+1) - CPCoil(icoil)%vect(iter_cp))*CPCoil(icoil)%basis_0(iter_pos,iter_cp) &
                         + 1.0*( - eval_points(iter_pos) + CPCoil(icoil)%vect(iter_cp+2))/ &
                         (CPCoil(icoil)%vect(iter_cp+2) - CPCoil(icoil)%vect(iter_cp+1))*CPCoil(icoil)%basis_0(iter_pos,iter_cp+1)
                enddo    
        enddo

    do iter_cp=0,CPCoil(icoil)%NCP
        do iter_pos=0,N -1
                if ( (CPCoil(icoil)%basis_1(iter_pos,iter_cp) .NE. 0) .OR. (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1) .NE. 0)) &
                         CPCoil(icoil)%basis_2(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp))/ & 
                        (CPCoil(icoil)%vect(iter_cp+2) - CPCoil(icoil)%vect(iter_cp))*CPCoil(icoil)%basis_1(iter_pos,iter_cp) &
                        + 1.0*( - eval_points(iter_pos) + CPCoil(icoil)%vect(iter_cp+3))/ &
                        (CPCoil(icoil)%vect(iter_cp+3) - CPCoil(icoil)%vect(iter_cp+1))*CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)
                enddo    
        enddo
    
    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if ( (CPCoil(icoil)%basis_2(iter_pos,iter_cp) .NE. 0) .OR. (CPCoil(icoil)%basis_2(iter_pos,iter_cp+1) .NE. 0))&
                         CPCoil(icoil)%basis_3(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - CPCoil(icoil)%vect(iter_cp))/ &
                        (CPCoil(icoil)%vect(iter_cp+3) - CPCoil(icoil)%vect(iter_cp))*CPCoil(icoil)%basis_2(iter_pos,iter_cp) &
                         + 1.0*( - eval_points(iter_pos) + CPCoil(icoil)%vect(iter_cp+4))/ &
                        (CPCoil(icoil)%vect(iter_cp+4) - CPCoil(icoil)%vect(iter_cp+1))*CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)
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

    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=vect(iter_cp).AND.eval_points(iter_pos)<vect(iter_cp+1)) then
                                                 CPCoil(icoil)%db_dt(iter_pos,iter_cp) = &
                                                 CPCoil(icoil)%basis_2(iter_pos,iter_cp)/(vect(iter_cp+3)-vect(iter_cp)) &
                                                + (eval_points(iter_pos)-vect(iter_cp))/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) + &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+1)-vect(iter_cp))))
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
                elseif(eval_points(iter_pos)>=vect(iter_cp+2).AND.eval_points(iter_pos)<vect(iter_cp+3)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) =  &
                                                 CPCoil(icoil)%basis_2(iter_pos,iter_cp)/(vect(iter_cp+3)-vect(iter_cp)) &
                                                + (eval_points(iter_pos)-vect(iter_cp))/(vect(iter_cp+3)-vect(iter_cp))* &
                                                 (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                - (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* & 
                                                (vect(iter_cp+3)-vect(iter_cp+2)))) &
                                                - CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) + &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2))) &
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+2)) + &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2))))
                elseif(eval_points(iter_pos)>=vect(iter_cp+3).AND.eval_points(iter_pos)<vect(iter_cp+4)) then 
                                                CPCoil(icoil)%db_dt(iter_pos,iter_cp) = &
                                                 -CPCoil(icoil)%basis_2(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+2)) - &
                                                (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+4)-vect(iter_cp+3))))
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

    vect = CPCoil(icoil)%vect
    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points
    
    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N-1 
                if(eval_points(iter_pos)>=vect(iter_cp).AND.eval_points(iter_pos)<vect(iter_cp+1)) then
                                                 CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 2/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) + &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+1)-vect(iter_cp)))) &
                                                + 2*(eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+3)-vect(iter_cp))* &
                                                (vect(iter_cp+1)-vect(iter_cp))*(vect(iter_cp+2)-vect(iter_cp)))
                else if(eval_points(iter_pos)>=vect(iter_cp+1).AND.eval_points(iter_pos)<vect(iter_cp+2)) then 
                                               CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 2/(vect(iter_cp+3)-vect(iter_cp))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp)/(vect(iter_cp+2)-vect(iter_cp)) - &
                                                (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+2)-vect(iter_cp))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1))) & 
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                + (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1)))) &
                                                - (eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+3)-vect(iter_cp))* &
                                                 (vect(iter_cp+2)-vect(iter_cp+1)))* &
                                                (2/(vect(iter_cp+2)-vect(iter_cp)) + 2/(vect(iter_cp+3)-vect(iter_cp+1))) &
                                                + 2/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) + &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+2)-vect(iter_cp+1)))) + &
                                                2*(vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+1)) * (vect(iter_cp+2)-vect(iter_cp+1))) 
                else if(eval_points(iter_pos)>=vect(iter_cp+2).AND.eval_points(iter_pos)<vect(iter_cp+3)) then
                                                 CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) = &
                                                 2/(vect(iter_cp+3)-vect(iter_cp))* &   
                                                 (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) &
                                                - (vect(iter_cp+3)-eval_points(iter_pos))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                 (vect(iter_cp+3)-vect(iter_cp+2)))) &
                                                + 2*(eval_points(iter_pos)-vect(iter_cp))/((vect(iter_cp+3)-vect(iter_cp))* &
                                                (vect(iter_cp+3)-vect(iter_cp+1))*(vect(iter_cp+3)-vect(iter_cp+2))) &
                                                - 2/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+3)-vect(iter_cp+1)) + &
                                                (eval_points(iter_pos)-vect(iter_cp+1))/((vect(iter_cp+3)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2))) &
                                                - CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+2)) + &
                                                + (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2)))) &
                                                - (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (vect(iter_cp+3)-vect(iter_cp+2)))* &
                                                (2/(vect(iter_cp+3)-vect(iter_cp+1)) + 2/(vect(iter_cp+4)-vect(iter_cp+2)))
                else if(eval_points(iter_pos)>=vect(iter_cp+3).AND.eval_points(iter_pos)<vect(iter_cp+4)) then 
                                                CPCoil(icoil)%db_dt_2(iter_pos,iter_cp) =  &
                                                 -2/(vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (-CPCoil(icoil)%basis_1(iter_pos,iter_cp+1)/(vect(iter_cp+4)-vect(iter_cp+2)) - &
                                                (vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+2))* &
                                                (vect(iter_cp+4)-vect(iter_cp+3)))) &
                                                + 2*(vect(iter_cp+4)-eval_points(iter_pos))/((vect(iter_cp+4)-vect(iter_cp+1))* &
                                                (vect(iter_cp+4)-vect(iter_cp+2))*(vect(iter_cp+4)-vect(iter_cp+3)))
                endif
        enddo
    enddo


    return 
END SUBROUTINE eval_basis2


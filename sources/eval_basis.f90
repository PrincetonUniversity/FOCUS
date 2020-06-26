SUBROUTINE eval_basis(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0: coil(icoil)%NS)
    integer :: iter_cp,iter_pos
    type(SplineCoil) :: Coil_temp

    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points
    Coil_temp = CPCoil(icoil)

    do iter_cp=0,Coil_temp%NCP+2
        do iter_pos=0,N 
                if(eval_points(iter_pos)>=Coil_temp%vect(iter_cp).AND.eval_points(iter_pos)<Coil_temp%vect(iter_cp+1))&
                 Coil_temp%basis_0(iter_pos,iter_cp) = 1
        enddo
    enddo

    do iter_cp=0,Coil_temp%NCP+1   
        do iter_pos=0,N 
                if( Coil_temp%basis_0(iter_pos,iter_cp) /= 0 .OR. Coil_temp%basis_0(iter_pos,iter_cp+1) /= 0)&
                         Coil_temp%basis_1(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - Coil_temp%vect(iter_cp)) / &
                        (Coil_temp%vect(iter_cp+1) - Coil_temp%vect(iter_cp))*Coil_temp%basis_0(iter_pos,iter_cp) &
                         + 1.0*( - eval_points(iter_pos) + Coil_temp%vect(iter_cp+2))/ &
                         (Coil_temp%vect(iter_cp+2) - Coil_temp%vect(iter_cp+1))*Coil_temp%basis_0(iter_pos,iter_cp+1)
                enddo    
        enddo

    do iter_cp=0,Coil_temp%NCP
        do iter_pos=0,N 
                if ( (Coil_temp%basis_1(iter_pos,iter_cp) .NE. 0) .OR. (Coil_temp%basis_1(iter_pos,iter_cp+1) .NE. 0)) &
                         Coil_temp%basis_2(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - Coil_temp%vect(iter_cp))/ & 
                        (Coil_temp%vect(iter_cp+2) - Coil_temp%vect(iter_cp))*Coil_temp%basis_1(iter_pos,iter_cp) &
                        + 1.0*( - eval_points(iter_pos) + Coil_temp%vect(iter_cp+3))/ &
                        (Coil_temp%vect(iter_cp+3) - Coil_temp%vect(iter_cp+1))*Coil_temp%basis_1(iter_pos,iter_cp+1)
                enddo    
        enddo
    
    do iter_cp=0,Coil_temp%NCP-1
        do iter_pos=0,N 
                if ( (Coil_temp%basis_2(iter_pos,iter_cp) .NE. 0) .OR. (Coil_temp%basis_2(iter_pos,iter_cp+1) .NE. 0))&
                         Coil_temp%basis_3(iter_pos,iter_cp) = 1.0*(eval_points(iter_pos) - Coil_temp%vect(iter_cp))/ &
                        (Coil_temp%vect(iter_cp+3) - Coil_temp%vect(iter_cp))*Coil_temp%basis_2(iter_pos,iter_cp) &
                         + 1.0*( - eval_points(iter_pos) + Coil_temp%vect(iter_cp+4))/ &
                        (Coil_temp%vect(iter_cp+4) - Coil_temp%vect(iter_cp+1))*Coil_temp%basis_2(iter_pos,iter_cp+1)
                enddo
        enddo


    return 
END SUBROUTINE eval_basis

SUBROUTINE eval_basis1(icoil)
    use globals  
    implicit none
    integer :: N,i,j,icoil
    REAL :: eval_points(0:N)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT)
       
    vect = CPCoil(icoil)%vect
    N = CPCoil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points

    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N 
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
    REAL :: eval_points(0:N)
    integer :: iter_cp,iter_pos
    REAL :: vect(0:CPCoil(icoil)%NT)

    vect = CPCoil(icoil)%vect(iter_cp)
    N = coil(icoil)%NS
    eval_points = CPCoil(icoil)%eval_points
    
    do iter_cp=0,CPCoil(icoil)%NCP-1
        do iter_pos=0,N 
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


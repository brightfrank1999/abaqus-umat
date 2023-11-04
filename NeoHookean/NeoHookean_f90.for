!************************************************************************
!
! UMAT for Neo-Hookean Hyperelastic Materials 
!
! Fangda Cui, May 2023, Implemented in Abaqus 2021
!
!************************************************************************
!    
!     Two Material Properties
!     Gshear, Kbulk 
!
!**********************************************************************
subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
& RPL,DDSDDT,DRPLDE,DRPLDT,&
& STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
& NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
& CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    include 'ABA_PARAM.INC'

    dimension STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),&
	 & DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),& 
	 & TIME(2), PREDEF(1), DPRED(1), PROPS(NPROPS),COORDS(3),&
	 & DROT(3,3),DFGRD0(3,3),DFGRD1(3,3), JSTEP(4)
	 
	!--DEFINE VARIABLES
    character*80 CMNAME, FILE1
    character*256 JOBNAME, OUTDIR, FILENAME

    integer I, J, K , L

	real*8 IDEN(3,3),F_T(3,3),F_TAU(3,3), DETF
     
    real*8 GSHEAR            !- shear modulus 
	real*8 KBULK             !- bulk modulus 
    real*8 BDIS(3,3)         
	real*8 TRBDIS
    real*8 BDIS0(3,3)      
    real*8 T_TAU(3,3)        !- Cauchy stress tensor 
	real*8 SPTANMOD(3,3,3,3) !- Spatial tangent modulus 

    !--DEFINE PARAMETERS
    real*8 ZERO,ONE,TWO,HALF,THREE,THIRD,NINE
    parameter(ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0,THREE=3.D0,THIRD=1.D0/3.D0,NINE=9.D0)
	
    !--IDENTITY MATRIX   
    call ONEM(IDEN)

    !--OBTAIN OLD AND NEW DEFORMATION GRADIENTS      
    F_T = DFGRD0
    F_TAU = DFGRD1

    !--OBTAIN MATERIAL PROPERTIES USED FOR INIT      
    GSHEAR = PROPS(1)
    KBULK  = PROPS(2)

    !--Calculate determinant of new deformation gradient
    call MDET(F_TAU,DETF)

    !--COMPUTE THE TOTAL DISTORTIONAL LEFT CAUCHY-GREEN TENSOR AND ITS DEVIATOR   
    BDIS = (DETF**(-TWO/THREE))*MATMUL(F_TAU,TRANSPOSE(F_TAU))
    TRBDIS = BDIS(1,1) + BDIS(2,2) + BDIS(3,3)
    BDIS0 = BDIS - THIRD*TRBDIS*IDEN 
      
    !--COMPUTE THE CAUCHY STRESS 
	T_TAU = (GSHEAR*BDIS0 + KBULK*(DETF-1)*IDEN)/DETF

	!--COMPUTE THE MATERIAL JACOBIAN 
    SPTANMOD = ZERO
    do I=1,3
		do J=1,3
            do K=1,3
                do L=1,3
                    SPTANMOD(I,J,K,L) = SPTANMOD(I,J,K,L)&
                       &+ (GSHEAR/DETF)*(HALF*IDEN(I,K)*BDIS(J,L)&
                       &+ HALF*IDEN(J,K)*BDIS(I,L)&
                       &+ HALF*IDEN(I,L)*BDIS(J,K)&
                       &+ HALF*IDEN(J,L)*BDIS(I,K)&
                       &- (TWO/THREE)*IDEN(I,J)*BDIS(K,L)&
                       &- (TWO/THREE)*IDEN(K,L)*BDIS(I,J)&
                       &+ (TWO/NINE)*TRBDIS*IDEN(I,J)*IDEN(K,L))&
                       &+ (KBULK/DETF)*IDEN(I,J)*IDEN(K,L)
                enddo 
            enddo 
        enddo 
    enddo

	!-- Return ABAQUS/STANDARD the Cauchy stress
	if(NTENS.eq.6) then 
		STRESS(1) = T_TAU(1,1)
		STRESS(2) = T_TAU(2,2)
		STRESS(3) = T_TAU(3,3)
		STRESS(4) = T_TAU(1,2)
		STRESS(5) = T_TAU(1,3)
		STRESS(6) = T_TAU(2,3)
	elseif(NTENS.eq.4) then
		STRESS(1) = T_TAU(1,1)
		STRESS(2) = T_TAU(2,2)
		STRESS(3) = T_TAU(3,3)
		STRESS(4) = T_TAU(1,2)
	endif 

	!--RETURN ABAQUS/STANDARD THE STRESS-DEFORMATION JACOBIAN     
    if(NTENS.eq.6) then 
		call JAC3D(SPTANMOD,DDSDDE)
	elseif(NTENS.eq.4) then  	
		call JAC2D(SPTANMOD,DDSDDE)
	endif 

return
end subroutine UMAT

subroutine JAC2D(SPTANMOD,DDSDDE)

	real*8 SPTANMOD(3,3,3,3),DDSDDE(4,4)

    DDSDDE(1,1) = SPTANMOD(1,1,1,1)
    DDSDDE(1,2) = SPTANMOD(1,1,2,2)
    DDSDDE(1,3) = SPTANMOD(1,1,3,3)
    DDSDDE(1,4) = SPTANMOD(1,1,1,2)

    DDSDDE(2,1) = SPTANMOD(2,2,1,1)
    DDSDDE(2,2) = SPTANMOD(2,2,2,2)
    DDSDDE(2,3) = SPTANMOD(2,2,3,3)
    DDSDDE(2,4) = SPTANMOD(2,2,1,2)

    DDSDDE(3,1) = SPTANMOD(3,3,1,1)
    DDSDDE(3,2) = SPTANMOD(3,3,2,2)
    DDSDDE(3,3) = SPTANMOD(3,3,3,3)
    DDSDDE(3,4) = SPTANMOD(3,3,1,2)

    DDSDDE(4,1) = SPTANMOD(1,2,1,1)
    DDSDDE(4,2) = SPTANMOD(1,2,2,2)
    DDSDDE(4,3) = SPTANMOD(1,2,3,3)
    DDSDDE(4,4) = SPTANMOD(1,2,1,2)

end subroutine JAC2D
	  
subroutine JAC3D(SPTANMOD,DDSDDE)

    real*8 SPTANMOD(3,3,3,3), DDSDDE(6,6)

    DDSDDE(1,1) = SPTANMOD(1,1,1,1)
    DDSDDE(1,2) = SPTANMOD(1,1,2,2)
    DDSDDE(1,3) = SPTANMOD(1,1,3,3)
    DDSDDE(1,4) = SPTANMOD(1,1,1,2)
    DDSDDE(1,5) = SPTANMOD(1,1,1,3)
    DDSDDE(1,6) = SPTANMOD(1,1,2,3)

    DDSDDE(2,1) = SPTANMOD(2,2,1,1)
    DDSDDE(2,2) = SPTANMOD(2,2,2,2)
    DDSDDE(2,3) = SPTANMOD(2,2,3,3)
    DDSDDE(2,4) = SPTANMOD(2,2,1,2)
    DDSDDE(2,5) = SPTANMOD(2,2,1,3)
    DDSDDE(2,6) = SPTANMOD(2,2,2,3)

    DDSDDE(3,1) = SPTANMOD(3,3,1,1)
    DDSDDE(3,2) = SPTANMOD(3,3,2,2)
    DDSDDE(3,3) = SPTANMOD(3,3,3,3)
    DDSDDE(3,4) = SPTANMOD(3,3,1,2)
    DDSDDE(3,5) = SPTANMOD(3,3,1,3)
    DDSDDE(3,6) = SPTANMOD(3,3,2,3)

    DDSDDE(4,1) = SPTANMOD(1,2,1,1)
    DDSDDE(4,2) = SPTANMOD(1,2,2,2)
    DDSDDE(4,3) = SPTANMOD(1,2,3,3)
    DDSDDE(4,4) = SPTANMOD(1,2,1,2)
    DDSDDE(4,5) = SPTANMOD(1,2,1,3)
    DDSDDE(4,6) = SPTANMOD(1,2,2,3)

    DDSDDE(5,1) = SPTANMOD(1,3,1,1)
    DDSDDE(5,2) = SPTANMOD(1,3,2,2)
    DDSDDE(5,3) = SPTANMOD(1,3,3,3)
    DDSDDE(5,4) = SPTANMOD(1,3,1,2)
    DDSDDE(5,5) = SPTANMOD(1,3,1,3)
    DDSDDE(5,6) = SPTANMOD(1,3,2,3)

    DDSDDE(6,1) = SPTANMOD(2,3,1,1)
    DDSDDE(6,2) = SPTANMOD(2,3,2,2)
    DDSDDE(6,3) = SPTANMOD(2,3,3,3)
    DDSDDE(6,4) = SPTANMOD(2,3,1,2)
    DDSDDE(6,5) = SPTANMOD(2,3,1,3)
    DDSDDE(6,6) = SPTANMOD(2,3,2,3)

end subroutine JAC3D

subroutine ONEM(A)
!**********************************************************************
!	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
!	3 BY 3 MATRIX [A]
!**********************************************************************
	REAL*8 A(3,3)

	dO I=1,3
		do J=1,3
			IF (I .EQ. J) THEN
				A(I,J) = 1.d0
			ELSE
				A(I,J) = 0.d0
			ENDIF
        end do 
	end do

return 
end subroutine ONEM

subroutine MDET(A,DET)
!**********************************************************************
! 	THIS SUBROUTINE CALCULATES THE DETERMINANT
! 	OF A 3 BY 3 MATRIX [A].
!**********************************************************************

	REAL*8  A(3,3), DET

	DET = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2)&
     	  &- A(3,1)*A(2,2)*A(1,3) - A(3,2)*A(2,3)*A(1,1) - A(3,3)*A(2,1)*A(1,2)

return 
end subroutine MDET

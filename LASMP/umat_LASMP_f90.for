!************************************************************************
!
! UMAT for Shape Light Acitvated Shape Memory Polymer
!
! Fangda Cui, Dec 2014, Implemented in Abaqus 6.13-1
!
!************************************************************************
!     State Variables
!     --------------------------------------------------------------
!     Twenty State Variables In Total
!
!     Twelve Material Properties
!**********************************************************************
subroutine UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,&
& RPL,DDSDDT,DRPLDE,DRPLDT,&
& STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,&
& NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,&
& ELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

    include 'ABA_PARAM.INC'

    dimension STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS),&
	 & DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),&
	 & TIME(2), PREDEF(1),DPRED(1), PROPS(NPROPS), COORDS(3),&
	 & DROT(3,3), DFGRD0(3,3), DFGRD1(3,3), JSTEP(4)

    !--DEFINE VARIABLES
    character*80 CMNAME

    integer I, J, K, L

    real*8 IDEN(3,3), F_T(3,3), F_TAU(3,3)
	real*8 FINV1(3,3), DETF
    real*8 F_L1(3,3), DETF_L1
    
    real*8 ALPHA,ALPHA_T,ALPHA_TAU,KF1,KR1
	
    real*8 GSHEAR, KBULK
    real*8 BDIS(3,3), TRBDIS, BDIS0(3,3)
    
    real*8 GSHEAR_L1, KBULK_L1
    real*8 BDIS_L1(3,3), TRBDIS_L1, BDIS0_L1(3,3)
    
    real*8 TA_TAU(3,3), TB_TAU(3,3), T_TAU(3,3), SPTANMOD(3,3,3,3)
	
    !--DEFINE PARAMETERS
    real*8 ZERO,ONE,TWO,HALF,THREE,THIRD,NINE
    PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0,THREE=3.D0,THIRD=1.D0/3.D0,NINE=9.D0)


    !--IDENTITY MATRIX
	CALL ONEM(IDEN)

    !--OBTAIN OLD AND NEW DEFORMATION GRADIENTS    
    F_T = DFGRD0
    F_TAU = DFGRD1
      
    !--OBTAIN MATERIAL PROPERTIES USED FOR INIT
    GSHEAR = PROPS(1)
    KBULK  = PROPS(2)
	
    GSHEAR_L1 = PROPS(3)
    KBULK_L1 = PROPS(4)
    ALPHA = PROPS(5)  
    KF1 = PROPS(6)
    KR1 = PROPS(7)
   
    !--INITIALIZE STATE VARIABLES AND NORMAL VARIABLES
    IF((TIME(1).EQ.0.0).and.(TIME(2).EQ.0.0)) THEN
		STATEV(1) = 1.0
		STATEV(2) = 0.0
		STATEV(3) = 0.0
		STATEV(4) = 0.0
		STATEV(5) = 1.0
		STATEV(6) = 0.0
		STATEV(7) = 0.0
		STATEV(8) = 0.0
		STATEV(9) = 1.0
		STATEV(10) = 0.0
    ELSEIF((KSTEP.EQ.2).and.(TIME(1).EQ.0.0)) THEN
    !--COMPUTE THE REFERENCE DEFORMATIONG GRADIENT
        CALL M3INV(F_TAU,FINV1)
        STATEV(1) = FINV1(1,1)
        STATEV(2) = FINV1(1,2)
        STATEV(3) = FINV1(1,3)
        STATEV(4) = FINV1(2,1)
        STATEV(5) = FINV1(2,2)
        STATEV(6) = FINV1(2,3)
        STATEV(7) = FINV1(3,1)
        STATEV(8) = FINV1(3,2)
        STATEV(9) = FINV1(3,3)
    ENDIF
	
    !--OBTAIN THE OLD STATE VARIABLE
    FINV1(1,1) = STATEV(1)
    FINV1(1,2) = STATEV(2)
    FINV1(1,3) = STATEV(3)
    FINV1(2,1) = STATEV(4)
    FINV1(2,2) = STATEV(5)
    FINV1(2,3) = STATEV(6)
    FINV1(3,1) = STATEV(7)
    FINV1(3,2) = STATEV(8)
    FINV1(3,3) = STATEV(9)
    ALPHA_T = STATEV(10)
	    
    IF(KSTEP.EQ.1) THEN
        ALPHA_TAU = 0.0
!        write(*,*) 'ALPHA_TAU=', ALPHA_TAU
    ENDIF
	
    !--Evolution equation I  
    IF(KSTEP.EQ.2) THEN
        ALPHA_TAU = ALPHA_T + DTIME*KF1*(ONE-ALPHA_T)**TWO
        write(*,*) 'ALPHA_TAU=',ALPHA_TAU
    ELSEIF(KSTEP.EQ.4) THEN
        ALPHA_TAU = ALPHA_T - DTIME*KR1*ALPHA_T
    ELSE
        ALPHA_TAU = ALPHA_T
    ENDIF
	
    IF(ALPHA_TAU.GT.0.8) THEN
        ALPHA_TAU = 0.8
    ELSEIF(ALPHA_TAU.LT.0) THEN
        ALPHA_TAU = 0.d0
    ENDIF  
      
    !--PRINT OUT THE FLAG       
!    if((KSTEP.EQ.1).AND.(TIME(1).EQ.0))THEN
!         print*, 'NOW WE ARE IN STEP ONE - Loading!\n'
!    endif
!    if((KSTEP.EQ.2).AND.(TIME(1).EQ.0))THEN
!         print*, 'NOW WE ARE IN STEP TWO - Xlinking!\n'
!    endif
!    if((KSTEP.EQ.3).AND.(TIME(1).EQ.0))THEN
!         print*, 'NOW WE ARE IN STEP THREE - Unloading!\n'
!    endif
!    if((KSTEP.EQ.4).AND.(TIME(1).EQ.0))THEN
!         print*, 'NOW WE ARE IN STEP FOUR - Cleavage!\n'
!    endif  
	
    !--COMPUTE DEFORMATION GRADIENT FOR THE FIRST LIGHT INDUCED NETWORK
    F_L1 = MATMUL(F_TAU,FINV1)
    CALL MDET(F_L1,DETF_L1) 
      
    !--COMPUTE THE RELATIVE VOLUME CHANGE
    CALL MDET(F_TAU,DETF)
	
    !--COMPUTE THE TOTAL DISTORTIONAL LEFT CAUCHY-GREEN TENSOR AND ITS DEVIATOR      
    BDIS = (DETF**(-TWO/THREE))*MATMUL(F_TAU,TRANSPOSE(F_TAU))
    TRBDIS = BDIS(1,1) + BDIS(2,2) + BDIS(3,3)
    BDIS0 = BDIS - THIRD*TRBDIS*IDEN 
       
    !--COMPUTE THE TOTAL DISTORTIONAL LEFT CAUCHY-GREEN TENSOR AND ITS DEVIATOR FOR THE FIRST LIGHT INDUCED NETWORK
    BDIS_L1 = (DETF_L1**(-TWO/THREE))*MATMUL(F_L1,TRANSPOSE(F_L1))
    TRBDIS_L1 = BDIS_L1(1,1) + BDIS_L1(2,2) + BDIS_L1(3,3)
    BDIS0_L1 = BDIS_L1 - THIRD*TRBDIS_L1*IDEN 
	
    !--COMPUTE THE CAUCHY STRESS 
    TA_TAU = (GSHEAR*BDIS0 + KBULK*(DETF-1)*IDEN)/DETF
    TB_TAU = (GSHEAR_L1*BDIS0_L1 + KBULK_L1*(DETF_L1-1)*IDEN)/DETF_L1   
    T_TAU = TA_TAU + ALPHA_TAU*TB_TAU
          
	!--COMPUTE THE MATERIAL JACOBIAN      
	SPTANMOD = ZERO
	DO I=1,3
		DO J=1,3
            DO K=1,3
				DO L=1,3
					SPTANMOD(I,J,K,L) = SPTANMOD(I,J,K,L)&
                      &+ (GSHEAR/DETF)*&
                      &(&
                      &HALF*IDEN(I,K)*BDIS(J,L)&
                      &+ HALF*IDEN(J,K)*BDIS(I,L)&
                      &+ HALF*IDEN(I,L)*BDIS(J,K)&
                      &+ HALF*IDEN(J,L)*BDIS(I,K)&
                      &- (TWO/THREE)*IDEN(I,J)*BDIS(K,L)&
                      &- (TWO/THREE)*IDEN(K,L)*BDIS(I,J)&
                      &+ (TWO/NINE)*TRBDIS*IDEN(I,J)*IDEN(K,L)&
                      &)&
                      &+ (KBULK/DETF)*IDEN(I,J)*IDEN(K,L)&
                      &+ ALPHA_TAU*(GSHEAR_L1/DETF_L1)*&
                      &(&
                      &HALF*IDEN(I,K)*BDIS_L1(J,L)&
                      &+ HALF*IDEN(J,K)*BDIS_L1(I,L)&
                      &+ HALF*IDEN(I,L)*BDIS_L1(J,K)&
                      &+ HALF*IDEN(J,L)*BDIS_L1(I,K)&
                      &- (TWO/THREE)*IDEN(I,J)*BDIS_L1(K,L)&
                      &- (TWO/THREE)*IDEN(K,L)*BDIS_L1(I,J)&
                      &+ (TWO/NINE)*TRBDIS_L1*IDEN(I,J)*IDEN(K,L)&
                      &)&
                      &+ ALPHA_TAU*(KBULK_L1/DETF_L1)*IDEN(I,J)*IDEN(K,L)
                ENDDO
            ENDDO
        ENDDO
    ENDDO
      
	 
    !--UPDATE STATE VARIABLES AT THE END OF THE INCREMENT
    STATEV(1) = FINV1(1,1)
    STATEV(2) = FINV1(1,2)
    STATEV(3) = FINV1(1,3)
    STATEV(4) = FINV1(2,1)
    STATEV(5) = FINV1(2,2)
    STATEV(6) = FINV1(2,3)
    STATEV(7) = FINV1(3,1)
    STATEV(8) = FINV1(3,2)
    STATEV(9) = FINV1(3,3)
    STATEV(10) = ALPHA_TAU
	
	
    !--RETURN ABAQUS/STANDARD THE CAUCHY STRESS
    IF(NTENS.EQ.6) THEN   
		!--3D PROBLEM
        STRESS(1) = T_TAU(1,1)
        STRESS(2) = T_TAU(2,2)
        STRESS(3) = T_TAU(3,3)
        STRESS(4) = T_TAU(1,2)
        STRESS(5) = T_TAU(1,3)
        STRESS(6) = T_TAU(2,3)
    ELSEIF(NTENS.EQ.4) THEN
        !--2D PROBLEM
        STRESS(1) = T_TAU(1,1)
        STRESS(2) = T_TAU(2,2)
        STRESS(3) = T_TAU(3,3)
        STRESS(4) = T_TAU(1,2)
    ENDIF
	
	
    !--RETURN ABAQUS/STANDARD THE STRESS-DEFORMATION JACOBIAN     
    IF(NTENS.EQ.6) THEN
       CALL JAC3D(SPTANMOD,DDSDDE)
    ELSEIF(NTENS.EQ.4) THEN
       CALL JAC2D(SPTANMOD,DDSDDE)
    ENDIF


return
end subroutine UMAT


SUBROUTINE JAC2D(SPTANMOD,DDSDDE)

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

END SUBROUTINE JAC2D


subroutine JAC3D(SPTANMOD, DDSDDE)

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

!**********************************************************************
!**********************************************************************
!	THE FOLLOWING SUBROUTINES ARE UTILITY ROUTINES
!**********************************************************************
!**********************************************************************
subroutine ONEM(A)
!**********************************************************************
!	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
!	3 BY 3 MATRIX [A]
!**********************************************************************

    real*8 A(3,3)
    data ZERO/0.D0/
    data ONE/1.D0/

	do I=1,3
		do J=1,3
			if (I .EQ. J) then
				A(I,J) = 1.0
			else 
				A(I,J) = 0.0
            endif
		enddo
	enddo

return 
end subroutine ONEM 


subroutine MDET(A, DET)
!********************************************************************** 
! 	THIS SUBROUTINE CALCULATES THE DETERMINANT
! 	OF A 3 BY 3 MATRIX [A].
!**********************************************************************

	real*8  A(3,3), DET

	DET = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
		&+ A(1,3)*A(2,1)*A(3,2) - A(3,1)*A(2,2)*A(1,3)&
    	&- A(3,2)*A(2,3)*A(1,1) - A(3,3)*A(2,1)*A(1,2)

return
end subroutine MDET


subroutine MCOFAC(A, ACOFAC)
!**********************************************************************
! 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
! 	AND PLACES THE RESULT IN [ACOFAC]. 
!**********************************************************************

	real*8  A(3,3), ACOFAC(3,3)

	ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
	ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
	ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
	ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
	ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
	ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
	ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
	ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
	ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

return 
end subroutine MCOFAC


subroutine MTRANS(A, ATRANS)
!********************************************************************** 
!	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
!	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
!**********************************************************************
	integer I, J
	real*8 A(3,3), ATRANS(3,3)

	do I=1,3
		do J=1,3
			ATRANS(J,I) = A(I,J)
		enddo 
	enddo 

return
end subroutine MTRANS


subroutine M3INV(A, AINV)
!**********************************************************************
! 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
!	[A] AND PLACES THE RESULT IN [AINV]. 
! 	IF DET(A) IS ZERO, THE CALCULATION
! 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
!**********************************************************************
	integer I, J 
	real*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

!	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
!	DET		-- THE COMPUTED DETERMINANT OF [A].
!	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
!			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
!			   IS CALLED THE COFACTOR OF A(I,J).
!	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
!			   OBTAINED BY REPLACING EACH ELEMENT OF
!			   [A] BY ITS COFACTOR, AND THEN TAKING
!			   TRANSPOSE OF THE RESULTING MATRIX.
!	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
!			   [AINV] = [AADJ]/DET.
!----------------------------------------------------------------------

	call MDET(A, DET)
	if ((DET - 0.D0).le.0.0000001d0) then 
		print *, "The determinant of the matrix is zero!"
		pause
	endif 
	call MCOFAC(A, ACOFAC)
	call MTRANS(ACOFAC, AADJ)
	do I = 1,3
		do J = 1,3
			AINV(I,J) = AADJ(I,J)/DET
		enddo 
	enddo 

return 
end subroutine M3INV

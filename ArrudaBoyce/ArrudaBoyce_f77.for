C************************************************************************
C
C UMAT for ArrudaBoyce
C
C Fangda Cui, May 2015, Implemented in Abaqus 6.13-1
C
C************************************************************************   
C     Material Properties Vector
C     --------------------------------------------------------------
C     C1 = props(1)
C     Lamda_L = props(2)
C     D1 = props(3) 
C     
C     ShearModules = 2*C1 
C     Lamda_L= sqrt(N)       
C     BulkModules = 2*D1
C
C**********************************************************************

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     + RPL,DDSDDT,DRPLDE,DRPLDT,
     + STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     + NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     + CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)

	  INCLUDE 'ABA_PARAM.INC'

	  DIMENSION STRESS(NTENS),STATEV(NSTATV),
     + DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     + STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     + PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

C
C     DEFINE VARIABLES
C
	  CHARACTER*80 CMNAME,FILE1
	  CHARACTER*256 JOBNAME,OUTDIR,FILENAME

	  INTEGER I,J,K,L,ITERERROR,LENJOBNAME,LENOUTDIR

	  REAL*8 IDEN(3,3),F_T(3,3),F_TAU(3,3),THETA_T,THETA_TAU,T_TAU(3,3)
	  REAL*8 T_A(3,3),  T_B(3,3)
	  REAL*8 GSHEAR,EFFSTR,DETF,DTRDF(3,3,3,3),SPTANMOD(3,3,3,3),TRB
	  REAL*8 FINV(3,3),BDIS(3,3),TRBDIS,BDIS0(3,3),B_TAU(3,3),LAMDA_CDIS
      REAL*8 KBULK,LAMDA_L
C
C     DEFINE PARAMETERS
C
      REAL*8 ZERO,ONE,TWO,HALF,THREE,THIRD,NINE
      PARAMETER(ZERO=0.D0,ONE=1.D0,TWO=2.D0,HALF=0.5D0,THREE=3.D0,
     +     FOUR=4.D0,THIRD=1.D0/3.D0,NINE=9.D0)
C
C     IDENTITY MATRIX
C      
      CALL ONEM(IDEN)
C
C     OBTAIN OLD AND NEW DEFORMATION GRADIENTS
C      
      F_T = DFGRD0
      F_TAU = DFGRD1
C      
C     OBTAIN OLD AND NEW TEMPERATURE
C
      THETA_T = TEMP
      THETA_TAU = TEMP + DTEMP
C	  print*, F_T
C      
C     OBTAIN MATERIAL PROPERTIES USED FOR INIT
C      
      C1 = PROPS(1)
      LAMDA_L = PROPS(2)
      D1 = PROPS(3)
      GSHEAR=TWO*C1
      KBULK=TWO*D1
C      PRINT*, GSHEAR
C	  PRINT*, C1
C      
C     COMPUTE THE RELATIVE VOLUME CHANGE
C      
      CALL MDET(F_TAU,DETF)
C
C     COMPUTE THE INVERSE OF THE DEFORMATION GRADIENT
C    
      CALL M3INV(F_TAU,FINV)
C      
C     COMPUTE THE TOTAL DISTORTIONAL LEFT CAUCHY-GREEN TENSOR
C     AND ITS DEVIATOR
C      
      B_TAU = MATMUL(F_TAU,TRANSPOSE(F_TAU))
      TRB = B_TAU(1,1) + B_TAU(2,2)+B_TAU(3,3)
      BDIS = (DETF**(-TWO/THREE))*B_TAU
      TRBDIS = BDIS(1,1) + BDIS(2,2) + BDIS(3,3)
      BDIS0 = BDIS - THIRD*TRBDIS*IDEN
C
C     COMPUTE DISTORTIONAL LAMDA_C
C
      LAMDA_CDIS=SQRT((ONE/THREE)*TRBDIS)
C      
C     COMPUTE THE CAUCHY STRESS 
C
	  T_A1 = THREE-(LAMDA_CDIS/LAMDA_L)**TWO
	  T_A2 = ONE-(LAMDA_CDIS/LAMDA_L)**TWO
      T_A = (ONE/THREE)*GSHEAR*(T_A1)/(T_A2)*BDIS0
	  T_B = KBULK*DLOG(DETF)*IDEN
      T_TAU = (ONE/DETF)*(T_A + T_B)
        
C	  PRINT*, T_TAU(3,3)
C      
C     COMPUTE THE MATERIAL JACOBIAN
C      
      SPTANMOD = ZERO
      DO I=1,3
	    DO J=1,3
		  DO K=1,3
		    DO L=1,3
			  SPTANMOD(I,J,K,L) = SPTANMOD(I,J,K,L)
     +            + (ONE/THREE)*(GSHEAR/DETF)*
     +              (
     +              (FOUR/THREE)*(ONE/(LAMDA_L)**2)*BDIS(K,L)
     +              -(FOUR/NINE)*(ONE/(LAMDA_L)**2)*TRBDIS*IDEN(L,K)
     +              )/
     +              (1-(LAMDA_CDIS/LAMDA_L)**2)**2*
     +              (
     +              BDIS(I,J)-(ONE/THREE)*TRBDIS*IDEN(I,J)
     +              )
C                  
     +              + (ONE/THREE)*(GSHEAR/DETF)*
     +              (
     +              HALF*IDEN(I,K)*BDIS(J,L)
     +              + HALF*IDEN(J,K)*BDIS(I,L)
     +              + HALF*IDEN(I,L)*BDIS(J,K)
     +              + HALF*IDEN(J,L)*BDIS(I,K)
     +              - (TWO/THREE)*IDEN(I,J)*BDIS(K,L)
     +              - (TWO/THREE)*IDEN(K,L)*BDIS(I,J)
     +              + (TWO/NINE)*TRBDIS*IDEN(I,J)*IDEN(K,L)
     +              )*
     +              (3-(LAMDA_CDIS/LAMDA_L)**2)/
     +              (1-(LAMDA_CDIS/LAMDA_L)**2)
C                  
     +              + (KBULK/DETF)*IDEN(I,J)*IDEN(K,L)
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C     ABAQUS/STANDARD THE CAUCHY STRESS
      
      IF(NTENS.EQ.6) THEN
C
C       3D PROBLEM
C
	    STRESS(1) = T_TAU(1,1)
		STRESS(2) = T_TAU(2,2)
		STRESS(3) = T_TAU(3,3)
		STRESS(4) = T_TAU(1,2)
		STRESS(5) = T_TAU(1,3)
		STRESS(6) = T_TAU(2,3)
      ELSEIF(NTENS.EQ.4) THEN
C     
C       2D PROBLEM
C         
		STRESS(1) = T_TAU(1,1)
		STRESS(2) = T_TAU(2,2)
		STRESS(3) = T_TAU(3,3)
		STRESS(4) = T_TAU(1,2)
	  ENDIF

C     ABAQUS/STANDARD THE STRESS-DEFORMATION JACOBIAN
      IF(NTENS.EQ.6) THEN
		CALL JAC3D(SPTANMOD,DDSDDE)
      ELSEIF(NTENS.EQ.4) THEN
		CALL JAC2D(SPTANMOD,DDSDDE)
      ENDIF

      RETURN
      END

C***********************************************************************

      SUBROUTINE JAC2D(SPTANMOD,DDSDDE)
	  
	  REAL*8 SPTANMOD(3,3,3,3),DDSDDE(4,4)
	  
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
	  
	  RETURN
      END 

C***********************************************************************

      SUBROUTINE JAC3D(SPTANMOD,DDSDDE)
	  
	  REAL*8 SPTANMOD(3,3,3,3),DDSDDE(6,6)
	  
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
	  
	  RETURN
      END 


C**********************************************************************
	  SUBROUTINE ONEM(A)

C	THIS SUBROUTINE STORES THE IDENTITY MATRIX IN THE 
C	3 BY 3 MATRIX [A]
C**********************************************************************

      REAL*8 A(3,3)
      DATA ZERO/0.D0/
      DATA ONE/1.D0/

	  DO 1 I=1,3
	    DO 1 J=1,3
	      IF (I .EQ. J) THEN
                A(I,J) = 1.0
              ELSE
                A(I,J) = 0.0
              ENDIF
1         CONTINUE

	  RETURN
	  END

C**********************************************************************
	  SUBROUTINE MTRANS(A,ATRANS)
 
C	THIS SUBROUTINE CALCULATES THE TRANSPOSE OF AN 3 BY 3 
C	MATRIX [A], AND PLACES THE RESULT IN ATRANS. 
C**********************************************************************

		REAL*8 A(3,3),ATRANS(3,3)

		DO 1 I=1,3
			DO 1 J=1,3
				ATRANS(J,I) = A(I,J)
1			CONTINUE
	  RETURN
	  END

C**********************************************************************
	  SUBROUTINE MDET(A,DET)
 
C 	THIS SUBROUTINE CALCULATES THE DETERMINANT
C 	OF A 3 BY 3 MATRIX [A].
C**********************************************************************

		REAL*8  A(3,3), DET

	  DET =	  A(1,1)*A(2,2)*A(3,3) 
     +	        + A(1,2)*A(2,3)*A(3,1)
     +	        + A(1,3)*A(2,1)*A(3,2)
     +		- A(3,1)*A(2,2)*A(1,3)
     +		- A(3,2)*A(2,3)*A(1,1)
     +		- A(3,3)*A(2,1)*A(1,2)

	  RETURN
	  END

C**********************************************************************
	  SUBROUTINE M3INV(A,AINV)

C 	THIS SUBROUTINE CALCULATES THE THE INVERSE OF A 3 BY 3 MATRIX
C	[A] AND PLACES THE RESULT IN [AINV]. 
C 	IF DET(A) IS ZERO, THE CALCULATION
C 	IS TERMINATED AND A DIAGNOSTIC STATEMENT IS PRINTED.
C**********************************************************************

	  REAL*8  A(3,3), AINV(3,3), DET, ACOFAC(3,3), AADJ(3,3)

C	A(3,3)	        -- THE MATRIX WHOSE INVERSE IS DESIRED.
C	DET		-- THE COMPUTED DETERMINANT OF [A].
C	ACOFAC(3,3)	-- THE MATRIX OF COFACTORS OF A(I,J).
C			   THE SIGNED MINOR (-1)**(I+J)*M_IJ
C			   IS CALLED THE COFACTOR OF A(I,J).
C	AADJ(3,3)	-- THE ADJOINT OF [A]. IT IS THE MATRIX
C			   OBTAINED BY REPLACING EACH ELEMENT OF
C			   [A] BY ITS COFACTOR, AND THEN TAKING
C			   TRANSPOSE OF THE RESULTING MATRIX.
C	AINV(3,3)	-- RETURNED AS INVERSE OF [A].
C			   [AINV] = [AADJ]/DET.
C----------------------------------------------------------------------

	  CALL MDET(A,DET)
	  IF ( DET .EQ. 0.D0 ) THEN
  	    write(*,10)
!	    STOP
	  ENDIF
	  CALL MCOFAC(A,ACOFAC)
	  CALL MTRANS(ACOFAC,AADJ)
  	  DO 1 I = 1,3
	  DO 1 J = 1,3
	       AINV(I,J) = AADJ(I,J)/DET
1	  CONTINUE
10	FORMAT(5X,'--ERROR IN M3INV--- THE MATRIX IS SINGULAR',/,
     +         10X,'PROGRAM TERMINATED')

	  RETURN
	  END


C**********************************************************************
	  SUBROUTINE MCOFAC(A,ACOFAC)
 
C 	THIS SUBROUTINE CALCULATES THE COFACTOR OF A 3 BY 3 MATRIX [A],
C 	AND PLACES THE RESULT IN [ACOFAC]. 
C**********************************************************************

		REAL*8  A(3,3), ACOFAC(3,3)

		ACOFAC(1,1) = A(2,2)*A(3,3) - A(3,2)*A(2,3)
		ACOFAC(1,2) = -(A(2,1)*A(3,3) - A(3,1)*A(2,3))
		ACOFAC(1,3) = A(2,1)*A(3,2) - A(3,1)*A(2,2)
		ACOFAC(2,1) = -(A(1,2)*A(3,3) - A(3,2)*A(1,3))
		ACOFAC(2,2) = A(1,1)*A(3,3) - A(3,1)*A(1,3)
		ACOFAC(2,3) = -(A(1,1)*A(3,2) - A(3,1)*A(1,2))
		ACOFAC(3,1) = A(1,2)*A(2,3)  - A(2,2)*A(1,3)
		ACOFAC(3,2) = -(A(1,1)*A(2,3) - A(2,1)*A(1,3))
		ACOFAC(3,3) = A(1,1)*A(2,2) - A(2,1)*A(1,2)

	  RETURN
	  END    
!*******************************************************************************
!	      PROGRAM FOR SIMULATING THE COUTTE POISEUILLE CHANNEL FLOW
!			USING MIXING LENGTH TURBULENCE MODEL 
!		   	  d((MEU+MEUT)*(dU/dY))/dY = dP/dX
!*******************************************************************************

PROGRAM MAIN
IMPLICIT NONE

!******************************VARIABLE DECLARATION*****************************

!****************************GLOBALLY USED VARIABLES****************************

  INTEGER		::	NUM,SAMPLE
  REAL*8		::	MEU, PRGR, VW, H, RHO,ALP1

!******************************MESH AND VELOCITY *******************************

  INTEGER		::	N,MID,R,S
  REAL*8		::	ALPHA1, ALPHA2,Y,VQ,U1,U2,u1old,u2old,vqold 
  REAL*8,ALLOCATABLE	::	X2(:),U(:),UNP1(:)
  
!*******************************THOMAS ALGORITHM********************************
 
  REAL*8		::	DELI,DELJ
  REAL*8,ALLOCATABLE	::	A(:),B(:),C(:),D(:)
  
!******************************VELOCITY GRADIENT********************************
  
  REAL*8,ALLOCATABLE	::	GRADV_LOW(:), GRADV_UP(:)
  
!*************************STRESS AND TURBULENT VISCOSITY************************

  REAL*8		::	REY_STRESS
  REAL*8,ALLOCATABLE	::	TAU_UP(:),TAU_LOW(:),MEUT_LOW(:), MEUT_UP(:), &
				MEUTP1_LOW(:),MEUTP1_UP(:)		 

!**********************************ERROR****************************************

  REAL*8		::	E
  REAL*8,PARAMETER	::	ERR_BAND = 0.0001

!****************************OTHER VARIABLES************************************ 
  
  INTEGER		::	I,J,K,ITER,LP		
  REAL*8		::	MAX_VAL
  CHARACTER(LEN = 8) 	::	fmt 
  CHARACTER*200		::	POINT
  fmt = '(I2)'
  
!*********************OPENING FILE FOR READING INPUT VARIABLES******************
  OPEN(UNIT = 100, FILE = "INPUT.csv", ACTION = "READ")				
  READ(100,*)
  
!******************************PREPARING OUTPUT FILE****************************

  OPEN(UNIT = 101, FILE = "RESULT.dat", ACTION = "WRITE")	
  WRITE(101,*)"No	ITER	dP/dx		Vq		u1		u2"	

!********************************READING MESH DATA******************************

  OPEN(UNIT = 1,FILE = "MESH.dat",ACTION = "READ")				
  READ(1,*)N,MID,ALPHA1,ALPHA2
  
  R = MID-1
  S = MID+1
  
  ALLOCATE(X2(N),UNP1(N),A(N),B(N),C(N),D(N),U(N))
  ALLOCATE(GRADV_LOW(R),GRADV_UP(S:N),MEUT_LOW(R),MEUT_UP(S:N), 	 &
  		MEUTP1_LOW(R),MEUTP1_UP(S:N),TAU_UP(S:N),TAU_LOW(R))

  DO J = 1,N
   READ(1,*)X2(J)
  END DO 
  
  CLOSE(1)

!********************LOOP FOR CHECKING EACH CASES IN SEQUENCE*******************

  DO SAMPLE = 1,18
  
!*****************************READING INPUT VARIABLES***************************

    READ(100,*)NUM,MEU,RHO,H,Vw,ALP1
    WRITE (POINT,fmt) SAMPLE
    LP = LEN_TRIM(POINT)
    
!**************************CALCULATING PRESSURE GRADIENT ***********************

    PRGR=RHO*ALP1
    
   WRITE(*,5)MEU,PRGR,VW,H
5 FORMAT("MEU:-",F10.9,"  PRESSURE GRADIENT:-",F10.3,"  WALL VELOCIT:-",F7.3,"  HEIGHT OF THE CHANNEL:-",F7.3)

!****************************INITIAL CONDITION *********************************

  U = 0.0D0
  U(N) = 2.0D0*VW
  MEUT_LOW = 0.00001
  MEUT_UP = 0.00001
  u1old=0.0
  u2old=0.0
  vqold=0.0
  
!************************PREPARING FILE FOR WRITING ERROR***********************

  OPEN(UNIT = 5,FILE = "error"//point(1:lp)//".dat",ACTION = "WRITE")
  WRITE(5,*)"#      ITERATON	    ERROR"	
  
!******************LOOP FOR EACH CASE TO CONVERGE THE RESULTS*******************

  ITER = 0
  E = 1.0d0
  DO WHILE(E >= ERR_BAND)
    ITER = ITER +1 
    CALL VELGRAD (N, MID, H, U, X2, GRADV_LOW, GRADV_UP)
  
    CALL TURBULENCE_MODEL (N, MID, R, S, VW, X2, H, RHO, MEU,MEUT_UP, MEUT_LOW,&
    			  GRADV_UP, GRADV_LOW,MEUTP1_UP,MEUTP1_LOW,U1,U2,      &
    			  TAU_LOW,TAU_UP)
  
!**************************COEFFICIENT FOR LOWER HALF***************************

    DO I = 2,R
      DELI = X2(I)-X2(I-1)
      DELJ = X2(I+1)-X2(I)
      
      A(I) = (MEU+MEUTP1_LOW(I-1))/(H*DELI)
      B(I) = -(((MEU+MEUTP1_LOW(I))/(H*DELJ))+(MEU+MEUTP1_LOW(I-1))/(H*DELI))
      C(I) = (MEU+MEUTP1_LOW(I))/(H*DELJ)
      D(I) = 0.5*PRGR*(H*DELJ+H*DELI)
    END DO
    
!*************************COEFFICIENT FOR MID POINT*****************************

    I = MID
    DELI = X2(I)-X2(I-1)
    DELJ = X2(I+1)-X2(I)
  
    A(I) = (MEU+MEUTP1_LOW(I-1))/(H*DELI)
    B(I) = -(((MEU+MEUTP1_UP(I+1))/(H*DELJ))+(MEU+MEUTP1_LOW(I-1))/(H*DELI))
    C(I) = (MEU+MEUTP1_UP(I+1))/(H*DELJ)
    D(I) = 0.5*PRGR*(H*DELJ+H*DELI)
    
!******************************COEFICIENT FOR UPPER HALF************************

    DO I = N-1,S,-1
      DELI = X2(I)-X2(I-1)
      DELJ = X2(I+1)-X2(I)
    
      A(I) = (MEU+MEUTP1_UP(I))/(H*DELI)
      B(I) = -(((MEU+MEUTP1_UP(I+1))/(H*DELJ))+(MEU+MEUTP1_UP(I))/(H*DELI))
      C(I) = (MEU+MEUTP1_UP(I+1))/(H*DELJ)
      D(I) = 0.5*PRGR*(H*DELJ+H*DELI)
    END DO
 
!*******************************BOUNDARY CONDITIONS*****************************

     B(1) = 1.0D0
     C(1) = 1.0D0
     D(1) = 0.0D0
    
     A(N) = 0.5
     B(N) = 0.5
     D(N) = Vw

    CALL THOMAS(VW,N,A,B,C,D,UNP1)	
    
    CALL FLOW_RATE_VELOCITY(N, X2, H, UNP1, Vq)					
  
    CALL ERROR(N,U,UNP1,u1old,U1,u2old,U2,vqold,Vq,E)
    
    u1old=U1
    u2old=U2
    vqold=Vq
    
!*************************WRITING ERROR VALUES FOR EACH CASES*******************
    
    WRITE(5,*)ITER,"	",E
    
!*******************************UNDER RELAXATION FACTOR*************************
    
    U = 0.7*UNP1+(1.0D0-0.7)*U
    
!**************PASSING VALUES TO OTHER VARIABLES FOR NEXT ITERATION*************

    MEUT_LOW = MEUTP1_LOW
    MEUT_UP = MEUTP1_UP
    
  END DO
  CLOSE(5)
    
     
!*****************************PREPARING OUTPUT**********************************

  X2(1) = (X2(1)+X2(2))*0.5
  UNP1(1) = (UNP1(1)+UNP1(2))*0.5
  X2(N) = (X2(N)+X2(N-1))*0.5
  UNP1(N) = (UNP1(N)+UNP1(N-1))*0.5
  
  OPEN(UNIT = 1,FILE = "velocity"//point(1:lp)//".dat",ACTION = "WRITE")	
  OPEN(UNIT = 2,FILE = "stress"//point(1:lp)//".dat",ACTION = "WRITE")			
  
!*****************************SCALING HEIGHT AND VEOCITY************************

  IF(VW.EQ.0.0) THEN
    MAX_VAL=MAXVAL(UNP1)
  ELSE
    MAX_VAL=ABS(VW)
  END IF

  IF(SAMPLE.GT.15) THEN
    UNP1=UNP1*VW/VQ
  ELSE
    CONTINUE
  END IF
  
!***********************WRITING VELOCITY RESULTS OF EACH CASES******************
  
  WRITE(1,*)"#		X2			V"
  
  DO I = 1,N
   WRITE(1,*)X2(I),UNP1(I)/MAX_VAL
  END DO
 
!************************WRITING STRESS RESULTS FOR EACH CASES******************
 
  WRITE(2,*)"#	X2			Reynolds Stress			Tau_tot"
  
  Y = X2(1)
  REY_STRESS= MEUTP1_LOW(1)*GRADV_LOW(1)
  WRITE(2,*)Y,REY_STRESS,TAU_LOW(1)
  
  DO I=2,R
  Y = 0.5*(X2(I)+X2(I+1))
  REY_STRESS= MEUTP1_LOW(I)*GRADV_LOW(I)
  WRITE(2,*)Y,REY_STRESS,TAU_LOW(I)
  END DO
  
  DO I=S,N
  Y= 0.5*(X2(I)+X2(I-1))
  REY_STRESS = MEUTP1_UP(I)*GRADV_UP(I)
  WRITE(2,*)Y,REY_STRESS,TAU_UP(I)
  END DO
  
  Y=X2(N)
  REY_STRESS = MEUTP1_UP(N)*GRADV_UP(N)
  WRITE(2,*)Y,REY_STRESS,TAU_UP(N)
  
!***********************WRITING OTHER RESULTS FOR EACH CASES********************

  CLOSE(1)
  CLOSE(2)
  WRITE(101,15)SAMPLE,iter,PRGR,VQ,U1,U2
15 FORMAT(I2,"	",I2,"	",F10.5,"	",F7.3,"		",F7.4,"		",F7.4)
  END DO
  CLOSE(100)
  CLOSE(101)
CONTAINS 



!*******************************************************************************
!		SUBROUTINE FOR CALCULATING THE VELOCITY GRADIENT
!*******************************************************************************

SUBROUTINE VELGRAD (N, MID, H, U, X2, GRADV_LOW, GRADV_UP)
IMPLICIT NONE

!******************************VARIABLE DECLARATION*****************************

  INTEGER, INTENT(IN)	::	N,MID
  REAL*8, INTENT(IN)	::	U(N),X2(N),H
  REAL*8, INTENT(OUT)	::	GRADV_LOW(MID-1),GRADV_UP((MID+1):N)
  INTEGER		::	I
  REAL*8		::	DELX2
  
!*************************VELOCITY GRADIENT FOR LOWER HALF**********************

  DO I=1,MID-1
    DELX2=X2(I+1)-X2(I)
    GRADV_LOW(I)=(U(I+1)-U(I)) / (H*DELX2)
  END DO
  
!***********************VELOCITY GRADIENT FOR UPPER HALF************************
  
  DO I=N,MID+1,-1
    DELX2=X2(I)-X2(I-1)
    GRADV_UP(I)=(U(I)-U(I-1))/(H*DELX2)
  END DO

END SUBROUTINE VELGRAD




!*******************************************************************************
!		SUBROUTINE FOR CALCULATING THE EDDY VISCOSITY
!*******************************************************************************

SUBROUTINE TURBULENCE_MODEL (N,MID,R,S,VW,X2,H,RHO,MEU,MEUT_UP, MEUT_LOW,      &
				GRADV_UP, GRADV_LOW,MEUTP1_UP, MEUTP1_LOW,     &
				U1,U2,TAU_LOW,TAU_UP)
IMPLICIT NONE

!********************************VARIABLE DECLARATION***************************


  INTEGER, INTENT(IN)	::	N,MID,R,S
  REAL*8, INTENT(IN)	::	VW,X2(N),H,RHO,MEU,MEUT_LOW(R), MEUT_UP(S:N),  &
  				GRADV_LOW(R),GRADV_UP(S:N)
   REAL*8,INTENT(OUT)	::	U1,U2,MEUTP1_LOW(R),MEUTP1_UP(S:N), TAU_LOW(R),&
  				TAU_UP(S:N)
  INTEGER		::	I
  REAL*8 		::	Y,L0_LOW(R), L0_UP(S:N),LM_LOW(R), LM_UP(S:N)
  
!***************************CALCULATION FOR LOWER WALL**************************

  DO I = 1,R
    Y = (X2(I+1)+X2(I))*0.5
    TAU_LOW(I) = ((MEU+MEUT_LOW(I))*GRADV_LOW(I))
    L0_LOW(I) = 0.5*H*(0.21-(0.43*(1.0D0-2.0D0*Y)**4)+(0.22*(1.0D0-2.0D0*Y)**6))
    LM_LOW(I) = L0_LOW(I)*(1.0D0-DEXP(-H*Y*DSQRT(DABS(TAU_LOW(I)/RHO))/ &
    							(30.0D0*MEU/RHO)))
    MEUTP1_LOW(I)=(LM_LOW(I)**2)*DABS(GRADV_LOW(I))*RHO
  END DO
 
 U1 = DSQRT(DABS(TAU_LOW(1)/RHO))
 
!**************************CALCULATION FOR UPPER WALL***************************
  DO I = N,S,-1
    Y = 1.0D0-(X2(I)+X2(I-1))*0.5
    TAU_UP(I) = ((MEU+MEUT_UP(I))*GRADV_UP(I))
    L0_UP(I) = 0.5*H*(0.21-(0.43*(1.0D0-2.0D0*Y)**4)+(0.22*(1.0D0-2.0D0*Y)**6))
    LM_UP(I) = L0_UP(I)*(1.0D0-DEXP(-H*Y*DSQRT(DABS(TAU_UP(I)/RHO))/	  &
    							(30.0D0*MEU/RHO)))
    MEUTP1_UP(I) = (LM_UP(I)**2)*DABS(GRADV_UP(I))*RHO
  END DO
 
  U2 = DSQRT(DABS(TAU_UP(N)/RHO))

 
END SUBROUTINE TURBULENCE_MODEL




!*******************************************************************************
!			SUBROUTINE FOR INVERTING THE MATRIX
!*******************************************************************************

SUBROUTINE THOMAS (VW, N, A, B, C, D, X)					
IMPLICIT NONE

!*****************************VARIABLE DECLARATION******************************

  INTEGER, INTENT(IN) 	:: 	N						
  REAL*8 		:: 	VW,A(N),B(N),C(N),D(N),X(N)
  REAL*8		::	V,Q
  INTEGER		::	I,J
 
!*****************************FORWARD ELIMINATION*******************************

  DO I = 2,N			
    Q = A(I)/B(I-1)
    B(I) = B(I)-C(I-1)*Q
    D(I) = D(I)-D(I-1)*Q
  END DO

!*****************************BACKWARD SUBSTITUTION*****************************
   
  X(N) = D(N)/B(N)
  DO I = N-1,1,-1								
    V = (D(I)-C(I)*X(I+1))/B(I)
    X(I) = V
  END DO
    			
  RETURN
  END SUBROUTINE THOMAS
  



!*******************************************************************************
!		SUBROUTINE FOR CALCULATING THE FLOW RATE VELOCITY
!*******************************************************************************

SUBROUTINE FLOW_RATE_VELOCITY (N, Y, H, U, Vq)
IMPLICIT NONE

!*****************************VARIABLE DECLARATION******************************

INTEGER, INTENT(IN) 	::	N
REAL*8, INTENT(IN)	::	U(N), H, Y(N)
REAL*8,INTENT(OUT)	::	Vq
   
!**********INTEGRATING VELOCITY PROFILE TO OBTAIN FLOWRATE VELOCITY*************

  Vq=0.
  DO I = 1,N-1
   Vq = Vq+((Y(I+1)-Y(I))*(U(I)+U(I+1))/2.0D0)
  END DO

END SUBROUTINE FLOW_RATE_VELOCITY


		
		
		
!*******************************************************************************
!			SUBROUTINE FOR CALCULATING THE ERROR
!*******************************************************************************	
	
SUBROUTINE ERROR (N, U, UNP1,u1old,U1,u2old,U2,vqold,Vq, E)
IMPLICIT NONE

!*****************************VARIABLE DECLARATION******************************

  INTEGER, INTENT(IN)	::	N
  REAL*8, INTENT(IN)	::	U(N),UNP1(N),u1old,u1,u2old,u2,vqold,vq
  REAL*8, INTENT(OUT)	::	E
  REAL*8		::	X(N),E1,Er_chck(4)
     
!****************************CALCULATING SECOND NORM****************************
  
  E1=0.0D0
  DO I=1,N
     E1=E1+(ABS(U(I))-ABS(UNP1(I)))**2
  END DO
  Er_chck(1)=DSQRT(E1)/N
  
!*************************ERROR ON OTHER PARAMATERS*****************************
  Er_chck(2)=DABS(u1old-U1)
  Er_chck(3)=DABS(u2old-U2)
  Er_chck(4)=DABS(vqold-Vq)
  E=MAXVAL(Er_chck)
  
  
    
END SUBROUTINE ERROR

END PROGRAM MAIN

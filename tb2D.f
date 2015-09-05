      PROGRAM TB2D

C This program computes the pressure distribution and several derived
C quantities for a finite tilting thrust bearing of rectangular shape. 
C Calculations can be performed for different L/B ratios.

C X(I)     = position along x-axis (ie. x coordinate)
C Y(J)     = position along y-axis (ie. y coordinate)
C NX       = number of mesh points along x-axis
C NY       = number of mesh points along y-axis 
C P(I,J)   = Pressure at node i,j at latest iteration level
C PPI(I,J) = Pressure at node i,j at the previous iteration
C H(I)     = film thickness at node i,j
C LB       = L/B ratio
C K        = H1/H0-1 
C H1       = leading edge (maximum) film thickness
C H0       = trailing edge (minimum) film thickness
C NDLOAD   = Non-dimensional load capacity
C QX       = Flow-rate per unit length in the x-direction
C QY       = Flow-rate per unit breadth in the y-directon
C QE       = End-leakage (evaluated at the trailing edge)
C QS       = Side-leakage (evaluated at y* = 0.0)
C DPDX     = Pressure gradient in the x-direction
C DPDY     = Pressure gradient in the y-direction
C TORX0    = Shear stress on the lower surface in the x-direction
C TORXH    = Shear stress on the upper surface in the x-direction
C FF0      = Friction force on the lower surface
C FFH      = Friction force on the upper surface
C XBAR     = Location of the pivot point 
C FOR      = Successive over-relaxation factor

C     Declaration of variables
      PARAMETER (N = 501, NTOTAL=N*N, PI=3.141592654, FOR=1.5)
      DOUBLE PRECISION K, LB, DX, DY, DX2, DY2, BL2
      DOUBLE PRECISION X(N), Y(N), H(N), PPI(N,N), P(N,N)
      DOUBLE PRECISION A(N), B(N), C(N), D(N)
      DOUBLE PRECISION NDLOAD, SGLINT, DBLINT
      DOUBLE PRECISION QX(N), QY(N), QE, QS 
      DOUBLE PRECISION DPDX(N,N), DPDY(N,N), TORX0(N,N), TORXH(N,N) 
      DOUBLE PRECISION XP(N,N), XBAR, FF0, FFH
      INTEGER NX, NXM1, NXM2, NY, NYM1, NYM2, ITER
      LOGICAL FLAG

C     Initialise variables
      DATA P /NTOTAL*0.0D0/, PPI /NTOTAL*0.0D0/
      DATA FLAG /.FALSE./, ITER /1/

C     Open results file
      OPEN (UNIT=1, FILE='tb2d_output.dat', STATUS='UNKNOWN')
      OPEN (UNIT=2, FILE='tb2d_plot.dat',   STATUS='UNKNOWN')
      OPEN (UNIT=3, FILE='tb2d_script.m',   STATUS='UNKNOWN')

C     Prompt user for input values
      WRITE(*,*) 'Enter value of L/B ratio'
      READ (*,*) LB
      WRITE(*,*) 'Enter value of K=(h1/h0)-1'
      READ (*,*) K
      WRITE(*,*) 'Enter value of NX'
      READ (*,*) NX
      WRITE(*,*) 'Enter value of NY'
      READ (*,*) NY

C     Dependent variables
      NXM1 = NX-1
      NXM2 = NX-2
      NYM1 = NY-1
      NYM2 = NY-2
      DX  = 1.0/NXM1
      DY  = 1.0/NYM1
      DX2 = DX**2
      DY2 = DY**2
      BL2 = (1/LB)**2

C     Initialise discrete grid points in the
C     (a) x-direction
      X(1)=0.0
      DO 5,I=2,NX
        X(I)=X(I-1)+DX
    5 CONTINUE

C     (b) y-direction
      Y(1)=0.0
      DO 10,J=2,NY
        Y(J)=Y(J-1)+DY
   10 CONTINUE

C     Initialise film thickness
      DO 15,I=1,NX
        H(I) = 1.0+K*(1-X(I))
   15 CONTINUE

       
C ---------------  ITERATION PROCESS  ---------------


      DO WHILE(FLAG .EQV. (.FALSE.))
           
C       Provide user feedback during iteration process
        WRITE(*,'(10X,A12,I5)') 'Iteration = ', ITER
        ITER = ITER + 1

C       Sweep all rows j=2,3,...,J-1
        DO 50,J=2,NYM1
          DO 35,I=2, NXM1
            A(I-1) = (1.0/DX2)+(1.5*K)/(DX*H(I))
            B(I-1) = -2.0*((1.0/DX2)+(BL2/DY2))
            C(I-1) = (1.0/DX2)-(1.5*K)/(DX*H(I))
            D(I-1) = -K/(H(I)**3)-(BL2/DY2)*(P(I,J-1)+P(I,J+1))
   35     CONTINUE
          
C         Use Thomas algorithm to solve tridiagonal system of algebraic equations          
          CALL THOMAS(NXM2,A,B,C,D,N)          
          
C         Capture solution from Thomas and apply SOR factor, For
          DO 40,I=2,NXM1            
            P(I,J)   = PPI(I,J)+FOR*(D(I-1)-PPI(I,J))  
            PPI(I,J) = P(I,J)
   40     CONTINUE             

   50   CONTINUE       
    
C       Sweep all columns i=2,3,...,I-1
        DO 65,I=2,NXM1          
          DO 55,J=2, NYM1            
            A(J-1) = (BL2/DY2)            
            B(J-1) = -2.0*((1.0/DX2)+(BL2/DY2))            
            C(J-1) = (BL2/DY2)            
            D(J-1) = -K/(H(I)**3)
     +               -((1.0/DX2)+(1.5*K)/(DX*H(I)))*P(I-1,J)
     +               -((1.0/DX2)-(1.5*K)/(DX*H(I)))*P(I+1,J)   
   55     CONTINUE                   
          
C         Use Thomas algorithm to solve tridiagonal system of algebraic equations
          CALL THOMAS(NYM2,A,B,C,D,N)
          
C         Capture solution from Thomas and apply SOR factor, For
          DO 60,J=2,NYM1
            P(I,J) = PPI(I,J)+FOR*(D(J-1)-PPI(I,J))
   60     CONTINUE

   65   CONTINUE

C       Convergence test
        CALL CVERGE(NX,NY,P,PPI,FLAG,N)
        
      END DO

C ---------------  END ITERATION PROCESS  --------------- 


C     Calculate the non-dimensional load
      NDLOAD = DBLINT(P,DX,DY,NXM1,NYM1,N)

C     Calculate location of the pivot point, XBAR
      DO 75,J=1,NY
        DO 70,I=1,NX
          XP(I,J)=X(I)*P(I,J)
   70   CONTINUE
   75 CONTINUE
      XBAR = DBLINT(XP,DX,DY,NXM1,NYM1,N)/NDLOAD

C     Calculate the pressure gradients
C     (a) in the x-direction
      DO 85,J=1,NY
        DPDX(1,J)  = (-3.0*P(1,J)+4.0*P(2,J)-P(3,J))/(2.0*DX)
        DPDX(NX,J) = (3.0*P(NX,J)-4.0*P(NXM1,J)+P(NXM2,J))/(2.0*DX)
        DO 80,I=2,NXM1
          DPDX(I,J) = (P(I+1,J)-P(I-1,J))/(2.0*DX)
   80   CONTINUE
   85 CONTINUE
      
C     (b) in the y-direction
      DO 90,I=1,NX
        DPDY(I,1) = (-3.0*P(I,1)+4.0*P(I,2)-P(I,3))/(2.0*DY)
   90 CONTINUE
     
C     Calculate the friction force by integration of the shear 
C     stress over the bearing area.  The friction force on the
C     lower surface is FFH, whereas that on the lower surface is FF0
      DO 100,J=1,NY
        DO 95,I=1,NX
          TORX0(I,J) = -3.0*H(I)*DPDX(I,J)-1.0/H(I)
          TORXH(I,J) =  3.0*H(I)*DPDX(I,J)-1.0/H(I)
   95   CONTINUE
  100 CONTINUE     
      FF0 = DBLINT(TORX0,DX,DY,NXM1,NYM1,N)
      FFH = DBLINT(TORXH,DX,DY,NXM1,NYM1,N)

C     Calculate the side leakage, Qs
      DO 105,I=1,NX
        QY(I) = -H(I)**3/LB*DPDY(I,1)
  105 CONTINUE
      QS = -SGLINT(QY,DX,NXM1,N)

C     Calculate the end leakage, Qe
      DO 110,J=1,NY
        QX(J) = 0.5*LB*(H(NX)-(H(NX)**3)*DPDX(NX,J))
  110 CONTINUE
      QE = SGLINT(QX,DY,NYM1,N)     


C --------------- PRINT RESULTS ---------------


C     Write results to end of existing results file      
  115 READ(1,*,END=120)
      GO TO 115
  120 CONTINUE
      WRITE(1,125) LB,K,6*NDLOAD,XBAR,FF0,FFH,QS,QE
  125 FORMAT(2F10.2,6F10.5)

C     Write the pressure profile and pressure gradient  
C     along the midplane to results file
      J=(NY+1)/2
      DO 130,I=1,NX 
        WRITE(2,'(3F10.6)') X(I), P(I,J), DPDX(I,J) 
  130 CONTINUE

C     Create Matlab script file to plot the non-dimensional film
C     thickness and the pressure distribution over the bearing
      WRITE(3,140) 'x=['
      DO 135,J=1,NY
        WRITE(3,145) (X(I),I=1,NX)
  135 CONTINUE
  140 FORMAT(A3)
  145 FORMAT(101F10.5)
      WRITE(3,*) '];'
      WRITE(3,*)

      WRITE(3,140) 'y=['
      DO 150,J=1,NY
        WRITE(3,145) (Y(J),I=1,NX)
  150 CONTINUE
      WRITE(3,*) '];'
      WRITE(3,*)

      WRITE(3,140) 'h=['
      DO 155,J=1,NY
        WRITE(3,145) (H(I),I=1,NX)
  155 CONTINUE
      WRITE(3,*) '];'
      WRITE(3,*)

      WRITE(3,140) 'p=['
      DO 160,J=1,NY
        WRITE(3,145) (P(I,J),I=1,NX)
  160 CONTINUE
      WRITE(3,*) '];'
      WRITE(3,*)

C     3D plot of the film thickness
      WRITE(3,*) 'figure(1);'
      WRITE(3,*) 'hold on;'
      WRITE(3,*) 'surf(x,y,h);'
      WRITE(3,*) 'view(-30,30);'
      WRITE(3,*) 'grid on;'
      WRITE(3,*) 'box on;'
      WRITE(3,*) 'axis([-Inf Inf -Inf Inf 0 Inf]);'
      WRITE(3,*) 'xlabel(''x^*'');'
      WRITE(3,*) 'ylabel(''y^*'');'
      WRITE(3,*) 'zlabel(''h^*'');'

C     3D plot of the pressure distribution
      WRITE(3,*) 'figure(2);'
      WRITE(3,*) 'hold on;'
      WRITE(3,*) 'surf(x,y,p);'
      WRITE(3,*) 'view(30,30);'
      WRITE(3,*) 'grid on;'
      WRITE(3,*) 'box on;'
      WRITE(3,*) 'xlabel(''x^*'');'
      WRITE(3,*) 'ylabel(''y^*'');'
      WRITE(3,*) 'zlabel(''p^*'');'

C     Contour plot of the pressure distribution
      WRITE(3,*) 'figure(3);'
      WRITE(3,*) 'hold on;'
      WRITE(3,*) 'contourf(x,y,p,10,''k'');'
      WRITE(3,*) 'box on;'
      WRITE(3,*) 'xlabel(''x^*'');'
      WRITE(3,*) 'ylabel(''y^*'');'

C     Close results files
      END FILE (UNIT=1)
      END FILE (UNIT=2)
      END FILE (UNIT=3)
      CLOSE (UNIT=1)
      CLOSE (UNIT=2)
      CLOSE (UNIT=3)

      STOP
      END


C ======================================================================      
      
      SUBROUTINE THOMAS(NN,A,B,C,D,N)      

C     This subroutine uses the Thomas algorithm for solving tri-diagonal
C     systems of algebraic equations with Dirichlet boundary conditions.  
C     It is used in iterative line-by-line sweeps to determine the pressure
C     distribution over the bearing surfaces. 
      
      DOUBLE PRECISION A(N),B(N),C(N),D(N)      
      INTEGER N, NN

      B(1)=1.0/B(1)      

      DO 5,I=2,NN        
        B(I)=1.0/( B(I)-A(I)*B(I-1)*C(I-1) )        
        D(I)=D(I)-A(I)*B(I-1)*D(I-1)    
    5 CONTINUE      

      D(NN)=D(NN)*B(NN)      

      DO 10,I=NN-1,1,-1        
        D(I)=( D(I)-D(I+1)*C(I) )*B(I)   
   10 CONTINUE      
  
      RETURN      
      END 


C ======================================================================

      SUBROUTINE CVERGE(NX,NY,P,PPI,FLAG,N)

C     This subroutine tests for convergence of the solution. If the 
C     residual is less than or equal to RESLIM then the solution is 
C     converged.  If the residual is greater than RESMAX, the solution 
C     is diverging and the program is stopped.  Otherwise, an over -
C     relaxation factor is applied to speed up the iterative process 
C     (This is performed in subroutines TB and CJB)

      PARAMETER (RESLIM=1.0D-6, RESMAX=1.0D2)
      DOUBLE PRECISION P(N,N), PPI(N,N), RESID
      LOGICAL FLAG

      RESID=0.0D0
      DO 10,J=1,NY
        DO 5,I=1,NX
          RESID = MAX(RESID,ABS(P(I,J)-PPI(I,J)))
    5 CONTINUE
   10 CONTINUE
      IF (RESID .LE. RESLIM) THEN
        FLAG = .TRUE.
      ELSE IF (RESID .GE. RESMAX) THEN
        STOP 'Solution not converging'
      ELSE
        DO 20,J=2,NY-1
          DO 15,I=2,NX-1
            PPI(I,J) = P(I,J)          
   15     CONTINUE
   20   CONTINUE
      END IF

      RETURN
      END


C ======================================================================            

      DOUBLE PRECISION FUNCTION SGLINT(FX,DX,NXM1,N)

C     This subroutine performs a single integration using the
C     trapezoidal rule of integration.

      DOUBLE PRECISION SUM, DX, FX(N)
      INTEGER I, NXM1, N

      SUM = 0.0
      DO 5,I=1,NXM1
        SUM = SUM + FX(I) + FX(I+1)
    5 CONTINUE

      SGLINT = 0.5*DX*SUM
      END


C ======================================================================            

      DOUBLE PRECISION FUNCTION DBLINT(FXY,DX,DY,NXM1,NYM1,N)            
 
C     This subroutine performs a double integration using the
C     trapezoidal rule of integration.

      DOUBLE PRECISION SUM, DX, DY, FXY(N,N)      
      INTEGER NXM1, NYM1, N      

      SUM=0.0      
      DO 10,J=1,NYM1        
        DO 5,I=1,NXM1          
          SUM = SUM + FXY(I,J) + FXY(I,J+1) +     
     +          FXY(I+1,J) + FXY(I+1,J+1)    
    5   CONTINUE   
   10 CONTINUE      

      DBLINT = 0.25*DX*DY*SUM      
      END

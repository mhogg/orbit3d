      PROGRAM ORBIT3D

C     This program was developed with the aim of determining the
C     stability of the VentrAssist implantable rotary blood pump.
C     Orbit3D determines the pressure distribution over the
C     thrust bearing and conical journal bearing surfaces of the
C     impeller, calculates the forces and moments on the impeller
C     resulting from this pressure and calculates the position 
C     (both translational and rotational) of the impeller within
C     pump housing. From observation of this 'orbit', the stability
C     of the impeller can be assessed.

C     Description of variables:
C
C I          = Index of grid points in the theta direction
C J          = Index of grid points in the radial (r and rho) directions
C M          = Index of impeller blades
C N          = Index of grid points in time
C NTH        = Number of grid points in the theta direction
C NR         = Number of grid points in the radial (r) direction
C NRO        = Number of grid points in the radial (rho) direction
C NP         = Number of impeller blades
C NT         = Number of temporal grid points
C DTH        = Uniform grid spacing in the theta direction (non-dim.)
C DR         = Uniform grid spacing in the radial (r) direction (non-dim.)
C DRO        = Uniform grid spacing in the radial (rho) direction (non-dim.)
C DT         = Uniform grid spacing in time (non-dim.)
C TH(,)      = Theta location of grid points on impeller m (non-dim.)
C R()        = Radial (r) location of grid points (non-dim.)
C RO()       = Radial (rho) location of grid points (non-dim.)
C T()        = Time of temporal grid points (non-dim.)
C NDT        = Time at the nth time step (non-dim.)
C REV        = Number of revolution of the impeller
C ALPHD      = Half-cone angle (degrees)
C ALP        = Half-cone angle (radians)
C TH0D       = The angular span of each impeller blade (deg.)
C TH0        = The angular span of each impeller blade (rad.)
C TH1        = Spacing of impeller blades (rad.) = 2*PI/NP
C R1         = The inner radius of the thrust bearing (m)
C R2         = The outer radius of the thrust bearing (m)
C LR         = R2 - R1 (m)
C RO1        = The inner radius of the conical journal bearing in the rho
C              direction (m)
C RO2        = The outer radius of the conical journal bearing in the rho
C              direction (m)
C LRO        = RO2 - RO1 (m)
C PHI1       = Angle used to specify the gravitational force (rad.)
C PHI2       = Angle used to specify the gravitational force (rad.)
C LOADP      = Load parameter (non-dim.)
C XH,YH,ZH   = Amplitude of the displacement of the pump housing (non-dim.)
C V1         = Frequency of "shaking" of the pump housing (non-dim.)
C V2         = Frequency of unbalance force (non-dim.)
C H0TB       = Minimum (outlet) film thickness of thrust bearing (m) 
C H1TB       = Maximum (inlet) film thickness of thrust bearing (m)
C KTB        = H1TB/H0TB-1 (non-dim.)
C PLTB       = Fraction of thrust bearing that is untapered land
C HVC        = Vertical clearance of conical journal bearing = H0TB (m)
C H0CJB      = Minimum (outlet) film thickness of conical journal bearing (m) 
C H1CJB      = Maximum (inlet) film thickness of conical journal bearing (m)
C KCJB       = H1CJB/H0CJB-1 (non-dim.)
C PLCJB      = Fraction of conical journal bearing that is untapered land
C IT         = Transverse moment of inertia (non-dim.)
C JP         = Polar moment of inertia (non-dim.)
C M1         = Mass of the impeller (non-dim.)
C M2         = Magnitude of the unbalance mass (non-dim.) 
C PTB(,,)    = Pressure value at node i,j on the impeller blade m of the 
C              thrust bearing at time step N (non-dim.) 
C PPITB(,,)  = Pressure value at node i,j on the impeller blade m of the 
C              thrust bearing at time step N-1 (non-dim.)
C PCJB(,,)   = Pressure value at node i,j on the impeller blade m of the 
C              conical journal bearing at time step N (non-dim.)
C PPICJB(,,) = Pressure value at node i,j on the impeller blade m of the 
C              conical journal bearing at time step N-1 (non-dim.)
C FXTB       = Force in the X-direction on the impeller from the thrust 
C              bearing (non-dim.)
C MXTB       = Moment about the X-axis on the impeller from the thrust
C              bearing (non-dim.)
C MYTB       = Moment about the Y-axis on the impeller from the thrust 
C              bearing (non-dim.)
C FXCJB      = Force in the X-direction on the impeller from the conical 
C              journal bearing (non-dim.)
C FYCJB      = Force in the Y-direction on the impeller from the conical
C              journal bearing (non-dim.)
C FZCJB      = Force in the Z-direction on the impeller from the conical
C              journal bearing (non-dim.)
C MXCJB      = Moment about the X-axis on the impeller from the conical
C              journal bearing (non-dim.)
C MYCJB      = Moment about the Y-axis on the impeller from the conical
C              journal bearing (non-dim.)
C FX(,)      = Total force in the X-direction on the impeller (non-dim.)
C FY(,)      = Total force in the Y-direction on the impeller (non-dim.)
C FZ(,)      = Total force in the Z-direction on the impeller (non-dim.)
C MX(,)      = Total moment about the X-axis on the impeller (non-dim.)
C MY(,)      = Total moment about the Y-axis on the impeller (non-dim.)
C K1()       = Variable used in the Runge-Kutta method
C K2()       = Variable used in the Runge-Kutta method
C K3()       = Variable used in the Runge-Kutta method
C K4()       = Variable used in the Runge-Kutta method
C Y1()       = Eccentricity ratio in the X-direction
C Y2()       = Instantaneous velocity in the X-direction (non-dim.) 
C Y3()       = Eccentricity ratio in the Y-direction
C Y4()       = Instantaneous velocity in the Y-direction (non-dim.)
C Y5()       = Eccentricity ratio in the Z-direction
C Y6()       = Instantaneous velocity in the Z-direction (non-dim.)
C Y7()       = Rotational position of the impeller about the X-axis (non-dim.)
C Y8()       = Rotational velocity of the impeller about the X-axis (non-dim.)
C Y9()       = Rotational position of the impeller about the Y-axis (non-dim.)
C Y10()      = Rotational velocity of the impeller about the Y-axis (non-dim.)
C FLAG1      = Set to true if touchdown of the thrust bearing occurs, otherwise false
C FLAG2      = Set to true if touchdown of the conical journal bearing occurs, 
C              otherwise false

C     Variable declaration
      PARAMETER (N1=101, N2=4, N3=5761)
      PARAMETER (NTH=21,NR=21,NRO=21,NP=4,NT=721,REV=1.0D0)
      PARAMETER (ALPD=45.0, TH0D=45.0, R1=0.0125, R2=0.025) 
      PARAMETER (PI=3.141592654, RTOD=180.0/PI, DTOR=PI/180.0)
      PARAMETER (PHI1=PI, PHI2=0.0, LOADP=152.03)
      PARAMETER (XH=10.0, YH=0.0, ZH=0.0, V1=1.0, V2=0.0)
      REAL H0TB, H1TB, KTB, PLTB, HVC, H0CJB, H1CJB, KCJB, PLCJB
      REAL LR, RO1, RO2, LRO, IT, JP, M1, M2
      REAL FZTB, MXTB, MYTB
      REAL FXCJB, FYCJB, FZCJB, MXCJB, MYCJB
      REAL FX(N3,4), FY(N3,4), FZ(N3,4), MX(N3,4), MY(N3,4)
      DOUBLE PRECISION K1(10), K2(10), K3(10), K4(10)
      DOUBLE PRECISION Y1(N3), Y2(N3), Y3(N3), Y4(N3), Y5(N3)
      DOUBLE PRECISION Y6(N3), Y7(N3), Y8(N3), Y9(N3), Y10(N3)
      DOUBLE PRECISION TH0, ALP, TH1, DTH, DR, DRO, DT 
      DOUBLE PRECISION TH(N1,N2), R(N1), RO(N1), T(N3)
      DOUBLE PRECISION NDT, SA, CA
      DOUBLE PRECISION PTB(N1,N1,N2), PPITB(N1,N1,N2) 
      DOUBLE PRECISION PCJB(N1,N1,N2), PPICJB(N1,N1,N2)
      INTEGER I, J, M, N, COUNT
      LOGICAL FLAG1, FLAG2
      DATA FLAG1/.FALSE./, FLAG2/.FALSE./, COUNT/1/

C     Open results files
      OPEN (UNIT=1, FILE='plottraj.m', STATUS='UNKNOWN')          
      OPEN (UNIT=2, FILE='plotcjb.m',  STATUS='UNKNOWN')          
      OPEN (UNIT=3, FILE='plottb.m',   STATUS='UNKNOWN')          

C     Initialised dependent variables for the
C     (a) thrust bearing
      H0TB = 100E-6
      H1TB = H0TB+50E-6
      KTB  = H1TB/H0TB-1.0
      PLTB = 0.2

C     (b) conical journal bearing
      ALP   = ALPD*DTOR
      SA    = DSIN(ALP)
      CA    = DCOS(ALP)
      RO1   = R1/SA
      RO2   = R2/SA
      ROG   = 1.05*RO2
      ROG   = ROG/(RO2*SA)
      HVC   = H0TB
      H0CJB = HVC*SA
      H1CJB = H0CJB+50E-6
      KCJB  = H1CJB/H0CJB-1.0 
      PLCJB = PLTB

C     Initialise the coordinates of the discrete grid points 
C     in the circumferential direction for the first impeller blade
      TH0=TH0D*DTOR
      DTH=TH0/(NTH-1)
      TH(1,1)=0.0
      DO 5,I=2,NTH
        TH(I,1)=TH(I-1,1)+DTH
    5 CONTINUE
 
C     Initialise the coordinates of the discrete grid points 
C     in the circumferential direction for the remaining blades
      TH1=2*PI/NP
      DO 15,M=2,NP
        DO 10,I=1,NTH
          TH(I,M)=TH(I,1)+(M-1)*TH1
   10   CONTINUE
   15 CONTINUE

C     Initialise the coordinates of the grid points in the radial
C     direction for the thrust bearing  
      LR=R2-R1
      DR=(LR/R2)/(NR-1)
      R(1)=R1/R2
      DO 20,J=2,NR
        R(J)=R(J-1)+DR
   20 CONTINUE

C     Initialise the coordinates of the grid points in the radial
C     direction for the conical journal bearing
      LRO=RO2-RO1
      DRO=LRO/R2/(NRO-1)
      RO(1)=RO1/R2
      DO 25,J=2,NRO
        RO(J)=RO(J-1)+DRO
   25 CONTINUE

C     Initialise the numerical grid in time. The non-dimensional
C     time step should be in the order of 1 degree, but is typically
C     0.5 or even 0.25 degrees depending on the operating conditions. 
      DT=2.0*PI*REV/(NT-1)
      T(1)=0.0D0
      DO 30,N=2,NT
        T(N)=T(N-1)+DT
   30 CONTINUE

C     Initialise pressure to ambient (ie. zero) at t=0 for the
C     (a) thrust bearing
      DO 45,M=1,NP
        DO 40,J=1,NR
          DO 35,I=1,NTH
            PTB(I,J,M)   = 0.0D0
            PPITB(I,J,M) = 0.0D0
   35     CONTINUE
   40   CONTINUE
   45 CONTINUE

C     (b) conical journal bearing
      DO 60,M=1,NP
        DO 55,J=1,NRO
          DO 50,I=1,NTH
            PCJB(I,J,M)   = 0.0D0
            PPICJB(I,J,M) = 0.0D0
   50     CONTINUE
   55   CONTINUE
   60 CONTINUE

C     Set properties of the impeller
      IT = 0.214097
      JP = 1.78578
      M1 = 0.69867
      M2 = 0.0

C     Set initial conditions
      Y1(1)  = 0.0
      Y2(1)  = 0.0
      Y3(1)  = 0.0
      Y4(1)  = 0.0
      Y5(1)  = 0.0
      Y6(1)  = 0.0
      Y7(1)  = 0.0
      Y8(1)  = 0.0
      Y9(1)  = 0.0
      Y10(1) = 0.0

C     The following DO WHILE loop uses the classical form of the fourth-order
C     Runge-Kutta method to solve the five second-order non-linear differential
C     equations that govern the translational and rotational movement of the
C     impeller.  This method requires four "estimates" before giving the 
C     actual position for each time step.  Thus four steps are involved.

C     DO WHILE loop control:
C     If FLAG1 is becomes true during the course of the run then touchdown 
C     has occurred in the thrust bearing.  Furthermore, if FLAG2 becomes true 
C     then touchdown has occurred in the conical journal bearing. If either 
C     of these occurs then the loop is exited and the program stopped.
      N=1

      DO WHILE ((FLAG1 .EQV. (.FALSE.)) .AND. (FLAG2 .EQV. (.FALSE.))
     +                                  .AND. (N .LE. NT))


C       *** FIRST STEP ***

        NDT = T(N)       

C       Call Subroutine TB (ie. Thrust Bearing) to calculate the pressure
C       distribution over the thrust bearing.  This subroutine also integrates
C       this pressure to determine the forces and moments on the impeller.
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT,     
     +          N1,N2,Y5(N),Y6(N),Y7(N),Y8(N),Y9(N),Y10(N),PTB,PPITB,      
     +          FZTB,MXTB,MYTB,FLAG1)                                    
                                                                         
C       Call Subroutine CJB (ie. Conical Journal Bearing) to calculate the
C       pressure distribution over the conical journal bearing.  This subroutine
C       also integrates this pressure to determine the forces and moments on the
C       impeller.                                                       
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG, 
     +       RO,DRO,NRO,NP,NDT,N1,N2,Y1(N),Y2(N),Y3(N),Y4(N),Y5(N),
     +       Y6(N),Y7(N),Y8(N),Y9(N),Y10(N),PCJB,PPICJB,FXCJB,FYCJB,
     +       FZCJB,MXCJB,MYCJB,FLAG2)

C       Calculate the net forces and moments on the impeller by summation
C       of the individual contributions of the thrust bearing and the 
C       conical journal bearing
        FX(N,1) = FXCJB
        FY(N,1) = FYCJB
        FZ(N,1) = FZTB + FZCJB
        MX(N,1) = MXTB + MXCJB
        MY(N,1) = MYTB + MYCJB

C       Write the pressure field at this time step to file (this is done 
C       during the first step of the Runge-Kutta method only) at the first
C       time step and for every twentieth time step thereafter. If you do
C       not want to print the pressure profile to file, then set 
C       (N .EQ. 1) to (N .EQ. 0) and (COUNT .EQ. 20) to (COUNT .EQ. 0) 
C       in the if statement below
        IF ( (N .EQ. 1) .OR. (COUNT .EQ. 40) ) THEN
          CALL PLTCJB(PCJB,N,NTH,NRO,NP,NT,SA,RO,TH,NDT,N1,N2)        
          CALL PLTTB(PTB,N,NTH,NR,NP,NT,R,TH,NDT,N1,N2)
          COUNT = 1
        ELSE
          COUNT = COUNT + 1
        END IF

        K1(1)  = DT*Y2(N)
        K1(2)  = DT*( (LOADP/M1)*FX(N,1) + (V1**2)*XH*SIN(V1*NDT) 
     +           -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K1(3)  = DT*Y4(N)
        K1(4)  = DT*( (LOADP/M1)*FY(N,1) + (V1**2)*YH*SIN(V1*NDT)
     +           -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K1(5)  = DT*Y6(N)
        K1(6)  = DT*( (LOADP/M1)*FZ(N,1) + (V1**2)*ZH*SIN(V1*NDT) 
     +           + COS(PHI1)/M1 )
        K1(7)  = DT*Y8(N)
        K1(8)  = DT*( JP*Y10(N)+(R2/H0TB)*MX(N,1)/IT )
        K1(9)  = DT*Y10(N)
        K1(10) = DT*( -JP*Y8(N)+(R2/H0TB)*MY(N,1)/IT )


C       *** SECOND STEP ***

        NDT = T(N)+0.5*DT       
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT,
     +       N1,N2,Y5(N)+K1(5)/2,Y6(N)+K1(6)/2,Y7(N)+K1(7)/2,Y8(N)+
     +       K1(8)/2,Y9(N)+K1(9)/2,Y10(N)+K1(10)/2,PTB,PPITB,FZTB,MXTB,
     +       MYTB,FLAG1)
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG,
     +       RO,DRO,NRO,NP,NDT,N1,N2,Y1(N)+K1(1)/2,Y2(N)+K1(2)/2,Y3(N)+
     +       K1(3)/2,Y4(N)+K1(4)/2,Y5(N)+K1(5)/2,Y6(N)+K1(6)/2,Y7(N)+
     +       K1(7)/2,Y8(N)+K1(8)/2,Y9(N)+K1(9)/2,Y10(N)+K1(10)/2,PCJB,
     +       PPICJB,FXCJB,FYCJB,FZCJB,MXCJB,MYCJB,FLAG2)
        FX(N,2) = FXCJB
        FY(N,2) = FYCJB
        FZ(N,2) = FZTB + FZCJB
        MX(N,2) = MXTB + MXCJB
        MY(N,2) = MYTB + MYCJB

        K2(1)  = DT*( Y2(N)+0.5*K1(2) )
        K2(2)  = DT*( (LOADP/M1)*FX(N,2) + (V1**2)*XH*SIN(V1*NDT) 
     +           -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K2(3)  = DT*( Y4(N)+0.5*K1(4) )
        K2(4)  = DT*( (LOADP/M1)*FY(N,2) + (V1**2)*YH*SIN(V1*NDT)
     +           -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K2(5)  = DT*( Y6(N)+0.5*K1(6) )
        K2(6)  = DT*( (LOADP/M1)*FZ(N,2) + (V1**2)*ZH*SIN(V1*NDT) 
     +           + COS(PHI1)/M1 )
        K2(7)  = DT*( Y8(N)+0.5*K1(8) )
        K2(8)  = DT*( JP*(Y10(N)+0.5*K1(10))+(R2/H0TB)*MX(N,2)/IT )
        K2(9)  = DT*( Y10(N)+0.5*K1(10) )
        K2(10) = DT*( -JP*(Y8(N)+0.5*K1(8))+(R2/H0TB)*MY(N,2)/IT )


C       *** THIRD STEP ***

        NDT = T(N)+0.5*DT       
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT,
     +       N1,N2,Y5(N)+K2(5)/2,Y6(N)+K2(6)/2,Y7(N)+K2(7)/2,Y8(N)+
     +       K2(8)/2,Y9(N)+K2(9)/2,Y10(N)+K2(10)/2,PTB,PPITB,FZTB,MXTB,
     +       MYTB,FLAG1)
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG,
     +       RO,DRO,NRO,NP,NDT,N1,N2,Y1(N)+K2(1)/2,Y2(N)+K2(2)/2,Y3(N)+
     +       K2(3)/2,Y4(N)+K2(4)/2,Y5(N)+K2(5)/2,Y6(N)+K2(6)/2,Y7(N)+
     +       K2(7)/2,Y8(N)+K2(8)/2,Y9(N)+K2(9)/2,Y10(N)+K2(10)/2,PCJB,
     +       PPICJB,FXCJB,FYCJB,FZCJB,MXCJB,MYCJB,FLAG2)
        FX(N,3) = FXCJB
        FY(N,3) = FYCJB
        FZ(N,3) = FZTB + FZCJB
        MX(N,3) = MXTB + MXCJB
        MY(N,3) = MYTB + MYCJB

        K3(1)  = DT*( Y2(N)+0.5*K2(2) )
        K3(2)  = DT*( (LOADP/M1)*FX(N,3) + (V1**2)*XH*SIN(V1*NDT) 
     +           -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K3(3)  = DT*( Y4(N)+0.5*K2(4) )
        K3(4)  = DT*( (LOADP/M1)*FY(N,3) + (V1**2)*YH*SIN(V1*NDT)
     +           -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K3(5)  = DT*( Y6(N)+0.5*K2(6) )
        K3(6)  = DT*( (LOADP/M1)*FZ(N,3) + (V1**2)*ZH*SIN(V1*NDT) 
     +           + COS(PHI1)/M1 )
        K3(7)  = DT*( Y8(N)+0.5*K2(8) )
        K3(8)  = DT*( JP*(Y10(N)+0.5*K2(10))+(R2/H0TB)*MX(N,3)/IT )
        K3(9)  = DT*( Y10(N)+0.5*K2(10) )
        K3(10) = DT*( -JP*(Y8(N)+0.5*K2(8))+(R2/H0TB)*MY(N,3)/IT )


C       *** FOURTH STEP ***

        NDT = T(N)+DT       
        CALL TB(H0TB,H1TB,KTB,PLTB,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT,
     +       N1,N2,Y5(N)+K3(5),Y6(N)+K3(6),Y7(N)+K3(7),Y8(N)+K3(8),
     +       Y9(N)+K3(9),Y10(N)+K3(10),PTB,PPITB,FZTB,MXTB,MYTB,FLAG1)
        CALL CJB(H0CJB,H1CJB,KCJB,PLCJB,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG,
     +       RO,DRO,NRO,NP,NDT,N1,N2,Y1(N)+K3(1),Y2(N)+K3(2),Y3(N)+
     +       K1(3),Y4(N)+K3(4),Y5(N)+K3(5),Y6(N)+K3(6),Y7(N)+K3(7),
     +       Y8(N)+K3(8),Y9(N)+K3(9),Y10(N)+K3(10),PCJB,PPICJB,FXCJB,
     +       FYCJB,FZCJB,MXCJB,MYCJB,FLAG2)
        FX(N,4) = FXCJB
        FY(N,4) = FYCJB
        FZ(N,4) = FZTB + FZCJB
        MX(N,4) = MXTB + MXCJB
        MY(N,4) = MYTB + MYCJB

        K4(1)  = DT*( Y2(N)+K3(2) )
        K4(2)  = DT*( (LOADP/M1)*FX(N,4) + (V1**2)*XH*SIN(V1*NDT) 
     +           -(M2/M1)*(V2**2)*COS(V2*NDT)+SIN(PHI1)*COS(PHI2)/M1 )
        K4(3)  = DT*( Y4(N)+K3(4) )
        K4(4)  = DT*( (LOADP/M1)*FY(N,4) + (V1**2)*YH*SIN(V1*NDT)
     +           -(M2/M1)*(V2**2)*SIN(V2*NDT)+SIN(PHI1)*SIN(PHI2)/M1 )
        K4(5)  = DT*( Y6(N)+K3(6) )
        K4(6)  = DT*( (LOADP/M1)*FZ(N,4) + (V1**2)*ZH*SIN(V1*NDT) 
     +           + COS(PHI1)/M1 )
        K4(7)  = DT*( Y8(N)+K3(8) )
        K4(8)  = DT*( JP*(Y10(N)+K3(10))+(R2/H0TB)*MX(N,4)/IT )
        K4(9)  = DT*( Y10(N)+K3(10) )
        K4(10) = DT*( -JP*(Y8(N)+K3(8))+(R2/H0TB)*MY(N,4)/IT )


C       Calculate the variables at the next time step based on the
C       results of these four steps. 
        Y1(N+1)  = Y1(N)+(K1(1)+2.*K2(1)+2.*K3(1)+K4(1))/6.
        Y2(N+1)  = Y2(N)+(K1(2)+2.*K2(2)+2.*K3(2)+K4(2))/6.
        Y3(N+1)  = Y3(N)+(K1(3)+2.*K2(3)+2.*K3(3)+K4(3))/6.
        Y4(N+1)  = Y4(N)+(K1(4)+2.*K2(4)+2.*K3(4)+K4(4))/6.
        Y5(N+1)  = Y5(N)+(K1(5)+2.*K2(5)+2.*K3(5)+K4(5))/6.
        Y6(N+1)  = Y6(N)+(K1(6)+2.*K2(6)+2.*K3(6)+K4(6))/6.
        Y7(N+1)  = Y7(N)+(K1(7)+2.*K2(7)+2.*K3(7)+K4(7))/6.
        Y8(N+1)  = Y8(N)+(K1(8)+2.*K2(8)+2.*K3(8)+K4(8))/6.
        Y9(N+1)  = Y9(N)+(K1(9)+2.*K2(9)+2.*K3(9)+K4(9))/6.
        Y10(N+1) = Y10(N)+(K1(10)+2.*K2(10)+2.*K3(10)+K4(10))/6.

C       Write results to the screen
        IF (N .EQ. 1) THEN
          WRITE(*,'(A10,5A12)') 't*','erx','ery','erz','xrot*','yrot*'
        END IF
        WRITE(*,65) T(N)*RTOD,Y1(N),Y3(N),Y5(N),Y7(N),Y9(N)
   65   FORMAT(F10.3,5F12.6)

C       Write results of all variables to file 
        CALL PLTTRA(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9,Y10,DT,T,
     +              N,NT,N3,FLAG1,FLAG2)

C       If touchdown occurs, then print message to screen
        IF (FLAG1 .EQV. (.TRUE.)) THEN
          WRITE(*,'(A40)') 'Touchdown reported by SUBROUTINE TB'
        ELSE IF (FLAG2 .EQV. (.TRUE.)) THEN
          WRITE(*,'(A40)') 'Touchdown reported by SUBROUTINE CJB'
        END IF

C       Increment time for next time step
        N = N + 1

      END DO

C     End and close result files
      END FILE (1)
      END FILE (2)
      END FILE (3)
      CLOSE (1)
      CLOSE (2)
      CLOSE (3)

      STOP
      END
     
      
C ======================================================================

      SUBROUTINE TB(H0,H1,K,PLAND,TH0,TH,DTH,NTH,R1,R2,R,DR,NR,NP,NDT,
     +              N1,N2,ERZ,ERZDT,XROT,XROTDT,YROT,YROTDT,P,PPI,FZT,
     +              MXT,MYT,FLAG2)

C     This subroutine first calculates the pressure field over the thrust
C     bearing surfaces of each impeller blade. This pressure field is then
C     integrated over the respective bearing surface to calculate the
C     forces and moments on the impeller from the thrust bearing alone. 

C     Description of variables (those not already described or different 
C     to those in the main program):
C
C FOR       = Successive over-relaxation factor to accelerate iterative process
C DTOR      = Converts degrees to radians 
C H0        = Minimum film thickness (m)
C H1        = Maximum film thickness (m)
C K         = H1/H0-1 (non-dim.)
C PLAND     = Fraction of angular span of bearing that is untapered land
C X(,,)     = X-coordinate of grid point i,j on impeller blade m (non-dim.)
C Y(,,)     = Y-coordinate of grid point i,j on impeller blade m (non-dim.)
C FZ()      = Force on the impeller in the Z-direction due to the thrust bearing
C             surface of impeller blade m (non-dim.)
C MX()      = Moment about the X-axis on impeller due to the thrust bearing
C             surface of impeller blade m (non-dim.)
C MY()      = Moment about the Y-axis on impeller due to the thrust bearing
C             surface of impeller blade m (non-dim.)
C FZT       = Total force on impeller in the Z-direction due to the thrust 
C             bearing (non-dim.)
C MXT       = Total moment on the impeller about the X-axis due to the thrust
C             bearing (non-dim.)
C MYT       = Total moment on the impeller about the Y-axis due to the thrust
C             bearing (non-dim.)
C TH0A      = Angular span of the land section (rad.)
C TH0B      = Angular span of the tapered section (rad.)
C ERZ       = Eccentricity of the impeller in the Z-direction (non-dim.)
C ERZDT     = Time derivative of ERZ (ie. Z-velocity) (non-dim.)
C XROT      = Rotational position of the impeller about the X-axis (non-dim.)
C XROTDT    = Time derivative of XROT (ie. X rotational velocity) (non-dim.)
C XROT      = Rotational position of the impeller about the Y-axis (non-dim.)
C XROTDT    = Time derivative of YROT (ie. Y rotational velocity) (non-dim.)
C A()-D()   = Coefficients for the Thomas algorithm subroutine
C H(,,)     = Film thickness at any point (non-dim.)
C DHDTH(,,) = Derivative of H with respect to theta (non-dim.)
C DHDR(,,)  = Derivative of H with respect to r (non-dim.)
C DHDT(,,)  = Derivative of H with respect to time (non-dim.)
C P(,,)     = Pressure at node i,j on impeller blade m at time step N (non-dim.)
C PPI(,,)   = Pressure at node i,j on impeller blade m at time step N-1 (non-dim.) 
C FLAG1     = Set to true if solution converges, otherwise false
C FLAG2     = Set to true if touchdown occurs, otherwise false

C     Declaration of variables
      PARAMETER (PI=3.141592654, DTOR=PI/180.0, FOR=1.5)
      REAL H0, H1, K, R1, R2, LR, NDR, PLAND, TH0D 
      REAL X(N1,N1,N2), Y(N1,N1,N2), DBLINT
      REAL PR(N1,N1,N2), PRX(N1,N1,N2), PRY(N1,N1,N2) 
      REAL FZ(N2), MX(N2), MY(N2), FZT, MXT, MYT
      DOUBLE PRECISION TH(N1,N2), TH0, TH0A, TH0B, TH1, R(N1)
      DOUBLE PRECISION DTH, DTH2, DR, DR2, NDT
      DOUBLE PRECISION ERZ, ERZDT, XROT, XROTDT, YROT, YROTDT 
      DOUBLE PRECISION A(N1), B(N1), C(N1), D(N1)
      DOUBLE PRECISION H(N1,N1,N2), DHDTH(N1,N1,N2), DHDR(N1,N1,N2)
      DOUBLE PRECISION DHDT(N1,N1,N2), P(N1,N1,N2), PPI(N1,N1,N2)
      INTEGER NTH, NTHM1, NTHM2, NR, NRM1, NRM2
      INTEGER I, J, M, N, NP, N1, N2
      LOGICAL FLAG1, FLAG2

C     Initialise dependent variables
      NTHM1 = NTH-1
      NTHM2 = NTH-2
      NRM1  = NR-1
      NRM2  = NR-2
      TH0A  = TH0*PLAND
      TH0B  = TH0-TH0A
      TH1   = 2*PI/NP
      LR    = R2-R1
      NDR   = (R2/LR)**2
      DTH2  = DTH**2
      DR2   = DR**2

C     Initialise variables
      FZT = 0.0
      MXT = 0.0
      MYT = 0.0

C     Calculate the pressure distribution and subsequently the forces and moments 
C     on the impeller as a result of this pressure over each pad m=1,2,...,Np                
      DO 65,M=1,NP
               
C       Calculate the local film thickness and spatial and temporal derivatives
        DO 10,J=1,NR
          DO 5,I=1,NTH

C           Film thickness for untapered land
            H(I,J,M)= 1.0+ERZ
     +                +XROT*R(J)*SIN(TH(I,M)+NDT)
     +                -YROT*R(J)*COS(TH(I,M)+NDT)  

C           Derivative with respect to (w.r.t.) theta for untapered land
            DHDTH(I,J,M) = XROT*R(J)*COS(TH(I,M)+NDT)
     +                    +YROT*R(J)*SIN(TH(I,M)+NDT)

C           Film thickess and derivative w.r.t. theta for taper
            IF ((TH(I,M)-(M-1)*TH1) .GE. TH0A) THEN
              H(I,J,M) = H(I,J,M)+K*(TH(I,M)-(TH0A+(M-1)*TH1))/TH0B
              DHDTH(I,J,M) = DHDTH(I,J,M) + K/TH0B
            END IF           

C           If the film thickness is less than or equal to zero, then
C           touchdown has occurred. Set FLAG2 to true.
            IF (H(I,J,M) .LE. 0.0D0) THEN
              FLAG2 = .TRUE.
              RETURN
            END IF

C           Radial and time derivatives
            DHDR(I,J,M) = XROT*SIN(TH(I,M)+NDT)
     +                   -YROT*COS(TH(I,M)+NDT)
            DHDT(I,J,M) = ERZDT
     +                    +(XROT-YROTDT)*R(J)*COS(TH(I,M)+NDT)
     +                    +(XROTDT+YROT)*R(J)*SIN(TH(I,M)+NDT)
    5     CONTINUE
   10   CONTINUE        

C       Initialise iterative process loop control variable FLAG1
        FLAG1 = .FALSE.

C       Iterative process to solve for the pressure distribution

C --------------- BEGIN ITERATIVE PROCESS ---------------
        
        DO WHILE(FLAG1 .EQV. (.FALSE.))

C         Sweep all rows j=2,3,...,J-1
          DO 25,J=2,NRM1
            DO 15,I=2,NTHM1
              A(I-1)=1./(R(J)*DTH2)-1.5*DHDTH(I,J,M)/(R(J)*H(I,J,M)*DTH)
              B(I-1)=-2.*(R(J)/DR2+1./(R(J)*DTH2))
              C(I-1)=1./(R(J)*DTH2)+1.5*DHDTH(I,J,M)/(R(J)*H(I,J,M)*DTH)
              D(I-1)=-(NDR*R(J))/(H(I,J,M)**3)*DHDTH(I,J,M)
     +               +(2.*NDR*R(J))/(H(I,J,M)**3)*DHDT(I,J,M)
     +               -(R(J)/DR2-(1.5*R(J))/(H(I,J,M)*DR)
     +               *DHDR(I,J,M)-0.5/DR)*P(I,J-1,M)
     +               -(R(J)/DR2+(1.5*R(J))/(H(I,J,M)*DR)
     +               *DHDR(I,J,M)+0.5/DR)*P(I,J+1,M)       
   15       CONTINUE
            
C           Use Thomas algorithm to solve tridiagonal system of algebraic equations
            CALL THOMAS(NTHM2,A,B,C,D,N1)
            
C           Capture solution from Thomas and apply SOR factor, For
            DO 20,I=2,NTHM1
              P(I,J,M)   = PPI(I,J,M)+FOR*(D(I-1)-PPI(I,J,M))             
              PPI(I,J,M) = P(I,J,M)
   20       CONTINUE
         
   25     CONTINUE
         
C         Sweep all columns i=2,3,...,I-1
          DO 40,I=2,NTHM1
            DO 30,J=2,NRM1
              A(J-1) = (R(J)/DR-(1.5*R(J)*DHDR(I,J,M))/H(I,J,M)-0.5)/DR
              B(J-1) = -2.*(R(J)/DR2+1./(R(J)*DTH2))
              C(J-1) = (R(J)/DR+(1.5*R(J)*DHDR(I,J,M))/H(I,J,M)+0.5)/DR
              D(J-1) = -(NDR*R(J))/(H(I,J,M)**3)*DHDTH(I,J,M)
     +                 +(2.*NDR*R(J))/(H(I,J,M)**3)*DHDT(I,J,M)
     +                 -(1./(R(J)*DTH2)-1.5/(R(J)*H(I,J,M)*DTH)
     +                 *DHDTH(I,J,M))*P(I-1,J,M)
     +                 -(1./(R(J)*DTH2)+1.5/(R(J)*H(I,J,M)*DTH)
     +                 *DHDTH(I,J,M))*P(I+1,J,M)
   30       CONTINUE
            
C           Use Thomas algorithm to solve tridiagonal system of algebraic equations
            CALL THOMAS(NRM2,A,B,C,D,N1)
            
C           Capture solution from Thomas and apply SOR factor, For
            DO 35,J=2,NRM1
              P(I,J,M) = PPI(I,J,M)+FOR*(D(J-1)-PPI(I,J,M))
   35       CONTINUE
            
   40     CONTINUE
         
C         Convergence test
          CALL CVERGE(NTH,NR,P,PPI,FLAG1,M,N1,N2)
          
        END DO
      
C --------------- END ITERATIVE PROCESS ---------------


C       Transform the position of each grid point on the lubricating surfaces of 
C       the rotating impeller in Cartesian x-y coordinates relative to the pump 
C       housing (which is fixed)
        DO 50,J=1,NR
          DO 45,I=1,NTH
            X(I,J,M)=R(J)*COS(TH(I,M)+NDT)
            Y(I,J,M)=R(J)*SIN(TH(I,M)+NDT)
   45     CONTINUE
   50   CONTINUE

C       Calculate the force (Fz) and moments (Mx and My) on the impeller 
C       due to each separate blade
        DO 60,J=1,NR
          DO 55,I=1,NTH
            PR(I,J,M)  = P(I,J,M)*R(J)
            PRX(I,J,M) = -P(I,J,M)*R(J)*X(I,J,M)
            PRY(I,J,M) = P(I,J,M)*R(J)*Y(I,J,M)
   55     CONTINUE
   60   CONTINUE
        FZ(M) = DBLINT(PR,DTH,DR,NTHM1,NRM1,N1,N2,M)
        MX(M) = DBLINT(PRY,DTH,DR,NTHM1,NRM1,N1,N2,M) 
        MY(M) = DBLINT(PRX,DTH,DR,NTHM1,NRM1,N1,N2,M)

C       Calculate the total force (Fz) and moment (Mx and My) on the centre
C       of mass of the impeller
        FZT = FZT + FZ(M)  
        MXT = MXT + MX(M)
        MYT = MYT + MY(M)

   65 CONTINUE

      RETURN
      END


C ======================================================================

      SUBROUTINE CJB(H0,H1,K,PLAND,ALP,TH0,TH,DTH,NTH,RO1,RO2,ROG,RO,
     +               DRO,NRO,NP,NDT,N1,N2,ERX,ERXDT,ERY,ERYDT,ERZ,ERZDT,
     +               XROT,XROTDT,YROT,YROTDT,P,PPI,FXT,FYT,FZT,MXT,MYT,
     +               FLAG2)

C     This subroutine first calculates the pressure field over the conical
C     journal bearing surfaces of each impeller blade. This pressure field 
C     is then integrated over the respective bearing surface to calculate
C     the forces and moments on the impeller from the conical journal bearing 
C     alone. 

C     Description of variables (those not already described or different 
C     to those in the main program):
C
C FOR       = Successive over-relaxation factor to accelerate iterative process
C DTOR      = Converts degrees to radians 
C H0        = Minimum film thickness (m)
C H1        = Maximum film thickness (m)
C K         = H1/H0-1 (non-dim.)
C PLAND     = Fraction of angular span of bearing that is untapered land
C FX()      = Force on the impeller in the X-direction due to the conical 
C             journal bearing surface of impeller blade m (non-dim.)
C FY()      = Force on the impeller in the Y-direction due to the conical
C             journal bearing surface of impeller blade m (non-dim.)
C FZ()      = Force on the impeller in the Z-direction due to the conical
C             journal bearing surface of impeller blade m (non-dim.)
C MX()      = Moment about the X-axis on impeller due to the conical journal
C             bearing surface of impeller blade m (non-dim.)
C MY()      = Moment about the Y-axis on impeller due to the conical journal
C             bearing surface of impeller blade m (non-dim.)
C FXT       = Total force on impeller in the X-direction due to the conical 
C             journal bearing (non-dim.)
C FYT       = Total force on impeller in the Y-direction due to the conical
C             journal bearing (non-dim.)
C FZT       = Total force on impeller in the Z-direction due to the conical
C             journal bearing (non-dim.)
C MXT       = Total moment on the impeller about the X-axis due to the conical
C             journal bearing (non-dim.)
C MYT       = Total moment on the impeller about the Y-axis due to the conical
C             journal bearing (non-dim.)
C TH0A      = Angular span of the land section (rad.)
C TH0B      = Angular span of the tapered section (rad.)
C ERX       = Eccentricity of the impeller in the X-direction (non-dim.)
C ERXDT     = Time derivative of ERX (ie. X-velocity) (non-dim.)
C ERY       = Eccentricity of the impeller in the Y-direction (non-dim.)
C ERYDT     = Time derivative of ERY (ie. Y-velocity) (non-dim.)
C ERZ       = Eccentricity of the impeller in the Z-direction (non-dim.)
C ERZDT     = Time derivative of ERZ (ie. Z-velocity) (non-dim.)
C XROT      = Rotational position of the impeller about the X-axis (non-dim.)
C XROTDT    = Time derivative of XROT (ie. X rotational velocity) (non-dim.)
C XROT      = Rotational position of the impeller about the Y-axis (non-dim.)
C XROTDT    = Time derivative of YROT (ie. Y rotational velocity) (non-dim.)
C A()-D()   = Coefficients for the Thomas algorithm subroutine
C H(,,)     = Film thickness at any point (non-dim.)
C DHDTH(,,) = Derivative of H with respect to theta (non-dim.)
C DHDRO(,,) = Derivative of H with respect to r (non-dim.)
C DHDT(,,)  = Derivative of H with respect to time (non-dim.)
C P(,,)     = Pressure at node i,j on impeller blade m at time step N (non-dim.)
C PPI(,,)   = Pressure at node i,j on impeller blade m at time step N-1 (non-dim.) 
C FLAG1     = Set to true if solution converges, otherwise false
C FLAG2     = Set to true if touchdown occurs, otherwise false

C     Declaration of variables
      PARAMETER (PI=3.141592654, DTOR=PI/180.0, FOR = 1.5)
      REAL H0, H1, K, RO1, RO2, ROG, LRO, NDR, PLAND
      REAL TH0D, ALPD, DBLINT
      REAL FX(N2), FY(N2), FZ(N2), MX(N2), MY(N2)
      REAL FXT, FYT, FZT, MXT, MYT
      REAL DFX(N1,N1,N2), DFY(N1,N1,N2), DFZ(N1,N1,N2)
      REAL DMX(N1,N1,N2), DMY(N1,N1,N2)
      DOUBLE PRECISION TH(N1,N2), TH0, TH0A, TH0B, TH1, RO(N1)
      DOUBLE PRECISION DTH, DTH2, DRO, DRO2
      DOUBLE PRECISION ERX, ERXDT, ERY, ERYDT, ERZ, ERZDT
      DOUBLE PRECISION XROT, XROTDT, YROT, YROTDT
      DOUBLE PRECISION A(N1), B(N1), C(N1), D(N1), NDT, ALP, SA, CA
      DOUBLE PRECISION H(N1,N1,N2), DHDRO(N1,N1,N2), DHDTH(N1,N1,N2)
      DOUBLE PRECISION DHDT(N1,N1,N2), P(N1,N1,N2), PPI(N1,N1,N2)
      INTEGER NTH, NRO, NTHM1, NTHM2, NROM1, NROM2
      INTEGER I, J, M, NP, N1, N2
      LOGICAL FLAG1, FLAG2

C     Independent variables
      SA   = DSIN(ALP)
      CA   = DCOS(ALP)

C     Dependent variables
      NTHM1 = NTH-1
      NTHM2 = NTH-2
      NROM1 = NRO-1
      NROM2 = NRO-2
      TH0A  = TH0*PLAND
      TH0B  = TH0-TH0A
      TH1   = 2*PI/NP
      LRO   = RO2-RO1
      NDR   = (RO2/LRO)**2
      DTH2  = DTH**2
      DRO2  = DRO**2

C     Initialise variables
      FXT = 0.0
      FYT = 0.0
      FZT = 0.0
      MXT = 0.0
      MYT = 0.0

C     Calculate the pressure distribution and subsequently the forces and moments 
C     on the impeller as a result of this pressure over each pad m=1,2,...,Np       
      DO 55,M=1,NP
      
C       Calculate local film thickness and the spatial and temporal derivatives 
        DO 10,J=1,NRO
          DO 5,I=1,NTH
            H(I,J,M) = ( 1.0-ERZ
     +                 -ERX*COS(TH(I,M)+NDT)
     +                 -ERY*SIN(TH(I,M)+NDT) )*SA                         
     +                 -XROT*(RO(J)-ROG)*SIN(TH(I,M)+NDT) 
     +                 +YROT*(RO(J)-ROG)*COS(TH(I,M)+NDT) 
            DHDTH(I,J,M) = ( ERX*SIN(TH(I,M)+NDT)                
     +                    -ERY*COS(TH(I,M)+NDT) )*SA                         
     +                    -XROT*(RO(J)-ROG)*COS(TH(I,M)+NDT)
     +                    -YROT*(RO(J)-ROG)*SIN(TH(I,M)+NDT) 
            IF ((TH(I,M)-(M-1)*TH1) .GE. TH0A) THEN           
              H(I,J,M) = H(I,J,M) + 
     +                   ( K*(TH(I,M)-(TH0A+(M-1)*TH1))/TH0B )*SA                  
              DHDTH(I,J,M) = DHDTH(I,J,M) + ( K/TH0B )*SA                                    
            END IF                                                       
            IF (H(I,J,M) .LE. 0.0D0) THEN
              FLAG2 = .TRUE.
              RETURN
            END IF
            DHDRO(I,J,M) = -XROT*SIN(TH(I,M)+NDT)            
     +                     +YROT*COS(TH(I,M)+NDT)            
            DHDT(I,J,M)  = ( -ERZDT
     +                     +(ERX-ERYDT)*SIN(TH(I,M)+NDT)           
     +                     -(ERXDT+ERY)*COS(TH(I,M)+NDT) )*SA                         
     +                     -(XROT-YROTDT)*(RO(J)-ROG)*COS(TH(I,M)+NDT)
     +                     -(XROTDT+YROT)*(RO(J)-ROG)*SIN(TH(I,M)+NDT)
    5     CONTINUE
   10   CONTINUE

C       Initialise iterative process loop control variable FLAG1
        FLAG1 = .FALSE.


C       Iterative process to solve for the pressure distribution

C --------------- BEGIN ITERATIVE PROCESS ---------------
        
        DO WHILE(FLAG1 .EQV. (.FALSE.))
        
C         Sweep all rows j=2,3,...,J-1
          DO 25,J=2,NROM1
            DO 15,I=2,NTHM1
              A(I-1) = (1./DTH-1.5/H(I,J,M)*DHDTH(I,J,M))
     +                 /(DTH*(RO(J)*SA)**2)
              B(I-1) = -2*(1./DRO2+1./(DTH*RO(J)*SA)**2)
              C(I-1) = (1/DTH+1.5/H(I,J,M)*DHDTH(I,J,M))
     +                 /(DTH*(RO(J)*SA)**2)
              D(I-1) = (2.0*NDR*DHDT(I,J,M))/H(I,J,M)**3
     +                 -NDR*DHDTH(I,J,M)/H(I,J,M)**3
     +                 -(1./DRO)*(1./DRO-1.5/H(I,J,M)*DHDRO(I,J,M)
     +                 -0.5/RO(J))*P(I,J-1,M)
     +                 -(1./DRO)*(1./DRO+1.5/H(I,J,M)*DHDRO(I,J,M)
     +                 +0.5/RO(J))*P(I,J+1,M)
   15       CONTINUE

C           Use Thomas algorithm to solve tridiagonal system of algebraic equations  
            CALL THOMAS(NTHM2,A,B,C,D,N1)
            
C           Capture solution from Thomas and apply SOR factor, For
            DO 20,I=2,NTHM1
              P(I,J,M)   = PPI(I,J,M)+FOR*(D(I-1)-PPI(I,J,M))
              PPI(I,J,M) = P(I,J,M)
   20       CONTINUE
          
   25     CONTINUE

          
C         Sweep all columns i=2,3,...,I-1
          DO 40,I=2,NTHM1
            DO 30,J=2,NROM1
              A(J-1) = (1./DRO)*(1./DRO-1.5/H(I,J,M)
     +                 *DHDRO(I,J,M)-0.5/RO(J))  
              B(J-1) = -2*(1./DRO2+1./(DTH*RO(J)*SA)**2)
              C(J-1) = (1./DRO)*(1./DRO+1.5/H(I,J,M)
     +                 *DHDRO(I,J,M)+0.5/RO(J))
              D(J-1) = (2.0*NDR*DHDT(I,J,M))/H(I,J,M)**3
     +                -NDR*DHDTH(I,J,M)/H(I,J,M)**3
     +                -(1./DTH-1.5/H(I,J,M)*DHDTH(I,J,M))
     +                 /(DTH*(RO(J)*SA)**2)*P(I-1,J,M)
     +                -(1./DTH+1.5/H(I,J,M)*DHDTH(I,J,M))
     +                 /(DTH*(RO(J)*SA)**2)*P(I+1,J,M)
   30       CONTINUE
            
C           Use Thomas algorithm to solve tridiagonal system of algebraic equations
            CALL THOMAS(NROM2,A,B,C,D,N1)
            
C           Capture solution from Thomas and apply SOR factor, For
            DO 35,J=2,NROM1
              P(I,J,M) = PPI(I,J,M)+FOR*(D(J-1)-PPI(I,J,M))
   35       CONTINUE
            
   40     CONTINUE

C         Convergence test
          CALL CVERGE(NTH,NRO,P,PPI,FLAG1,M,N1,N2)


        END DO       
      
C --------------- END ITERATIVE PROCESS ---------------


C       Calculate forces (Fx, Fy and Fz) and moments (Mx and My) on 
C       the impeller due to each separate blade
        DO 50,J=1,NRO
          DO 45,I=1,NTH
            DFX(I,J,M)=P(I,J,M)*RO(J)*COS(TH(I,M)+NDT)
            DFY(I,J,M)=P(I,J,M)*RO(J)*SIN(TH(I,M)+NDT) 
            DFZ(I,J,M)=P(I,J,M)*RO(J)
            DMX(I,J,M)=P(I,J,M)*RO(J)*(RO(J)-ROG)*SIN(TH(I,M)+NDT)
            DMY(I,J,M)=P(I,J,M)*RO(J)*(RO(J)-ROG)*COS(TH(I,M)+NDT)
   45     CONTINUE
   50   CONTINUE
        FX(M) = -SA*CA*DBLINT(DFX,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        FY(M) = -SA*CA*DBLINT(DFY,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        FZ(M) = -SA*SA*DBLINT(DFZ,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        MX(M) = -SA*DBLINT(DMX,DRO,DTH,NROM1,NTHM1,N1,N2,M)
        MY(M) = SA*DBLINT(DMY,DRO,DTH,NROM1,NTHM1,N1,N2,M)

C       Calculate the total force (Fx, Fy and Fz) and moment (Mx and My) 
C       on the centre of mass of the impeller
        FXT = FXT + FX(M)  
        FYT = FYT + FY(M)
        FZT = FZT + FZ(M)
        MXT = MXT + MX(M)
        MYT = MYT + MY(M)

   55 CONTINUE


      RETURN
      END


C ======================================================================

      SUBROUTINE THOMAS(NXM2,A,B,C,D,N1)

C     This subroutine uses the Thomas algorithm for solving tri-diagonal
C     systems of algebraic equations with Dirichlet boundary conditions.  
C     It is used in iterative line-by-line sweeps to determine the pressure
C     distribution over the bearing surfaces. 

      DOUBLE PRECISION A(N1),B(N1),C(N1),D(N1)
      INTEGER I, N1, NXM2

      B(1)=1.0/B(1)

      DO 5,I=2,NXM2
        B(I)=1.0/(B(I)-A(I)*B(I-1)*C(I-1))
        D(I)=D(I)-A(I)*B(I-1)*D(I-1)
    5 CONTINUE

      D(NXM2)=D(NXM2)*B(NXM2)

      DO 10,I=NXM2-1,1,-1
        D(I)=(D(I)-D(I+1)*C(I))*B(I)
   10 CONTINUE

      RETURN
      END


C ======================================================================

      SUBROUTINE CVERGE(NX,NY,P,PPI,FLAG,M,N1,N2)

C     This subroutine tests for convergence of the solution. If the 
C     residual is less than or equal to RESLIM then the solution is 
C     converged.  If the residual is greater than RESMAX, the solution 
C     is diverging and the program is stopped.  Otherwise, an over -
C     relaxation factor is applied to speed up the iterative process 
C     (This is performed in subroutines TB and CJB)

      PARAMETER (RESLIM=1.0D-4, RESMAX=1.0D4)
      DOUBLE PRECISION P(N1,N1,N2), PPI(N1,N1,N2), RESID
      LOGICAL FLAG

      RESID=0.0D0
      DO 10,J=1,NY
        DO 5,I=1,NX
          RESID = MAX(RESID,ABS(P(I,J,M)-PPI(I,J,M)))
    5 CONTINUE
   10 CONTINUE
      IF (RESID .LE. RESLIM) THEN
        FLAG = .TRUE.
      ELSE IF (RESID .GE. RESMAX) THEN
        STOP 'Solution not converging'
      ELSE
        DO 20,J=2,NY-1
          DO 15,I=2,NX-1
            PPI(I,J,M) = P(I,J,M)          
   15     CONTINUE
   20   CONTINUE
      END IF

      RETURN
      END


C ======================================================================

      REAL FUNCTION DBLINT(FXY,DX,DY,NXM1,NYM1,N1,N2,M)

C     This subroutine performs a double integration using the
C     trapezoidal rule of integration.

      REAL SUM, FXY(N1,N1,N2)
      DOUBLE PRECISION DX, DY
      INTEGER I, J, M, NXM1, NYM1

      SUM = 0.0
      DO 10,J=1,NYM1
        DO 5,I=1,NXM1
          SUM = SUM + FXY(I,J,M) + FXY(I,J+1,M) +
     +          FXY(I+1,J,M) + FXY(I+1,J+1,M)
    5   CONTINUE
   10 CONTINUE

      DBLINT = 0.25*DX*DY*SUM

      END


C ======================================================================

      SUBROUTINE PLTTRA(ERX,ERXDT,ERY,ERYDT,ERZ,ERZDT,XROT,XROTDT,YROT,
     +                  YROTDT,DT,T,N,NT,N3,FLAG1,FLAG2)

C     This subroutine creates a Matlab script file that plots the
C     trajectory of the impeller.  It doubles as a results file.

      PARAMETER(PI=3.141592654, RTOD=180.0/PI)
      DOUBLE PRECISION ERX(N3), ERY(N3), ERZ(N3), XROT(N3), YROT(N3)
      DOUBLE PRECISION ERXDT(N3), ERYDT(N3), ERZDT(N3)
      DOUBLE PRECISION XROTDT(N3), YROTDT(N3), DT, T(N3)
      INTEGER N, NT, N3
      LOGICAL FLAG1, FLAG2

      IF (N .EQ. 1) THEN 
        WRITE(1,*)'% Where C refers to the column no. of matrix data:'
        WRITE(1,*)'% C1:  Non-dim. time'
        WRITE(1,*)'% C2:  Eccentricity ratio in the X-direction'
        WRITE(1,*)'% C3:  Translational velocity in the X-direction'
        WRITE(1,*)'% C4:  Eccentricity ratio in the Y-direction'
        WRITE(1,*)'% C5:  Translational velocity in the Y-direction'
        WRITE(1,*)'% C6:  Eccentricity ratio in the Z-direction'
        WRITE(1,*)'% C7:  Translational velocity in the X-direction'
        WRITE(1,*)'% C8:  Non-dim. rotation about the X-axis'
        WRITE(1,*)'% C9:  Non-dim. rotational velocity about the X-axis'
        WRITE(1,*)'% C10: Non-dim. rotation about the X-axis'
        WRITE(1,*)'% C11: Non-dim. rotational velocity about the Y-axis'
        WRITE(1,*)
        WRITE(1,*) 'data=['
      END IF

      WRITE(1,5) T(N)*RTOD,ERX(N),ERXDT(N),ERY(N),ERYDT(N),ERZ(N), 
     +           ERZDT(N),XROT(N),XROTDT(N),YROT(N),YROTDT(N)

      IF ((N .EQ. NT) .OR. (FLAG1 .EQV. (.TRUE.))
     +                  .OR. (FLAG2 .EQV. (.TRUE.))) THEN
        WRITE(1,*) '];'
        WRITE(1,*)
        WRITE(1,*) 'sf=', 2*PI/DT, ';'
        WRITE(1,*) 'K=', K, ';'
        WRITE(1,*) 'erx=data(:,2);' 
        WRITE(1,*) 'ery=data(:,4);'
        WRITE(1,*) 'erz=data(:,6);'
        WRITE(1,*) 'xrot=data(:,8);'
        WRITE(1,*) 'yrot=data(:,10);'
C       A correction is made to the eccentricities in the X and Y
C       directions where the erz term below accounts for increased / 
C       decreased clearance as the impeller moves down/up inside the
C       cone shaped upper pump housing.
        WRITE(1,*) 'erx=erx./(1.0-erz);'
        WRITE(1,*) 'ery=ery./(1.0-erz);'
      END IF
    5 FORMAT(F10.3,10F12.6)

      RETURN
      END

      
C ======================================================================

      SUBROUTINE PLTCJB(P,N,NTH,NRO,NP,NT,SA,RO,TH,NDT,N1,N2)

C     This subroutine plots the pressure distribution over the conical
C     journal bearing of the impeller at each time step.  A Matlab movie
C     is created by playing these plots in the correct sequence.  To run
C     this movie load and run the Matlab script file "Orbit3D.m" and
C     click the appropriate selection on the pop-up menu.

      DOUBLE PRECISION P(N1,N1,N2), TH(N1,N2), RO(N1), NDT, SA
      INTEGER I, J, M, N, NTH, NRO, NP, NT 

      DO 65,M=1,NP

        IF ((N .EQ. 1) .AND. (M .EQ. 1)) THEN
          WRITE(2,*)  'dth=2*pi/100;'
          WRITE(2,*)  'theta=0:dth:2*pi;'
          WRITE(2,40) 'xro1=', RO(1),   '*cos(theta);'
          WRITE(2,40) 'yro1=', RO(1),   '*sin(theta);'
          WRITE(2,40) 'xro2=', RO(NRO), '*cos(theta);'
          WRITE(2,40) 'yro2=', RO(NRO), '*sin(theta);'
          WRITE(2,*)  'MOV1=moviein(', NT, ');'
          WRITE(2,*)
        END IF

        IF (N .LT. 10) THEN
          WRITE(2,10) 'xcjb', N, M, '=['
        ELSE IF (N .LT. 100) THEN
          WRITE(2,15) 'xcjb', N, M, '=['
        ELSE
          WRITE(2,20) 'xcjb', N, M, '=['
        END IF
        DO 5,J=1,NRO
          WRITE(2,25) ((RO(J)*COS(TH(I,M)+NDT)),I=1,NTH)
    5   CONTINUE
   10   FORMAT(A4,I1,I1,A2)
   15   FORMAT(A4,I2,I1,A2)
   20   FORMAT(A4,I3,I1,A2)
   25   FORMAT(101F10.5)
        WRITE(2,*) '];'
        WRITE(2,*)

        IF (N .LT. 10) THEN
          WRITE(2,10) 'ycjb', N, M, '=['
        ELSE IF (N .LT. 100) THEN
          WRITE(2,15) 'ycjb', N, M, '=['
        ELSE
          WRITE(2,20) 'ycjb', N, M, '=['
        END IF
        DO 30,J=1,NRO
          WRITE(2,25) ((RO(J)*SIN(TH(I,M)+NDT)),I=1,NTH)
   30   CONTINUE
        WRITE(2,*) '];'
        WRITE(2,*)

        IF (N .LT. 10) THEN
          WRITE(2,10) 'pcjb', N, M, '=['
        ELSE IF (N .LT. 100) THEN
          WRITE(2,15) 'pcjb', N, M, '=['
        ELSE
          WRITE(2,20) 'pcjb', N, M, '=['
        END IF
        DO 35,J=1,NRO
          WRITE(2,25) (P(I,J,M),I=1,NTH)
   35   CONTINUE
        WRITE(2,*) '];'
        WRITE(2,*)     

        IF (N .LT. 10) THEN
          WRITE(2,45) 'surf(xcjb',N,M,',ycjb',N,M,',pcjb',N,M,');'
        ELSE IF (N .LT. 100) THEN
          WRITE(2,50) 'surf(xcjb',N,M,',ycjb',N,M,',pcjb',N,M,');' 
        ELSE
          WRITE(2,55) 'surf(xcjb',N,M,',ycjb',N,M,',pcjb',N,M,');' 
        END IF
        IF (M .EQ. 1) THEN           
          WRITE(2,*) 'hold on;'
          WRITE(2,60) 'axis([',-1./SA,1./SA,-1./SA,1./SA,0,0.5,']);'
          WRITE(2,*) 'set(gca,''PlotBoxAspectRatio'',[2 2 1]);'
          WRITE(2,*) 'view(25,30);'
          WRITE(2,*) 'plot(xro1,yro1,''--'');'
          WRITE(2,*) 'plot(xro2,yro2,''--'');' 
          WRITE(2,*) 'xlabel(''X / sin\ alpha'');'   
          WRITE(2,*) 'ylabel(''Y / sin\ alpha'');'
          WRITE(2,*) 'zlabel(''Pressure, p^*'');'
          WRITE(2,*)
        ELSE IF (M .EQ. NP) THEN
          WRITE(2,*) 'shading interp;'
          WRITE(2,*) 'hold off;'
          WRITE(2,*) 'MOV1(', N ,') = getframe;'
          WRITE(2,*)
        END IF
   40   FORMAT(A5,F6.4,A12)
   45   FORMAT(A9,2I1,A5,2I1,A5,2I1,A2)
   50   FORMAT(A9,I2,I1,A5,I2,I1,A5,I2,I1,A2)
   55   FORMAT(A9,I3,I1,A5,I3,I1,A5,I3,I1,A2)
   60   FORMAT(A6,6F10.6,A3)

   65 CONTINUE

      RETURN    
      END 


C ======================================================================

      SUBROUTINE PLTTB(P,N,NTH,NR,NP,NT,R,TH,NDT,N1,N2)

C     This subroutine plots the pressure distribution over the thrust
C     bearing of the VentrAssist impeller at each time step.  A Matlab 
C     movie is created by playing these plots in the correct sequence.
C     To run this movie load and run the Matlab script file "Orbit3D.m" 
C     and click the appropriate selection on the pop-up menu.

      DOUBLE PRECISION P(N1,N1,N2), TH(N1,N2), R(N1), NDT
      INTEGER I, J, M, N, NTH, NR, NP, NT 

      DO 60,M=1,NP

        IF ((N .EQ. 1) .AND. (M .EQ. 1)) THEN
          WRITE(3,*) 'dth=2*pi/100;'
          WRITE(3,*) 'theta=0:dth:2*pi;'
          WRITE(3,40) 'xr1=', R(1),  '*cos(theta);'
          WRITE(3,40) 'yr1=', R(1),  '*sin(theta);'
          WRITE(3,40) 'xr2=', R(NR), '*cos(theta);'
          WRITE(3,40) 'yr2=', R(NR), '*sin(theta);'
          WRITE(3,*)  'MOV2=moviein(', NT, ');'
          WRITE(3,*)
        END IF

        IF (N .LT. 10) THEN
          WRITE(3,10) 'xtb', N, M, '=['
        ELSE IF (N .LT. 100) THEN
          WRITE(3,15) 'xtb', N, M, '=['
        ELSE
          WRITE(3,20) 'xtb', N, M, '=['
        END IF
        DO 5,J=1,NR
          WRITE(3,25) ((R(J)*COS(TH(I,M)+NDT)),I=1,NTH)
    5   CONTINUE
   10   FORMAT(A4,I1,I1,A2)
   15   FORMAT(A4,I2,I1,A2)
   20   FORMAT(A4,I3,I1,A2)
   25   FORMAT(101F10.5)
        WRITE(3,*) '];'
        WRITE(3,*)

        IF (N .LT. 10) THEN
          WRITE(3,10) 'ytb', N, M, '=['
        ELSE IF (N .LT. 100) THEN
          WRITE(3,15) 'ytb', N, M, '=['
        ELSE
          WRITE(3,20) 'ytb', N, M, '=['
        END IF
        DO 30,J=1,NR
          WRITE(3,25) ((R(J)*SIN(TH(I,M)+NDT)),I=1,NTH)
   30   CONTINUE
        WRITE(3,*) '];'
        WRITE(3,*)

        IF (N .LT. 10) THEN
          WRITE(3,10) 'ptb', N, M, '=['
        ELSE IF (N .LT. 100) THEN
          WRITE(3,15) 'ptb', N, M, '=['
        ELSE
          WRITE(3,20) 'ptb', N, M, '=['
        END IF
        DO 35,J=1,NR
          WRITE(3,25) (P(I,J,M),I=1,NTH)
   35   CONTINUE
        WRITE(3,*) '];'
        WRITE(3,*)     
     
        IF (N .LT. 10) THEN
          WRITE(3,45) 'surf(xtb',N,M,',ytb',N,M,',ptb',N,M,');'
        ELSE IF (N .LT. 100) THEN
          WRITE(3,50) 'surf(xtb',N,M,',ytb',N,M,',ptb',N,M,');' 
        ELSE
          WRITE(3,55) 'surf(xtb',N,M,',ytb',N,M,',ptb',N,M,');' 
        END IF
        IF (M .EQ. 1) THEN
          WRITE(3,*) 'hold on;'
          WRITE(3,*) 'axis([-1 1 -1 1 0 0.5]);'
          WRITE(3,*) 'set(gca,''PlotBoxAspectRatio'',[2 2 1]);'
          WRITE(3,*) 'view(25,30);'
          WRITE(3,*) 'plot(xr1,yr1,''--'');'
          WRITE(3,*) 'plot(xr2,yr2,''--'');'
          WRITE(3,*) 'xlabel(''X'');'
          WRITE(3,*) 'ylabel(''Y'');'
          WRITE(3,*) 'zlabel(''Pressure, p^*'');'
          WRITE(3,*)
        ELSE IF (M .EQ. NP) THEN
          WRITE(3,*) 'shading interp;'
          WRITE(3,*) 'hold off;'
          WRITE(3,*) 'MOV2(', N ,') = getframe;'
          WRITE(3,*)
        END IF
   40   FORMAT(A4,F6.4,A12)
   45   FORMAT(A9,2I1,A5,2I1,A5,2I1,A2)
   50   FORMAT(A9,I2,I1,A5,I2,I1,A5,I2,I1,A2)
   55   FORMAT(A9,I3,I1,A5,I3,I1,A5,I3,I1,A2)

   60 CONTINUE

      RETURN    
      END

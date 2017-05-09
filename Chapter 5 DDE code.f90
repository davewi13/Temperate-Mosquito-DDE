MODULE define_DDEs

  IMPLICIT NONE
  
  ! Set number of equations, delays and event functions
  INTEGER, PARAMETER :: NEQN=21,NLAGS=7,NEF=1
  INTEGER :: K
  DOUBLE PRECISION :: PI=3.1415927D0
  ! Set parameter values for predation and egg raft size
  DOUBLE PRECISION :: densa,densr,densh,UPS,SHP,VOL,MAXEGG
  ! Set WNV transmission parameters
  DOUBLE PRECISION :: DEATHbird=0.000685D0,DEATHbirdWNV=0.167D0
  DOUBLE PRECISION :: RECOVERY=0.25D0,PHH=0.33D0,VERTICAL=0.004D0
  DOUBLE PRECISION :: tranMB=0.88D0,tranBM=0.4D0,CAR=400D0,INOCB=2D0
  DOUBLE PRECISION :: M=9.8D0,A=6.4D0,PHASE=28.9D0,DTR=9.4D0,WARM=0D0,INOCT=151D0
  
CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

  	! Parameters for DDE solver
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY

    ! Define variable names
    DOUBLE PRECISION :: TEMPnowair,TEMPGCair,TEMPEIPair,TEMPEair,TEMPLair,TEMPPair,TEMPELair,TEMPLPair,TEMPELPair
    DOUBLE PRECISION :: TEMPnowwater,TEMPEwater,TEMPLwater,TEMPPwater,TEMPELwater,TEMPLPwater,TEMPELPwater
    DOUBLE PRECISION :: BIRTHnow,BIRTHE,BIRTHEL,BIRTHELP,BIRTHbird
    DOUBLE PRECISION :: EGGMATnow,EGGMATE,EGGMATL,EGGMATEL,EGGMATLP,EGGMATELP
    DOUBLE PRECISION :: LARMATnow,LARMATL,LARMATP,LARMATLP
    DOUBLE PRECISION :: REt,RLt,RPt,RSAt,REAt,RIAt,MEt,MLt,MPt,MSAt,MEAt,DEt,DLt,DPt,DSAt,DEAt,DIAt
    DOUBLE PRECISION :: dEdt,dLdt,dPdt,dSAdt,dEAdt,dIAdt,dSEdt,dSLdt,dSPdt,dDEdt,dDLdt,dDELdt,dDPdt,dDLPdt,dDELPdt,dGCdt
    DOUBLE PRECISION :: dSBdt,dIBdt,dRBdt,dEIPdt,dSEAdt
    DOUBLE PRECISION :: E,LAR,PUP,SA,EA,IA,SE,SL,SP,DE,DL,DP,GC,EIP,SEA,NB,SB,IB,RB,NBEIP
    DOUBLE PRECISION :: INOCM,BITINGnow,BITINGEIP
    DOUBLE PRECISION :: DEATHeggnow,DEATHeggE,DEATHlarnow,DEATHlarL,DEATHpupnow,DEATHpupP,DEATHadunow,DEATHaduEIP
    DOUBLE PRECISION :: PUPMATnow,PUPMATP
    DOUBLE PRECISION :: DIAPAUSEE,DIAPAUSEEL,DIAPAUSEELP,DIAPAUSEnow,DIAPAUSEEIP
    DOUBLE PRECISION :: PPnow,PPE,PPEL,PPELP,PPEIP
    DOUBLE PRECISION :: GCnow,GCC,EIPnow,EIPdelay

    ! Give names to responses of system of equations    
	E = Y(1)
	LAR = Y(2)
    PUP = Y(3)
	SA = Y(4)
    EA = Y(5)
    IA = Y(6)
    SE = Y(7)
    SL = Y(8)
    SP = Y(9)
    SEA = Y(10)
    DE = Y(11)
    DL = Y(12)
    DP = Y(13)
    GC = Y(17)
    EIP = Y(18)
    SB = Y(19)
    IB = Y(20)
    RB = Y(21)
    NB = (SB+IB+RB)
    NBEIP = Z(19,7)+Z(20,7)+Z(21,7)

    ! Extract temperature values at various time points
	TEMPnowair = TEMPAIR(T)
    TEMPnowwater = TEMPWATER(T,TEMPnowair)
    
    TEMPEair = TEMPAIR(T-DE)
    TEMPEwater = TEMPWATER(T-DE,TEMPEair)

    TEMPLair = TEMPAIR(T-DL)
    TEMPLwater = TEMPWATER(T-DL,TEMPLair)

    TEMPPair = TEMPAIR(T-DP)
    TEMPPwater = TEMPWATER(T-DP,TEMPPair)

    TEMPELair = TEMPAIR(T-DL-Z(11,4))
    TEMPELwater = TEMPWATER(T-DL-Z(11,4),TEMPELair)

    TEMPELPair = TEMPAIR(T-DP-Z(12,6)-Z(11,5))
    TEMPELPwater = TEMPWATER(T-DP-Z(12,6)-Z(11,5),TEMPELPair)

    TEMPLPair = TEMPAIR(T-DP-Z(12,6))
    TEMPLPwater = TEMPWATER(T-DP-Z(12,6),TEMPLPair)
    
    TEMPGCair = TEMPAIR(T-GC)
    TEMPEIPair = TEMPAIR(T-EIP)

	! Calculate photoperiod values for a range of time points
    PPnow = DAYLIGHT(T)
    PPE = DAYLIGHT(T-DE)
    PPEL = DAYLIGHT(T-DL-Z(11,4))
    PPELP = DAYLIGHT(T-DP-Z(12,6)-Z(11,5))
    PPEIP = DAYLIGHT(T-EIP)

	! Calculate gonotrophic cycle length for various time points
    GCnow = GONOTROPHIC(TEMPnowair)
    GCC = GONOTROPHIC(TEMPGCair)

    EIPnow = EXTRINSIC_INCUBATION(TEMPnowair)
    EIPdelay = EXTRINSIC_INCUBATION(TEMPEIPair)

	! Calculate diapause percentage for given photoperiod values.
    ! Use these diapause percentages to caluclate birth rates.
    ! Do this for both spring and autumn diapause thresholds.
    IF (DAYLIGHT(T) > DAYLIGHT(T-1)) THEN
      DIAPAUSEnow = DIAPAUSE_SPRING(PPnow)
      DIAPAUSEE = DIAPAUSE_SPRING(PPE)
      DIAPAUSEEL = DIAPAUSE_SPRING(PPEL)
      DIAPAUSEELP = DIAPAUSE_SPRING(PPELP)
      DIAPAUSEEIP = DIAPAUSE_SPRING(PPEIP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GC)
      BIRTHE = BIRTH(DIAPAUSEE,Z(17,1))
      BIRTHEL = BIRTH(DIAPAUSEEL,Z(17,2))
      BIRTHELP = BIRTH(DIAPAUSEELP,Z(17,3))
    ELSE
      DIAPAUSEnow = DIAPAUSE_AUTUMN(PPnow)
      DIAPAUSEE = DIAPAUSE_AUTUMN(PPE)
      DIAPAUSEEL = DIAPAUSE_AUTUMN(PPEL)
      DIAPAUSEELP = DIAPAUSE_AUTUMN(PPELP)
      DIAPAUSEEIP = DIAPAUSE_AUTUMN(PPEIP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GC)
      BIRTHE = BIRTH(DIAPAUSEE,Z(17,1))
      BIRTHEL = BIRTH(DIAPAUSEEL,Z(17,2))
      BIRTHELP = BIRTH(DIAPAUSEELP,Z(17,3))
    END IF

	! Calculate death rates for each life stage
	DEATHeggnow = DEATHegg(TEMPnowwater)
    DEATHeggE = DEATHegg(TEMPEwater)

    DEATHlarnow = DEATHlar(TEMPnowwater)
    DEATHlarL = DEATHlar(TEMPLwater)

    DEATHpupnow = DEATHpup(TEMPnowwater)
    DEATHpupP = DEATHpup(TEMPPwater)
    
    DEATHadunow = DEATHadu(TEMPnowair,GC,T)
    DEATHaduEIP = DEATHadu(TEMPEIPair,Z(17,7),T-EIP)

    BIRTHbird = BIRD_BIRTH_FUNC(T)
    
	! Calculate development rates for the immature stages
    LARMATnow = LARMATURATION(TEMPnowwater)
    LARMATL = LARMATURATION(TEMPLwater)
    LARMATP = LARMATURATION(TEMPPwater)
    LARMATLP = LARMATURATION(TEMPLPwater)

    EGGMATnow = EGGMATURATION(TEMPnowwater)
    EGGMATE = EGGMATURATION(TEMPEwater)
    EGGMATL = EGGMATURATION(TEMPLwater)
    EGGMATEL = EGGMATURATION(TEMPELwater)
	EGGMATLP = EGGMATURATION(TEMPLPwater)
	EGGMATELP = EGGMATURATION(TEMPELPwater)

    PUPMATnow = PUPMATURATION(TEMPnowwater)
    PUPMATP = PUPMATURATION(TEMPPwater)

    BITINGnow = 0.5D0*DIAPAUSEnow*GONOTROPHIC(TEMPnowair)
    BITINGEIP = 0.5D0*DIAPAUSEEIP*GONOTROPHIC(TEMPEIPair)

	! Innoculate the system with a given number of adults.
    INOCM = INOCCULATEM(T)

	! Delay differential equations for immature stage durations
    dDEdt = 1D0 - EGGMATnow/EGGMATE
    dDLdt = 1D0 - LARMATnow/LARMATL
    dDPdt = 1D0 - PUPMATnow/PUPMATP

    ! Delay differential equation for gonotrophic cycle duration 
    dGCdt = 1D0 - GCnow/GCC
    dEIPdt = 1D0 - EIPnow/EIPdelay

    ! Delay differential equations for stage durations referenced back through previous stages
    dDELdt = (1D0 - dDLdt) * (1D0 - EGGMATL/EGGMATEL)
    dDLPdt = (1D0 - dDPdt) * (1D0 - LARMATP/LARMATLP)
    dDELPdt = (1D0 - dDPdt - dDLPdt) * (1D0 - EGGMATLP/EGGMATELP)
          
	! Recruitment equations for each stage
	REt = BIRTHnow * (SA + EA + IA)
	RLt = BIRTHE * (Z(4,1) +Z(5,1) +Z(6,1)) * SE * EGGMATnow/EGGMATE
    RPt = BIRTHEL * (Z(4,2) +Z(5,2) +Z(6,2)) * Z(7,4) * SL * LARMATnow/LARMATL * EGGMATL/EGGMATEL
    
    RSAt = BIRTHELP * (Z(4,3) +Z(5,3) +Z(6,3) * (1D0-VERTICAL)) * Z(7,5) * Z(8,6) * SP * PUPMATnow/PUPMATP *&
    LARMATP/LARMATLP * EGGMATLP/EGGMATELP + INOCM
    
    REAt = BITINGnow * SA * (tranBM * IB) / NB
    
    RIAt = BITINGEIP * Z(4,7) * SEA * EIPnow/EIPdelay * (tranBM * Z(20,7)) /&
    NBEIP + BIRTHELP * Z(6,3) * VERTICAL * Z(7,5) * Z(8,6) * SP * PUPMATnow/PUPMATP *&
    LARMATP/LARMATLP * EGGMATLP/EGGMATELP

    
	! Maturation equations for each stage
    MEt = RLt
    MLt = RPt
    MPt = BIRTHELP * Z(4,3) * Z(7,5) * Z(8,6) * SP * PUPMATnow/PUPMATP * LARMATP/LARMATLP * EGGMATLP/EGGMATELP
    MSAt = REAt
    MEAt = BITINGEIP * Z(4,7) * SEA * EIPnow/EIPdelay * (tranBM * Z(20,7)) / NBEIP

	! Death equations for each stage
	DEt = DEATHeggnow * E
    DLt = (densa*densr*(((1D0+COS(2D0*PI*(T-182.5D0-UPS)/365D0))/2D0)**SHP)*LAR/(VOL+densh*LAR)+DEATHlarnow)*LAR
    DPt = DEATHpupnow * PUP
    DSAt = DEATHadunow * SA
    DEAt = DEATHadunow * EA
    DIAt = DEATHadunow * IA

    ! Balance equations for each life stage
	dEdt = REt - MEt - DEt
	dLdt = RLt - MLt - DLt
    dPdt = RPt - MPt - DPt
    
	dSAdt = RSAt - MSAt - DSAt
    dEAdt = REAt - MEAt - DEAt
    dIAdt = RIAt - DIAt
    
    dSBdt = BIRTHbird * NB - BITINGnow * tranMB * IA * SB / NB - PHH * IB * SB / NB - (DEATHbird + BIRTHbird*NB/CAR) * SB
    
    dIBdt = BITINGnow*tranMB*IA*SB/NB + PHH*IB*SB/NB - (DEATHbird+DEATHbirdWNV+BIRTHbird*NB/CAR)*IB - RECOVERY*IB
    
    dRBdt = RECOVERY * IB - (DEATHbird + BIRTHbird*NB/CAR) * RB

	! Survival equations for the immature stages
    dSEdt = SE * ((EGGMATnow * DEATHeggE / EGGMATE) - DEATHeggnow)
	dSLdt =SL*(((densa*densr*(((1D0+COS(2D0*PI*(T-DL-182.5D0-UPS)/365D0))/2D0)**SHP)*Z(2,4)/(VOL+densh*Z(2,4)))+DEATHlarL)&
    * (1-dDLdt) - (densa*densr*(((1D0+COS(2D0*PI*(T-182.5D0-UPS)/365D0))/2D0)**SHP)*LAR / (VOL+densh*LAR)) - DEATHlarnow)
    dSPdt = SP * ((PUPMATnow * DEATHpupP / PUPMATP) - DEATHpupnow)
    dSEAdt = SEA * ((EIPnow * DEATHaduEIP / EIPdelay) - DEATHadunow)

    ! Derivatives for the integrator:
    DY = (/ dEdt, dLdt, dPdt, dSAdt, dEAdt, dIAdt, dSEdt, dSLdt, dSPdt, dSEAdt, dDEdt, dDLdt, dDPdt, dDELdt, dDLPdt, dDELPdt,&
    dGCdt, dEIPdt, dSBdt, dIBdt, dRBdt /)

    RETURN
    END SUBROUTINE DDES

  SUBROUTINE BETA(T,Y,BVAL)
    
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION, DIMENSION(NLAGS) :: BVAL
    INTENT(IN)  :: T,Y
    INTENT(OUT) :: BVAL

    ! Set the delay values
	! T - Eggdelay(T)
	BVAL(1) = T-Y(11)
    ! T - Lardelay(T) - Eggdelay(T-Lardelay(T))				
    BVAL(2) = T-Y(12)-Y(14)
    !T - Pupdelay(T) - Lardelay(T-Pupdelay(T)) - Eggdelay(T-Pupdelay(T)-Lardelay(T-Pupdelay(T)))
    BVAL(3) = T-Y(13)-Y(15)-Y(16)	
    !T - Lardelay(T)
    BVAL(4) = T-Y(12)				
    !T - Pupdelay(T) - Lardelay(T-Pupdelay(T))
    BVAL(5) = T-Y(13)-Y(15)			
    !T - Pupdelay(T)
    BVAL(6) = T-Y(13)				
    !T - EIPdelay(T)
    BVAL(7) = T-Y(18)
   	  
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE HISTORY(T,Y)
    DOUBLE PRECISION :: T,TEMPhist,TEMPhistair
    DOUBLE PRECISION, DIMENSION(2*NEQN) :: Y
    INTENT(IN)  :: T
    INTENT(OUT) :: Y

	! Set the temperatures for T < 0
    TEMPhistair = TEMPAIR(T)
    TEMPhist = TEMPWATER(T,TEMPhistair)
    
	! Set historical values for all stages to be zero
    Y(1) = 0D0
    Y(2) = 0D0
    Y(3) = 0D0
    Y(4) = 0D0
    Y(5) = 0D0
    Y(6) = 0D0

    ! Calculate historical survival rates based on temperature
    Y(7) = EXP(-DEATHegg(TEMPhist)*(1D0/EGGMATURATION(TEMPhist)))
    Y(8) = EXP(-DEATHlar(TEMPhist)*(1D0/LARMATURATION(TEMPhist)))
    Y(9) = EXP(-DEATHpup(TEMPhist)*(1D0/PUPMATURATION(TEMPhist)))

	! Calculate historical development rates based on temperature
    Y(11) = 1D0/EGGMATURATION(TEMPhist)
    Y(12) = 1D0/LARMATURATION(TEMPhist)
    Y(13) = 1D0/PUPMATURATION(TEMPhist)
    Y(14) = 1D0/EGGMATURATION(TEMPhist)
    Y(15) = 1D0/LARMATURATION(TEMPhist)
    Y(16) = 1D0/EGGMATURATION(TEMPhist)
    Y(17) = 1D0/GONOTROPHIC(TEMPhistair)

    Y(10) = EXP(-DEATHadu(TEMPhistair,Y(17),T)*(1D0/EXTRINSIC_INCUBATION(TEMPhistair)))

    ! Set historical EIP
    Y(18) = 1D0/EXTRINSIC_INCUBATION(TEMPhistair)
    
	! Set historical values for birds to zero
	Y(19) = 0.875D0*CAR     
    Y(20) = 0D0
    Y(21) = 0D0
    
    RETURN
  END SUBROUTINE HISTORY


  SUBROUTINE EF(T,Y,DY,Z,G)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    DOUBLE PRECISION, DIMENSION(NEF) :: G
    INTENT(IN) :: T,Y,DY,Z
    INTENT(OUT) :: G

    ! Set events as turning points in adult time series and diapause entry/exit
    G(1) = T-INOCT

     RETURN
  END SUBROUTINE EF

  SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
                  DIRECTION,ISTERMINAL,QUIT)
  ! Function to change a flag so that the DDE model will
  ! be evaluated in subroutine DDES instead of the ODE model.

     INTEGER :: NEVENT
     INTEGER, DIMENSION(NEF) :: DIRECTION
     DOUBLE PRECISION :: TEVENT,HINIT
     DOUBLE PRECISION, DIMENSION(NEQN) :: YEVENT,DYEVENT
     LOGICAL :: QUIT
     LOGICAL, DIMENSION(NEF) :: ISTERMINAL
     INTENT(IN) :: NEVENT,TEVENT
     INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,&
                      ISTERMINAL,QUIT

     YEVENT(20) = 0.0025D0*CAR

    RETURN
  END SUBROUTINE CHNG

  DOUBLE PRECISION FUNCTION TEMPAIR(T)
     
    DOUBLE PRECISION :: T
    
	! Calculate air temperature for time, T
    TEMPAIR = M + A * COS(2D0 * Pi * (T-PHASE-182.5) / 365D0) - DTR/2 * COS(2D0* Pi * T) + WARM
    IF (T<0D0) THEN
      TEMPAIR = M + A * COS(2D0 * Pi * (0D0-PHASE-182.5) / 365D0) - DTR/2 * COS(2D0* Pi * 0D0) + WARM
    END IF
    
  RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION TEMPWATER(T,AIRTEMP)
     
    DOUBLE PRECISION :: T,AIRTEMP
    
	! Calculate water temperature for time, T
    TEMPWATER = 0.9505D0*AIRTEMP + 3.8887D0
    
  RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION INOCCULATEM(T)

    DOUBLE PRECISION :: T

	! Set inoculation value
    IF (T < 1D0 .AND. T > 0D0) THEN
	   INOCCULATEM = 5000D0
    ELSE
       INOCCULATEM = 0D0
    END IF

	RETURN
  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DAYLIGHT(T)
     
    DOUBLE PRECISION :: T,EPS,NUM,DEN
    REAL, PARAMETER :: Pi = 3.1415927D0, L = 51D0

    ! Calculate photperiod value
    EPS = ASIN(0.39795D0 * COS(0.2163108D0 + 2 * ATAN(0.9671396D0 * TAN(0.00860D0 * (T-185.5D0)))))
    NUM = SIN(0.8333D0*Pi/180D0) + (SIN(L*Pi/180D0) * SIN(EPS))
    DEN = COS(L*Pi/180D0) * COS(EPS)
    DAYLIGHT = 24D0 - (24D0/Pi) * ACOS(NUM / DEN)
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DIAPAUSE_SPRING(PP)
     
    DOUBLE PRECISION :: PP

    ! Set spring photoperiod threshold
    DIAPAUSE_SPRING = 1D0 / (1D0 + EXP(5D0*(13.7D0-PP)))
    RETURN
    
  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DIAPAUSE_AUTUMN(PP)
     
    DOUBLE PRECISION :: PP

    ! Set autumn photoperiod threshold
    DIAPAUSE_AUTUMN = 1D0 / (1D0 + EXP(3.5D0*(15D0-PP)))
    RETURN
    
  END FUNCTION

  DOUBLE PRECISION FUNCTION BIRTH(DIAPAUSE,GONOTROPHICtime)
     
    DOUBLE PRECISION :: GONOTROPHICtime,EGGRAFT,DIAPAUSE

	! Set birth rate for mosquitoes
    EGGRAFT = DIAPAUSE*MAXEGG*0.5D0
    BIRTH = EGGRAFT/GONOTROPHICtime
  
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION BIRD_BIRTH_FUNC(T)
     
    DOUBLE PRECISION :: T
    DOUBLE PRECISION :: K_BIRD=0.15D0,PHI_BIRD=50D0,S_BIRD=10D0

	! Set birth rate for birds
    BIRD_BIRTH_FUNC = 0.5D0 * K_BIRD * EXP(-S_BIRD * (COS((PI*T+PHI_BIRD)/365D0))**2)
  
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION GONOTROPHIC(TEMP)
     
	DOUBLE PRECISION :: TEMP
	DOUBLE PRECISION :: KG=0.2024D0,QG=74.48D0,BG=0.2456D0

	! Set gonotrophic cycle rate
	IF (TEMP < 0D0) THEN
      GONOTROPHIC = 0.0333D0
    ELSE
	  GONOTROPHIC = KG / (1+QG*EXP(-BG*TEMP))
    END IF
	IF(GONOTROPHIC < 0.0333D0) THEN
      GONOTROPHIC = 0.0333D0
    END IF
    
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION EXTRINSIC_INCUBATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    
	! Set EIP development rate
	EXTRINSIC_INCUBATION = 0.0092D0 * TEMP - 0.132D0
    IF (EXTRINSIC_INCUBATION .LE. 0.005D0) THEN
      EXTRINSIC_INCUBATION = 0.005D0
    END IF

    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHegg(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	! Calculate egg death rate
	DEATHegg = U3 * EXP(((TEMP-U4)/U5)**2)
	IF (DEATHegg > 1D0) THEN
      DEATHegg = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHlar(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	! Calculate larval death rate
	DEATHlar = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHlar > 1D0) THEN
      DEATHlar = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHpup(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	! Calculate pupal death rate
	DEATHpup = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHpup > 1D0) THEN
      DEATHpup = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHadu(TEMP,GONOTROPHICtime,T)
    
	DOUBLE PRECISION :: TEMP,GONOTROPHICtime,T
    DOUBLE PRECISION :: ALPHA=2.166D-8,BETA=4.483D0,PI=3.1415927D0,MULTIPLIER=8D0,SIGMASQ=4D0
	
	! Calculate adult death rate
	IF (TEMP < 0D0) THEN
    	DEATHadu = 0.003D0
    ELSE
		DEATHadu = ALPHA*(TEMP**BETA)
    END IF
	IF (DEATHadu < 0.003D0) THEN
      DEATHadu = 0.003D0
    END IF
   	DEATHadu = DEATHadu + (MULTIPLIER/SQRT(SIGMASQ*2D0*PI))*EXP((-1D0/(SIGMASQ*2D0))*&
    (MOD(T,365D0)-GONOTROPHICtime-109D0)**2)
   
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION EGGMATURATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.0022D0,BETA=1.77D0
	
	! Calculate egg development rate
	IF (TEMP < 0D0) THEN
    	EGGMATURATION = 0.016667D0
    ELSE
		EGGMATURATION = ALPHA*(TEMP**BETA)
	END IF
	IF (EGGMATURATION < 0.016667D0) THEN
    	EGGMATURATION = 0.016667D0
    END IF     

    RETURN

  END FUNCTION

  DOUBLE PRECISION FUNCTION LARMATURATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.00315D0,BETA=1.12D0
	
	! Calculate larval development rate
	IF (TEMP < 0D0) THEN
    	LARMATURATION = 0.016667D0
    ELSE
		LARMATURATION = ALPHA*(TEMP**BETA)
	END IF
	IF (LARMATURATION < 0.016667D0) THEN
    	LARMATURATION = 0.016667D0
    END IF    

    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION PUPMATURATION(TEMP)
     
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.0007109D0,BETA=1.8865648D0
	
	! Calculate pupal maturation rate
	IF (TEMP < 0D0) THEN
    	PUPMATURATION = 0.016667D0
    ELSE
		PUPMATURATION = ALPHA*(TEMP**BETA)
	END IF
    IF (PUPMATURATION < 0.016667D0) THEN
    	PUPMATURATION = 0.016667D0
    END IF 
    
    RETURN
    
  END FUNCTION

END MODULE define_DDEs

!******************************************************************

PROGRAM basicmodel

  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  INTEGER :: I,J ! Local variables

  INTEGER, DIMENSION(3) :: NVAR = (/NEQN,NLAGS,NEF/)

  ! Set length of solution and output points
  INTEGER, PARAMETER :: NOUT=366D0
  DOUBLE PRECISION, PARAMETER :: T0=0D0,TFINAL=365D0
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN= &
  (/ (T0+(I-1)*((TFINAL - T0)/(NOUT-1)), I=1,NOUT) /)

  TYPE(DDE_SOL) :: SOL
  TYPE(DDE_OPTS) :: OPTS
  
  ! Set length of maximum delay
  DOUBLE PRECISION :: MAXDELAY = 200D0
    
  ! Set options for DDE solver
  OPTS = DDE_SET(RE=1D-5,AE=1D-5,MAX_STEPS=100000000,MAX_DELAY=MAXDELAY,TRIM_FREQUENCY=10000)  

  ! Run DDE solver code
  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,TSPAN,EVENT_FCN=EF,CHANGE_FCN=CHNG,OPTIONS=OPTS)
  
  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN
     ! Output results
	 OPEN(UNIT=10, FILE='WallingfordCosine.dat')
     DO I = 1,SOL%NPTS
       WRITE(UNIT=10,FMT='(28E14.5E3)') INOCT, WARM, SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
     END DO
     CLOSE(10)

  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  END IF


  STOP
END PROGRAM basicmodel
!!!!!!!!!!!!!
!!!!!!!!!!!!!
! NOTE: The file "Chapter 2 DDE code.f90" contains a more up-to-date version of this code and should be used for any future analyses!!!

MODULE define_DDEs

  IMPLICIT NONE

  ! Set numbers of equations and parameter values
  INTEGER, PARAMETER :: NEQN=13,NLAGS=6,NEF=6
  DOUBLE PRECISION :: MAXEGG=100D0
  DOUBLE PRECISION :: M=10D0,PHASE=30D0,POWER=1D0,A=6D0,PI=3.1415927D0
  DOUBLE PRECISION :: B1=0.5D0,B2=100000D0
  
CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY

    ! Local variables
    DOUBLE PRECISION :: TEMPnow,TEMPE,TEMPL,TEMPP,TEMPEL,TEMPLP,TEMPELP
	DOUBLE PRECISION :: BIRTHnow,BIRTHE,BIRTHEL,BIRTHELP
    DOUBLE PRECISION :: EGGMATnow,EGGMATE,EGGMATL,EGGMATEL,EGGMATLP,EGGMATELP
	DOUBLE PRECISION :: LARMATnow,LARMATL,LARMATP,LARMATLP
	DOUBLE PRECISION :: PUPMATnow,PUPMATP
	DOUBLE PRECISION :: PPnow,PPE,PPEL,PPELP
	DOUBLE PRECISION :: REt,RLt,RPt,RAt
	DOUBLE PRECISION :: MEt,MLt,MPt
	DOUBLE PRECISION :: DEt,DLt,DPt,DAt
	DOUBLE PRECISION :: dEdt,dLdt,dPdt,dAdt,dSEdt,dSLdt,dSPdt,dDEdt,dDLdt,dDELdt,dDPdt,dDLPdt,dDELPdt
	DOUBLE PRECISION :: E,LAR,PUP,ADU,SE,SL,SP,DE,DL,DP,DLP,DELP
	DOUBLE PRECISION :: INOC
	DOUBLE PRECISION :: DEATHeggnow,DEATHeggE,DEATHlarnow,DEATHlarL,DEATHpupnow,DEATHpupP,DEATHadunow
	DOUBLE PRECISION :: GONOTROPHICnow,GONOTROPHICE,GONOTROPHICEL,GONOTROPHICELP
	DOUBLE PRECISION :: DIAPAUSEnow,DIAPAUSEE,DIAPAUSEEL,DIAPAUSEELP

    ! Set number of individuals for initial innoculation
    INOC = INOCCULATE(T)
	
    ! Local solution variables and delayed solution values:
	E = Y(1)
	LAR = Y(2)
    PUP = Y(3)
	ADU = Y(4)
    SE = Y(5)
    SL = Y(6)
    SP = Y(7)
    DE = Y(8)
    DL = Y(9)
    DP = Y(10)
    DLP = Y(12)
    DELP = Y(13)

	! Define temperature at a range of time points
    TEMPnow = TEMP(T)
    TEMPE = TEMP(T-DE)
    TEMPL = TEMP(T-DL)
    TEMPP = TEMP(T-DP)
    TEMPEL = TEMP(T-DL-Z(8,4))
    TEMPELP = TEMP(T-DP-Z(9,6)-Z(8,5))
    TEMPLP = TEMP(T-DP-Z(9,6))
	
	! Define photoperiod at a range of time points
	PPnow = DAYLIGHT(T)
    PPE = DAYLIGHT(T-DE)
    PPEL = DAYLIGHT(T-DL-Z(8,4))
    PPELP = DAYLIGHT(T-DP-Z(9,6)-Z(8,5))

	! Define gonotrophic cycle length at a range of time points
    GONOTROPHICnow = GONOTROPHIC(TEMPnow)
    GONOTROPHICE = GONOTROPHIC(TEMPE)
    GONOTROPHICEL = GONOTROPHIC(TEMPEL)
    GONOTROPHICELP = GONOTROPHIC(TEMPELP)
	
	! Define diapause proportions and birth rates for spring and then autumn
	IF (DAYLIGHT(T) > DAYLIGHT(T-1)) THEN
      DIAPAUSEnow = DIAPAUSE_SPRING(PPnow)
      DIAPAUSEE = DIAPAUSE_SPRING(PPE)
      DIAPAUSEEL = DIAPAUSE_SPRING(PPEL)
      DIAPAUSEELP = DIAPAUSE_SPRING(PPELP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GONOTROPHICnow,ADU,T)
      BIRTHE = BIRTH(DIAPAUSEE,GONOTROPHICE,ADU,T)
      BIRTHEL = BIRTH(DIAPAUSEEL,GONOTROPHICEL,ADU,T)
      BIRTHELP = BIRTH(DIAPAUSEELP,GONOTROPHICELP,ADU,T)
    ELSE
      DIAPAUSEnow = DIAPAUSE_AUTUMN(PPnow)
      DIAPAUSEE = DIAPAUSE_AUTUMN(PPE)
      DIAPAUSEEL = DIAPAUSE_AUTUMN(PPEL)
      DIAPAUSEELP = DIAPAUSE_AUTUMN(PPELP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GONOTROPHICnow,ADU,T)
      BIRTHE = BIRTH(DIAPAUSEE,GONOTROPHICE,ADU,T)
      BIRTHEL = BIRTH(DIAPAUSEEL,GONOTROPHICEL,ADU,T)
      BIRTHELP = BIRTH(DIAPAUSEELP,GONOTROPHICELP,ADU,T)
    END IF

	! Set death rates at a range of time points for each life stage
	DEATHeggnow = DEATHegg(TEMPnow)
    DEATHeggE = DEATHegg(TEMPE)

    DEATHlarnow = DEATHlar(TEMPnow)
    DEATHlarL = DEATHlar(TEMPL)

    DEATHpupnow = DEATHpup(TEMPnow)
    DEATHpupP = DEATHpup(TEMPP)

    DEATHadunow = DEATHadu(TEMPnow)

	! Set larval development rates at a range of time points
    LARMATnow = LARMATURATION(TEMPnow)
    LARMATL = LARMATURATION(TEMPL)
    LARMATP = LARMATURATION(TEMPP)
    LARMATLP = LARMATURATION(TEMPLP)

	! Set egg development rates at a range of time points
    EGGMATnow = EGGMATURATION(TEMPnow)
    EGGMATE = EGGMATURATION(TEMPE)
    EGGMATL = EGGMATURATION(TEMPL)
    EGGMATEL = EGGMATURATION(TEMPEL)
	EGGMATLP = EGGMATURATION(TEMPLP)
	EGGMATELP = EGGMATURATION(TEMPELP)

	! Set pupal development rates at a range of time points
    PUPMATnow = PUPMATURATION(TEMPnow)
    PUPMATP = PUPMATURATION(TEMPP)

	! Define stage duration (delay) equations
    dDEdt = 1D0 - EGGMATnow/EGGMATE
    dDLdt = 1D0 - LARMATnow/LARMATL
    dDPdt = 1D0 - PUPMATnow/PUPMATP
    dDELdt = (1D0 - dDLdt) * (1D0 - EGGMATL/EGGMATEL)
    dDLPdt = (1D0 - dDPdt) * (1D0 - LARMATP/LARMATLP)
    dDELPdt = (1D0 - dDPdt - dDLPdt) * (1D0 - EGGMATLP/EGGMATELP)
          
	! Define recruitment for each stage
	REt = BIRTHnow * ADU
	RLt = BIRTHE * Z(4,1) * SE * EGGMATnow/EGGMATE
    RPt = BIRTHEL * Z(4,2) * Z(5,4) * SL * LARMATnow/LARMATL * EGGMATL/EGGMATEL
    RAt = BIRTHELP * Z(4,3) * Z(5,5) * Z(6,6) * SP * PUPMATnow/PUPMATP * LARMATP/LARMATLP * EGGMATLP/EGGMATELP + INOC
    
	! Define maturation for each stage
    MEt = RLt
    MLt = RPt
    MPt = BIRTHELP * Z(4,3) * Z(5,5) * Z(6,6) * SP * PUPMATnow/PUPMATP * LARMATP/LARMATLP * EGGMATLP/EGGMATELP	

	! Define death for each stage
	DEt = DEATHeggnow * E
    DLt = (B1*LAR/(B2+LAR) + DEATHlarnow) * LAR
    DPt = DEATHpupnow * PUP
    DAt = DEATHadunow * ADU
        
    !Local derivatives
    dEdt = REt - MEt - DEt
	dLdt = RLt - MLt - DLt
    dPdt = RPt - MPt - DPt
	dAdt = RAt - DAt
    
    dSEdt = SE * ((EGGMATnow * DEATHeggE / EGGMATE) - DEATHeggnow)
    dSLdt = SL * (((B1*Z(2,4) / (B2+Z(2,4))) + DEATHlarL) * (1-dDLdt) - (B1*LAR / (B2+LAR)) - DEATHlarnow)
    dSPdt = SP * ((PUPMATnow * DEATHpupP / PUPMATP) - DEATHpupnow)
    
    
    ! Derivatives for the integrator:
    DY = (/ dEdt, dLdt, dPdt, dAdt, dSEdt, dSLdt, dSPdt, dDEdt, dDLdt, dDPdt, dDELdt, dDLPdt, dDELPdt /)
    

    RETURN
    END SUBROUTINE DDES

  SUBROUTINE BETA(T,Y,BVAL)
    
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION, DIMENSION(NLAGS) :: BVAL
    INTENT(IN)  :: T,Y
    INTENT(OUT) :: BVAL
	  
	! Set the delays for the solver 
	BVAL(1) = T-Y(8)				!T - Eggdelay(T)
    BVAL(2) = T-Y(9)-Y(11)			!T - Lardelay(T) - Eggdelay(T-Lardelay(T))
    BVAL(3) = T-Y(10)-Y(12)-Y(13)	!T - Pupdelay(T) - Lardelay(T-Pupdelay(T)) - Eggdelay(T-Pupdelay(T)-Lardelay(T-Pupdelay(T)))
    BVAL(4) = T-Y(9)				!T - Lardelay(T)
    BVAL(5) = T-Y(10)-Y(12)			!T - Pupdelay(T) - Lardelay(T-Pupdelay(T))
    BVAL(6) = T-Y(10)				!T - Pupdelay(T)
   	  
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE HISTORY(T,Y)
    DOUBLE PRECISION :: T,TEMPhist,TEMPhistL,TEMPhistP,TEMPhistLP
    DOUBLE PRECISION, DIMENSION(2*NEQN) :: Y
    INTENT(IN)  :: T
    INTENT(OUT) :: Y
      
	! Set historical values for time < 0.
	! Note that whilst this will reproduce the results in the paper further work
	! has highlighted that temperatures should be fixed for T<0 -though the difference is small in this case.  Please use the file
	! "Chapter 2 code.f90" instead.  See Chapter 2 of the thesis (http://theses.gla.ac.uk/8450/1/2017ewingphd.pdf).
	TEMPhist = TEMP(T)
      
    Y(1) = 0D0
    Y(2) = 0D0
    Y(3) = 0D0
    Y(4) = 0D0

    Y(8) = 1D0/EGGMATURATION(TEMPhist)
    Y(9) = 1D0/LARMATURATION(TEMPhist)
    Y(10) = 1D0/PUPMATURATION(TEMPhist)
     
    Y(5) = EXP(-DEATHegg(TEMPhist)*(1D0/EGGMATURATION(TEMPhist)))
    Y(6) = EXP(-DEATHlar(TEMPhist)*(1D0/LARMATURATION(TEMPhist)))
    Y(7) = EXP(-DEATHpup(TEMPhist)*(1D0/PUPMATURATION(TEMPhist)))
      
    TEMPhistL = TEMP(T-Y(9))
    TEMPhistP = TEMP(T-Y(10))
      
    Y(11) = 1D0/EGGMATURATION(TEMPhistL)
    Y(12) = 1D0/LARMATURATION(TEMPhistP)
      
    TEMPhistLP = TEMP(T-Y(10)-Y(12))
    
    Y(13) = 1D0/EGGMATURATION(TEMPhistLP)

    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE EF(T,Y,DY,Z,G)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    DOUBLE PRECISION, DIMENSION(NEF) :: G
    INTENT(IN) :: T,Y,DY,Z
    INTENT(OUT) :: G

	! Event function - only required if searching for turning points or other events. See DDE_solver documentation.
    G(1) = DY(4)
	G(2) = 14D0 - DAYLIGHT(T)
	G(3) = 13D0 - DAYLIGHT(T)
	G(4) = 10D0 - TEMP(T)
	G(5) = 10D0 - TEMP(T)
	G(6) = DY(6)

     RETURN
  END SUBROUTINE EF

  DOUBLE PRECISION FUNCTION INOCCULATE(T)

    DOUBLE PRECISION :: T

	! Function to innoculate within 0<T<1.
    IF (T < 1D0 .AND. T > 0D0) THEN
       INOCCULATE = 100D0
    ELSE
       INOCCULATE = 0D0
    END IF

	RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION TEMP(T)
     
	! Temperature function
    DOUBLE PRECISION :: T
	TEMP = (M-A) + A * 2 * (0.5D0 * (1D0 + COS(2D0 * Pi * (T-PHASE) / 365D0)))**POWER
    RETURN

  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DAYLIGHT(T)
     
	! Function used to calculate photoperiod given latitude
    DOUBLE PRECISION :: T,EPS,NUM,DEN
    REAL, PARAMETER :: Pi = 3.1415927D0, L = 51D0
    EPS = ASIN(0.39795D0 * COS(0.2163108D0 + 2 * ATAN(0.9671396D0 * TAN(0.00860D0 * (T-3.5D0)))))
    NUM = SIN(0.8333D0*Pi/180D0) + (SIN(L*Pi/180D0) * SIN(EPS))
    DEN = COS(L*Pi/180D0) * COS(EPS)
    DAYLIGHT = 24D0 - (24D0/Pi) * ACOS(NUM / DEN)
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DIAPAUSE_SPRING(PP)
     
	! Function to calculate diapause proportion in spring
    DOUBLE PRECISION :: PP
    DIAPAUSE_SPRING = 1D0 / (1D0 + EXP(5D0*(14D0-PP)))
    RETURN
    
  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DIAPAUSE_AUTUMN(PP)
     
	! Function to calculate diapause proportion in autumn
    DOUBLE PRECISION :: PP
    DIAPAUSE_AUTUMN = 1D0 / (1D0 + EXP(5D0*(13D0-PP)))
    RETURN
    
  END FUNCTION

  DOUBLE PRECISION FUNCTION BIRTH(DIAPAUSE,GONOTROPHICtime,ADU,T)
     
	! Function to calculate the birth rate
    DOUBLE PRECISION :: GONOTROPHICtime,EGGRAFT,ADU,T,DIAPAUSE
    EGGRAFT = DIAPAUSE*MAXEGG
    IF(ADU < 1D0 .AND. T > 1D0) THEN
      EGGRAFT = 0D0
    END IF
    BIRTH = EGGRAFT/GONOTROPHICtime
  
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION GONOTROPHIC(TEMP)
     
	! Function to calculate the length of the gonotrophic cycle
	DOUBLE PRECISION :: TEMP,GONOTROPHICRATE
	DOUBLE PRECISION :: ALPHA=0.1215D0,BETA=1.7612D0
    IF(TEMP .LE. 10) THEN
      GONOTROPHICRATE = 0D0
    ELSE
	  GONOTROPHICRATE = ALPHA * LOG10(TEMP-9)**BETA
    END IF
    IF (GONOTROPHICRATE > 0D0) THEN
      GONOTROPHIC = 1/GONOTROPHICRATE
    ELSE
      GONOTROPHIC = 1000D0
    END IF

    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHegg(TEMP)
    
	! Function to calculate the egg death rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0952D0,U4=22.88D0,U5=7D0
    
	DEATHegg = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHegg > 1D0) THEN
      DEATHegg = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHlar(TEMP)
    
	! Function to calculate the larval death rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0135D0,U4=20.125D0,U5=7D0
    
    DEATHlar = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHlar > 1D0) THEN
      DEATHlar = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHpup(TEMP)
    
	! Function to calculate the pupal death rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0135D0,U4=20.125D0,U5=7D0
    
	DEATHpup = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHpup > 1D0) THEN
      DEATHpup = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHadu(TEMP)
    
	!Function to calculate the adult death rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=2.166D-8,BETA=4.483D0
	
	DEATHadu = ALPHA*(TEMP**BETA)
    IF (DEATHadu < 0.01D0) THEN
      DEATHadu = 0.01D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION EGGMATURATION(TEMP)
     
	! Function to calculate the egg development rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.002621D0,BETA=1.718531D0
	
	EGGMATURATION = ALPHA*(TEMP**BETA)
	IF (EGGMATURATION < 0.016667D0) THEN
    EGGMATURATION = 0.016667D0
    END IF     

    RETURN

  END FUNCTION

  DOUBLE PRECISION FUNCTION LARMATURATION(TEMP)
     
	! Function to calculate the larval development rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.0013154D0,BETA=1.3882843D0
	
	LARMATURATION = ALPHA*(TEMP**BETA)
	IF (LARMATURATION < 0.016667D0) THEN
    LARMATURATION = 0.016667D0
    END IF    

    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION PUPMATURATION(TEMP)
     
	! Function to calculate the pupal development rate
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: ALPHA=0.0007109D0,BETAPAR=1.8865648D0
	
	PUPMATURATION = ALPHA*(TEMP**BETAPAR)
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

  ! Set the timespan over which to run the simulation and the frequency with which to output solution values.
  INTEGER, PARAMETER :: NOUT=913D0
  DOUBLE PRECISION, PARAMETER :: T0=0D0,TFINAL=912D0
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN= &
  (/ (T0+(I-1)*((TFINAL - T0)/(NOUT-1)), I=1,NOUT) /)
  
  TYPE(DDE_SOL) :: SOL
  TYPE(DDE_OPTS) :: OPTS
  
  ! Sets the maximum length of the solution queue which the solver retains in memory
  DOUBLE PRECISION :: MAXDELAY = 200D0
  CHARACTER (len=90) :: filename
  
  ! Options for the solver code (see solver code documentation)
  OPTS = DDE_SET(RE=1D-11,AE=1D-20,MAX_STEPS=1000000000,MAX_DELAY=MAXDELAY,TRIM_FREQUENCY=10000,DIRECTION=(/ 0,-1,1,-1,1,0 /))

  ! Calls the solver code
  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,TSPAN,EVENT_FCN=EF,OPTIONS=OPTS)

  ! Was the solver successful?
		IF (SOL%FLAG == 0) THEN

			! output solutions to file
			WRITE(filename, '( "output.dat" )' )
			OPEN(unit=11,file=filename)
			DO I = 1,SOL%NPTS
				WRITE(UNIT=11,FMT='(16E14.5E3)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
			END DO
			CLOSE(11) 

		ELSE
			
			! return an error
			PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
			SOL%FLAG

		END IF

  STOP
END PROGRAM basicmodel

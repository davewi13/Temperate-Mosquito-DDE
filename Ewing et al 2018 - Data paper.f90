MODULE define_DDEs

  IMPLICIT NONE

  ! Set number of equations, delays and event functions
  INTEGER, PARAMETER :: NEQN=14,NLAGS=6,NEF=3
  ! Set length of temperature datasets
  INTEGER, PARAMETER :: TEMPNwater = 17520, TEMPNair = 1460
  INTEGER :: K
  DOUBLE PRECISION :: PI=3.1415927D0,VOL=20D0
  ! Set predation parameters and egg raft size
  DOUBLE PRECISION :: densa=1.03D0,densr=0.021298D0,densh=0.043D0,UPS=19.841D0,SHP=2.4537D0,MAXEGG=200D0,B0=0.0031946D0,B1=0.0046884D0
  ! Set time points for temperature values
  double precision, DIMENSION(TEMPNwater) :: Xtempwater= (/ (K/24D0, K=1,TEMPNwater) /)
  double precision, DIMENSION(TEMPNair) :: Xtempair= (/ (K/2D0, K=1,TEMPNair) /)
  
CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

	! Set variables for DDE solver
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY
    INTEGER :: K,L
    ! Set values used in extracting temperature data and interpolating between values
    DOUBLE PRECISION :: ytempp1 = 1D0, ytemppn = 1D0
    DOUBLE PRECISION, DIMENSION(TEMPNwater) :: watertemp, watertemp2
    DOUBLE PRECISION, DIMENSION(TEMPNair) :: airtemp, airtemp2

    ! Define variable names
    DOUBLE PRECISION :: TEMPnowair,TEMPGCair
    DOUBLE PRECISION :: TEMPnowwater,TEMPEwater,TEMPLwater,TEMPPwater,TEMPELwater,TEMPLPwater,TEMPELPwater
    DOUBLE PRECISION :: BIRTHnow,BIRTHE,BIRTHEL,BIRTHELP
    DOUBLE PRECISION :: EGGMATnow,EGGMATE,EGGMATL,EGGMATEL,EGGMATLP,EGGMATELP
    DOUBLE PRECISION :: LARMATnow,LARMATL,LARMATP,LARMATLP
    DOUBLE PRECISION :: REt,RLt,RPt,RAt,MEt,MLt,MPt,DEt,DLt,DPt,DAt
    DOUBLE PRECISION :: dEdt,dLdt,dPdt,dAdt,dSEdt,dSLdt,dSPdt,dDEdt,dDLdt,dDELdt,dDPdt,dDLPdt,dDELPdt,dGCdt
    DOUBLE PRECISION :: E,LAR,PUP,ADU,SE,SL,SP,DE,DL,DP,GC
    DOUBLE PRECISION :: INOC
    DOUBLE PRECISION :: DEATHeggnow,DEATHeggE,DEATHlarnow,DEATHlarL,DEATHpupnow,DEATHpupP,DEATHadunow
    DOUBLE PRECISION :: PUPMATnow,PUPMATP
    DOUBLE PRECISION :: DIAPAUSEE,DIAPAUSEEL,DIAPAUSEELP,DIAPAUSEnow
    DOUBLE PRECISION :: PPnow,PPE,PPEL,PPELP
    DOUBLE PRECISION :: GCnow,GCC

	! Read in temperature data
	OPEN (UNIT=7, FILE="Butt4hourlytemps-winterstart-twoyears.csv", STATUS="OLD", ACTION="READ")
    READ(7,*)
    DO K = 1, TEMPNwater
      READ (7,*) watertemp(K)
    END DO
    CLOSE(7)

    OPEN (UNIT=7, FILE="Airminmaxmetsite-winterstart-twoyears.csv", STATUS="OLD", ACTION="READ")
    READ(7,*)
    DO L = 1, TEMPNair
      READ (7,*) airtemp(L)
    END DO
    CLOSE(7)

	! Interpolate between temperature values using spline
    CALL cubic_spline_air(xtempair, airtemp, TEMPnair, ytempp1, ytemppn, airtemp2)
    CALL cubic_spline_water(xtempwater, watertemp, TEMPnwater, ytempp1, ytemppn, watertemp2)

    ! Give names to responses of system of equations    
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
    GC = Y(14)

    ! Extract temperature values at various time points
	TEMPnowwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T)
    TEMPnowair = splint(xtempair,airtemp,airtemp2,TEMPnair,T)
    TEMPEwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T-DE)
    TEMPLwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T-DL)
    TEMPPwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T-DP)
    TEMPELwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T-DL-Z(8,4))
    TEMPELPwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T-DP-Z(9,6)-Z(8,5))
    TEMPLPwater = splint(xtempwater,watertemp,watertemp2,TEMPnwater,T-DP-Z(9,6))
    TEMPGCair = splint(xtempair,airtemp,airtemp2,TEMPnair,T-GC)

    ! Define temperature values for time points before T=0
    IF ((T-DE) .LE. 0) THEN
      TEMPEwater = 5D0
    END IF
    IF ((T-DL) .LE. 0) THEN
      TEMPLwater = 5D0
    END IF
    IF ((T-DP) .LE. 0) THEN
      TEMPPwater = 5D0
    END IF
    IF ((T-DL-Z(8,4)) .LE. 0) THEN
      TEMPELwater = 5D0
    END IF
    IF ((T-DP-Z(9,6)-Z(8,5)) .LE. 0) THEN
      TEMPELPwater = 5D0
    END IF
    IF ((T-DP-Z(9,6)) .LE. 0) THEN
      TEMPLPwater = 5D0
    END IF
    IF ((T-GC) .LE. 0) THEN
      TEMPGCair = 5D0
    END IF
    
	! Calculate photoperiod values for a range of time points
    PPnow = DAYLIGHT(T)
    PPE = DAYLIGHT(T-DE)
    PPEL = DAYLIGHT(T-DL-Z(8,4))
    PPELP = DAYLIGHT(T-DP-Z(9,6)-Z(8,5))

	! Calculate gonotrophic cycle length for various time points
    GCnow = GONOTROPHIC(TEMPnowair)
    GCC = GONOTROPHIC(TEMPGCair)

	! Calculate diapause percentage for given photoperiod values.
    ! Use these diapause percentages to caluclate birth rates.
    ! Do this for both spring and autumn diapause thresholds.
    IF (DAYLIGHT(T) > DAYLIGHT(T-1)) THEN
      DIAPAUSEnow = DIAPAUSE_SPRING(PPnow)
      DIAPAUSEE = DIAPAUSE_SPRING(PPE)
      DIAPAUSEEL = DIAPAUSE_SPRING(PPEL)
      DIAPAUSEELP = DIAPAUSE_SPRING(PPELP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GC)
      BIRTHE = BIRTH(DIAPAUSEE,Z(14,1))
      BIRTHEL = BIRTH(DIAPAUSEEL,Z(14,2))
      BIRTHELP = BIRTH(DIAPAUSEELP,Z(14,3))
    ELSE
      DIAPAUSEnow = DIAPAUSE_AUTUMN(PPnow)
      DIAPAUSEE = DIAPAUSE_AUTUMN(PPE)
      DIAPAUSEEL = DIAPAUSE_AUTUMN(PPEL)
      DIAPAUSEELP = DIAPAUSE_AUTUMN(PPELP)
      BIRTHnow = BIRTH(DIAPAUSEnow,GC)
      BIRTHE = BIRTH(DIAPAUSEE,Z(14,1))
      BIRTHEL = BIRTH(DIAPAUSEEL,Z(14,2))
      BIRTHELP = BIRTH(DIAPAUSEELP,Z(14,3))
    END IF

	! Calculate death rates for each life stage
	DEATHeggnow = DEATHegg(TEMPnowwater)
    DEATHeggE = DEATHegg(TEMPEwater)

    DEATHlarnow = DEATHlar(TEMPnowwater)
    DEATHlarL = DEATHlar(TEMPLwater)

    DEATHpupnow = DEATHpup(TEMPnowwater)
    DEATHpupP = DEATHpup(TEMPPwater)
    
    DEATHadunow = DEATHadu(TEMPnowair,GC,T)

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

	! Innoculate the system with a given number of adults.
    INOC = INOCCULATE(T)

	! Delay differential equations for immature stage durations
    dDEdt = 1D0 - EGGMATnow/EGGMATE
    dDLdt = 1D0 - LARMATnow/LARMATL
    dDPdt = 1D0 - PUPMATnow/PUPMATP

    ! Delay differential equation for gonotrophic cycle duration 
    dGCdt = 1D0 - GCnow/GCC

    ! Delay differential equations for stage durations referenced back through previous stages
    dDELdt = (1D0 - dDLdt) * (1D0 - EGGMATL/EGGMATEL)
    dDLPdt = (1D0 - dDPdt) * (1D0 - LARMATP/LARMATLP)
    dDELPdt = (1D0 - dDPdt - dDLPdt) * (1D0 - EGGMATLP/EGGMATELP)
          
	! Recruitment equations for each stage
	REt = BIRTHnow * ADU
	RLt = BIRTHE * Z(4,1) * SE * EGGMATnow/EGGMATE
    RPt = BIRTHEL * Z(4,2) * Z(5,4) * SL * LARMATnow/LARMATL * EGGMATL/EGGMATEL
    RAt = BIRTHELP * Z(4,3) * Z(5,5) * Z(6,6) * SP * PUPMATnow/PUPMATP * LARMATP/LARMATLP * EGGMATLP/EGGMATELP + INOC
    
	! Maturation equations for each stage
    MEt = RLt
    MLt = RPt
    MPt = BIRTHELP * Z(4,3) * Z(5,5) * Z(6,6) * SP * PUPMATnow/PUPMATP * LARMATP/LARMATLP * EGGMATLP/EGGMATELP	

	! Death equations for each stage
	DEt = DEATHeggnow * E
    DLt = (densa*densr*(((1D0+COS(2D0*PI*(T-182.5D0-UPS)/365D0))/2D0)**SHP)*(LAR/(VOL/5D0))/(1D0+densh*densa*LAR/(VOL/5D0))&
    + DEATHlarnow + B0*EXP(B1*LAR/VOL))*LAR
    DPt = DEATHpupnow * PUP
    DAt = DEATHadunow * ADU

    ! Balance equations for each life stage
	dEdt = REt - MEt - DEt
	dLdt = RLt - MLt - DLt
    dPdt = RPt - MPt - DPt
	dAdt = RAt - DAt
    
	! Survival equations for the immature stages
    dSEdt = SE * ((EGGMATnow * DEATHeggE / EGGMATE) - DEATHeggnow)
	dSLdt =SL*(((densa*densr*(((1D0+COS(2D0*PI*(T-DL-182.5D0-UPS)/365D0))/2D0)**SHP)*Z(2,4)/(VOL/5D0)/(1D0+densa*densh*Z(2,4)&
    /(VOL/5D0)))+DEATHlarL + B0*EXP(B1*Z(2,4)/VOL)) * (1-dDLdt) - (densa*densr*(((1D0+COS(2D0*PI*(T-182.5D0-UPS)/365D0))/2D0)**SHP)&
    *LAR/(VOL/5D0) / (1D0+densa*densh*LAR/(VOL/5D0))) - DEATHlarnow - B0*EXP(B1*LAR/VOL))
    dSPdt = SP * ((PUPMATnow * DEATHpupP / PUPMATP) - DEATHpupnow)

    ! Derivatives for the integrator:
    DY = (/ dEdt, dLdt, dPdt, dAdt, dSEdt, dSLdt, dSPdt, dDEdt, dDLdt, dDPdt, dDELdt, dDLPdt, dDELPdt, dGCdt /)

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
	BVAL(1) = T-Y(8)
    ! T - Lardelay(T) - Eggdelay(T-Lardelay(T))				
    BVAL(2) = T-Y(9)-Y(11)
    !T - Pupdelay(T) - Lardelay(T-Pupdelay(T)) - Eggdelay(T-Pupdelay(T)-Lardelay(T-Pupdelay(T)))
    BVAL(3) = T-Y(10)-Y(12)-Y(13)	
    !T - Lardelay(T)
    BVAL(4) = T-Y(9)				
    !T - Pupdelay(T) - Lardelay(T-Pupdelay(T))
    BVAL(5) = T-Y(10)-Y(12)			
    !T - Pupdelay(T)
    BVAL(6) = T-Y(10)				
   	  
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE HISTORY(T,Y)
    DOUBLE PRECISION :: T,TEMPhist
    DOUBLE PRECISION, DIMENSION(2*NEQN) :: Y
    INTENT(IN)  :: T
    INTENT(OUT) :: Y
      
	! Set the temperatures for T < 0
    TEMPhist = 5D0
    
	! Set historical values for all stages to be zero
    Y(1) = 0D0
    Y(2) = 0D0
    Y(3) = 0D0
    Y(4) = 0D0

    ! Calculate historical survival rates based on temperature
    Y(5) = EXP(-DEATHegg(TEMPhist)*(1D0/EGGMATURATION(TEMPhist)))
    Y(6) = EXP(-DEATHlar(TEMPhist)*(1D0/LARMATURATION(TEMPhist)))
    Y(7) = EXP(-DEATHpup(TEMPhist)*(1D0/PUPMATURATION(TEMPhist)))

	! Calculate historical development rates based on temperature
    Y(8) = 1D0/EGGMATURATION(TEMPhist)
    Y(9) = 1D0/LARMATURATION(TEMPhist)
    Y(10) = 1D0/PUPMATURATION(TEMPhist)
    Y(11) = 1D0/EGGMATURATION(TEMPhist)
    Y(12) = 1D0/LARMATURATION(TEMPhist)
    Y(13) = 1D0/EGGMATURATION(TEMPhist)
    Y(14) = 1D0/GONOTROPHIC(TEMPhist)       
    

    RETURN
  END SUBROUTINE HISTORY

   SUBROUTINE cubic_spline_air(x,airtemp,n,airtempp1,airtemppn,airtemp2)
    INTEGER n,NMAX
    double precision airtempp1,airtemppn,x(n),airtemp(n),airtemp2(n)
    PARAMETER (NMAX=1460)
    
	!Code taken from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
    !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
    !x1 < x2 < ::: < xN, and given values yp1 and ypn for the first derivative of the interpolating
    !function at points 1 and n, respectively, this routine returns an array y2(1:n) of
    !length n which contains the second derivatives of the interpolating function at the tabulated
    !points xi. If yp1 and/or ypn are equal to 1  1030 or larger, the routine is signaled to set
    !the corresponding boundary condition for a natural spline, with zero second derivative on
    !that boundary.
    
    !Parameter: NMAX is the largest anticipated value of n.
    INTEGER i,k
    double precision p,qn,sig,un,u(NMAX)
    if (airtempp1.gt..99e30) then !The lower boundary condition is set either to be
    airtemp2(1)=0. !\natural"
    u(1)=0.
    else !or else to have a specified first derivative.
    airtemp2(1)=-0.5
    u(1)=(3./(x(2)-x(1)))*((airtemp(2)-airtemp(1))/(x(2)-x(1))-airtempp1)
    endif
    do i=2,n-1 
    !This is the decomposition loop of the tridiagonal
    !algorithm. y2 and u are used for temporary
    !storage of the decomposed factors.
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*airtemp2(i-1)+2.
    airtemp2(i)=(sig-1.)/p
    u(i)=(6.*((airtemp(i+1)-airtemp(i))/(x(i+1)-x(i))-(airtemp(i)-airtemp(i-1)) &
    /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (airtemppn.gt..99e30) then !The upper boundary condition is set either to be
    qn=0. !\natural"
    un=0.
    else !or else to have a specified first derivative.
    qn=0.5
    un=(3./(x(n)-x(n-1)))*(airtemppn-(airtemp(n)-airtemp(n-1))/(x(n)-x(n-1)))
    endif
    airtemp2(n)=(un-qn*u(n-1))/(qn*airtemp2(n-1)+1.)
    do k=n-1,1,-1 !This is the backsubstitution loop of the tridiagonal algorithm.
    airtemp2(k)=airtemp2(k)*airtemp2(k+1)+u(k) 
    enddo
    return

   END SUBROUTINE

   SUBROUTINE cubic_spline_water(x,watertemp,n,watertempp1,watertemppn,watertemp2)
    INTEGER n,NMAX
    double precision watertempp1,watertemppn,x(n),watertemp(n),watertemp2(n)
    PARAMETER (NMAX=131400)
    
	!Code taken from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
    !Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
    !x1 < x2 < ::: < xN, and given values yp1 and ypn for the rst derivative of the interpolating
    !function at points 1 and n, respectively, this routine returns an array y2(1:n) of
    !length n which contains the second derivatives of the interpolating function at the tabulated
    !points xi. If yp1 and/or ypn are equal to 1  1030 or larger, the routine is signaled to set
    !the corresponding boundary condition for a natural spline, with zero second derivative on
    !that boundary.
    
    !Parameter: NMAX is the largest anticipated value of n.
    INTEGER i,k
    double precision p,qn,sig,un,u(NMAX)
    if (watertempp1.gt..99e30) then !The lower boundary condition is set either to be
    watertemp2(1)=0. !\natural"
    u(1)=0.
    else !or else to have a specified first derivative.
    watertemp2(1)=-0.5
    u(1)=(3./(x(2)-x(1)))*((watertemp(2)-watertemp(1))/(x(2)-x(1))-watertempp1)
    endif
    do i=2,n-1 
    !This is the decomposition loop of the tridiagonal
    !algorithm. y2 and u are used for temporary
    !storage of the decomposed factors.
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*watertemp2(i-1)+2.
    watertemp2(i)=(sig-1.)/p
    u(i)=(6.*((watertemp(i+1)-watertemp(i))/(x(i+1)-x(i))-(watertemp(i)-watertemp(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    if (watertemppn.gt..99e30) then !The upper boundary condition is set either to be
    qn=0. !\natural"
    un=0.
    else !or else to have a specified first derivative.
    qn=0.5
    un=(3./(x(n)-x(n-1)))*(watertemppn-(watertemp(n)-watertemp(n-1))/(x(n)-x(n-1)))
    endif
    watertemp2(n)=(un-qn*u(n-1))/(qn*watertemp2(n-1)+1.)
    do k=n-1,1,-1 !This is the backsubstitution loop of the tridiagonal algorithm.
    watertemp2(k)=watertemp2(k)*watertemp2(k+1)+u(k) 
    enddo
    return
   
   END SUBROUTINE

   double precision FUNCTION SPLINT(xa,ytempa,ytemp2a,tempn,x)
    INTEGER tempn
    double precision xa(tempn),ytemp2a(tempn),ytempa(tempn)
    double precision x
    !Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
    !xai 's in order), and given the array y2a(1:n), which is the output from spline above,
    !and given a value of x, this routine returns a cubic-spline interpolated value y.
    INTEGER k,khi,klo
    double precision a1,b,h1

    !We will find the right place in the table by means of bisection.
    !This is optimal if sequential calls to this routine are at random
    !values of x. If sequential calls are in order, and closely
    !spaced, one would do better to store previous values of
    !klo and khi and test if they remain appropriate on the
    !next call.
    
    klo=1 
    khi=tempn
    do
    if (khi-klo.le.1) exit
    k=(khi+klo)/2
    if(xa(k).gt.x) then
        khi=k
    else
        klo=k
    endif
    enddo
    !klo and khi now bracket the input value of x.
    
    h1=xa(khi)-xa(klo)
    if (h1.eq.0) pause 'bad xa input in splint' !The xa's must be distinct.
    a1=(xa(khi)-x)/h1 !Cubic spline polynomial is now evaluated.
    b=(x-xa(klo))/h1
    if (x .le. 0) then
      splint=0.1d0
      else
    splint=a1*ytempa(klo)+b*ytempa(khi)+ ((a1**3-a1)*ytemp2a(klo)+(b**3-b)*ytemp2a(khi))*(h1**2)/6.
    end if
    return
  END FUNCTION


  DOUBLE PRECISION FUNCTION INOCCULATE(T)

    DOUBLE PRECISION :: T

	! Set inoculation value
    IF (T < 1D0 .AND. T > 0D0) THEN
       INOCCULATE = 5000D0
    ELSE
       INOCCULATE = 0D0
    END IF

	RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DAYLIGHT(T)
     
    DOUBLE PRECISION :: T,EPS,NUM,DEN
    REAL, PARAMETER :: Pi = 3.1415927D0, L = 51.6D0

    ! Calculate photoperiod for time, T
    EPS = ASIN(0.39795D0 * COS(0.2163108D0 + 2 * ATAN(0.9671396D0 * TAN(0.00860D0 * (T-185.5D0)))))
    NUM = SIN(0.8333D0*Pi/180D0) + (SIN(L*Pi/180D0) * SIN(EPS))
    DEN = COS(L*Pi/180D0) * COS(EPS)
    DAYLIGHT = 24D0 - (24D0/Pi) * ACOS(NUM / DEN)
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DIAPAUSE_SPRING(PP)
     
    DOUBLE PRECISION :: PP
    ! Set spring diapause threshold
    DIAPAUSE_SPRING = 1D0 / (1D0 + EXP(5D0*(13.7D0-PP)))
    RETURN
    
  END FUNCTION
  
  DOUBLE PRECISION FUNCTION DIAPAUSE_AUTUMN(PP)
     
    DOUBLE PRECISION :: PP
    !Set autumn diapause threshold
    DIAPAUSE_AUTUMN = 1D0 / (1D0 + EXP(3.5D0*(15D0-PP)))
    RETURN
    
  END FUNCTION

  DOUBLE PRECISION FUNCTION BIRTH(DIAPAUSE,GONOTROPHICtime)
     
    DOUBLE PRECISION :: GONOTROPHICtime,EGGRAFT,DIAPAUSE

	!Set birth rate
    EGGRAFT = DIAPAUSE*MAXEGG*0.5D0
    BIRTH = EGGRAFT/GONOTROPHICtime
  
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION GONOTROPHIC(TEMP)
     
	DOUBLE PRECISION :: TEMP
	DOUBLE PRECISION :: KG=0.2024D0,QG=74.48D0,BG=0.2456D0
    
	!Set gonotrophic cycle rate
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

  DOUBLE PRECISION FUNCTION DEATHegg(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	! Set death rate of eggs
	DEATHegg = U3 * EXP(((TEMP-U4)/U5)**2)
	IF (DEATHegg > 1D0) THEN
      DEATHegg = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHlar(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	! Set death rate of larvae
	DEATHlar = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHlar > 1D0) THEN
      DEATHlar = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHpup(TEMP)
    
	DOUBLE PRECISION :: TEMP
    DOUBLE PRECISION :: U3=0.0157D0,U4=20.5D0,U5=7D0
    
	! Set death rate of pupae
	DEATHpup = U3 * EXP(((TEMP-U4)/U5)**2)
    IF (DEATHpup > 1D0) THEN
      DEATHpup = 1D0
    END IF
        
    RETURN
  END FUNCTION

  DOUBLE PRECISION FUNCTION DEATHadu(TEMP,GONOTROPHICtime,T)
    
	DOUBLE PRECISION :: TEMP,GONOTROPHICtime,T
    DOUBLE PRECISION :: ALPHA=2.166D-8,BETA=4.483D0,PI=3.1415927D0,MULTIPLIER=8D0,SIGMASQ=4D0
	
	!Set adult death rate
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
	
	! Set egg development rate
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
	
	! Set larval development rate
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
	
	! Set pupal development rate
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
  INTEGER, PARAMETER :: NOUT=731D0
  DOUBLE PRECISION, PARAMETER :: T0=0D0,TFINAL=730D0
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN= &
  (/ (T0+(I-1)*((TFINAL - T0)/(NOUT-1)), I=1,NOUT) /)

  TYPE(DDE_SOL) :: SOL
  TYPE(DDE_OPTS) :: OPTS
  
  ! Set length of maximum delay
  DOUBLE PRECISION :: MAXDELAY = 200D0
  CHARACTER (len=90) :: filename
  
  ! Set options for DDE solver
  OPTS = DDE_SET(RE=1D-5,AE=1D-5,MAX_STEPS=100000000,MAX_DELAY=MAXDELAY,TRIM_FREQUENCY=10000)  

  ! Call DDE solver code
  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,TSPAN,OPTIONS=OPTS)
  
  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN
	 ! Output solution values
     WRITE(filename, '( "abc_post_median.dat" )' )
     OPEN(unit=11,file=filename)
     DO I = 1,SOL%NPTS
       WRITE(UNIT=11,FMT='(16E14.5E3)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
     END DO
     CLOSE(11)

  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  END IF


  STOP
END PROGRAM basicmodel

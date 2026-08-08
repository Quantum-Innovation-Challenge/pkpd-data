$PROBLEM NIELSEN 2007 semimech_PKPD_antibiotics

$INPUT ID C TYPE EVID CMT TIME DV AMT L2 DIL

$DATA Simulated_semimech_PKPD_antibiotics.tab IGNORE=@

$SUBS  ADVAN9 TOL7

$MODEL
NCOMP=4
COMP=(S)    ;Bacteria in growth state (Sucseptible)
COMP=(R)    ;Bacteria in resting state (Resistant)
COMP=(ABS)  ;Total drug concentration
COMP=(PK)   ;Active drug concentration

$PK
" FIRST
"  COMMON /PRCOMG/ IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"  INTEGER IDUM1,IDUM2,IMAX,IDUM4,IDUM5
"  IMAX=10000000

;Store the estimated subpopulation of the current individual in EST
EST=MIXEST

TVKGS   = THETA(1)
TVKK    = THETA(2)

IF(TYPE.EQ.2) TVEMAX  = THETA(3) ;PCG
IF(TYPE.EQ.3) TVEMAX  = THETA(4) ;ERYTRO
IF(TYPE.EQ.4) TVEMAX  = THETA(5) ;CEFUR
IF(TYPE.EQ.5) TVEMAX  = THETA(6) ;VANCO
IF(TYPE.EQ.6) TVEMAX  = THETA(7) ;MOXI

IF(TYPE.EQ.2) TVEC50  = THETA(8)*0.001  ;PCG
IF(TYPE.EQ.3) TVEC50  = THETA(9)*0.001  ;ERYTRO
IF(TYPE.EQ.4) TVEC50  = THETA(10)*0.001 ;CEFUR
IF(TYPE.EQ.5) TVEC50  = THETA(11)*0.001 ;VANCO
IF(TYPE.EQ.6) TVEC50  = THETA(12)*0.001 ;MOXI

IF(TYPE.EQ.2) TVGAM   = THETA(13) ;PCG
IF(TYPE.EQ.3) TVGAM   = THETA(14) ;ERYTRO
IF(TYPE.EQ.4) TVGAM   = THETA(15) ;CEFUR
IF(TYPE.EQ.5) TVGAM   = THETA(16) ;VANCO
IF(TYPE.EQ.6) TVGAM   = THETA(17) ;MOXI

IF(TYPE.EQ.2) TVKD  = THETA(22)  ;PCG
IF(TYPE.EQ.3) TVKD  = THETA(23)  ;ERYTRO
IF(TYPE.EQ.4) TVKD  = THETA(24)  ;CEFUR
IF(TYPE.EQ.5) TVKD  = THETA(25)  ;VANCO
IF(TYPE.EQ.6) TVKD  = THETA(26)  ;MOXI

IF(TYPE.EQ.2) TVKEO  = THETA(27) ;PCG
IF(TYPE.EQ.3) TVKEO  = THETA(28) ;ERYTRO
IF(TYPE.EQ.4) TVKEO  = THETA(29) ;CEFUR
IF(TYPE.EQ.5) TVKEO  = THETA(30) ;VANCO
IF(TYPE.EQ.6) TVKEO  = THETA(31) ;MOXI

TVBMAX  = THETA(18)*1000000

;PERS = Pre Existing Resistant, resting phase
;Estimate the fraction of bacteria in resting phase in the current mixture state
IF(MIXNUM.EQ.1) PERS=AMT*THETA(19) ;Log growth, no bacteria (0 FIX) in resting state at time=0
IF(MIXNUM.EQ.2) PERS=AMT*THETA(20) ;Stationary phase, estimate fraction of bacteria in resting state at time=0

;Set initial conditions
A_0(2)=PERS ;Bacteria in resting state
A_0(3)=C    ;Total drug concentration
A_0(4)=0    ;Active drug concentration

IF(MIXNUM.EQ.1) F1=(1-THETA(19))  ;Log growth, no bacteria (0 FIX) in resting state at time=0
IF(MIXNUM.EQ.2) F1=(1-THETA(20))  ;Stationary phase, estimate fraction of bacteria in resting state at time=0

KGS    = TVKGS*EXP(ETA(1))
KK     = TVKK
EMAX   = TVEMAX
EC50   = TVEC50
GAM    = TVGAM
BMAX   = TVBMAX
KD     = TVKD
KEO    = TVKEO

;Estimate the mixture probabillities of the experiment for the two subcategories; 1 Log growth, 2 Stationary phase
$MIX
P(1)=THETA(21)
P(2)=1.-THETA(21)
NSPOP=2

$DES

CABS=(A(3)) ;Total drug
IF(CABS.LT.0.001) CABS=0

CEFF=(A(4)) ;Active drug
IF(CEFF.LT.0.00001) CEFF=0

;Drug effect
DRUGS = EMAX*(CEFF)**GAM/(EC50**GAM+(CEFF)**GAM)

;Tranfer rate from growing to resting state
FEED=(KGS-KK)/BMAX

;Tranfer rate from growing to resting state scaled to the total number of bacteria
SR=FEED*((A(1))+(A(2)))

;Tranfer rate from resting to growing state is not experimentally observed
RS=0

;Differential equation system
DADT(1)=KGS*(A(1))-(KK+DRUGS)*(A(1)) - SR*(A(1)) + RS*(A(2)) ;Bacteria in growing state
DADT(2)=-KK*(A(2)) + SR*(A(1)) - RS*(A(2))                   ;Bacteria in resting state
DADT(3)=-KD*(A(3))                                           ;Total drug concentration
DADT(4)=KEO*(A(3))-KEO*(A(4))                                ;Active drug concentration

$ERROR
A1=(A(1))  ;Bacteria in growing state
A2=(A(2))  ;Bacteria in resting state
ATOT=A1+A2 ;Total number of bacteria
A3=(A(3))  ;Total drug concentration
A4=(A(4))  ;Active drug concentration

DEL=0
IF(F.EQ.0) DEL=1 ;Guard against log(0)

IPRED=LOG(ATOT+DEL)                ;Log transform data

IRES  = DV-IPRED                   ;Calculate residual

W=1                                ;Constant additative error
IWRES = IRES/W

IF(DIL.EQ.1) Y=IPRED+EPS(1)+EPS(2) ;Error replicate 1
IF(DIL.EQ.2) Y=IPRED+EPS(1)+EPS(3) ;Error replicate 2
IF(DIL.EQ.3) Y=IPRED+EPS(1)+EPS(4) ;Error replicate 3

$THETA  (0.5,1.35,5)  ; GS
$THETA  (0.05,0.18,1) ; KK

$THETA  (0.5,2.4,5)   ; EMAX_PCG
$THETA  (0.5,2.1,5)   ; EMAX_ERYTRO
$THETA  (0.5,3.3,5)   ; EMAX_CEFUR
$THETA  (0.5,1.4,5)   ; EMAX_VANCO
$THETA  (0.5,3.2,5)   ; EMAX_MOXI

$THETA  (1,4.6,10)    ; E50_PCG
$THETA  (1,28,50)     ; E50_ERYTRO
$THETA  (1,8.2,50)    ; E50_CEFUR
$THETA  (1,384,400)   ; E50_VANCO
$THETA  (1,75,200)    ; E50_MOXI

$THETA  (0,1,5)       ; GAM_PCG
$THETA  (0,1,5)       ; GAM_ERYTRO
$THETA  (0,1,5)       ; GAM_CEFUR
$THETA  20 FIX        ; GAM_VANCO
$THETA  (0,1,5)       ; GAM_MOXI
 
$THETA  (300,420,600)      ; BMAX
$THETA  0 FIX              ; FRAKTION_PERS_1
$THETA  (0.0001,0.05,0.99) ; FRAKTION_PERS_2
$THETA  (0.001,0.74,1)     ; FRAKTION_1

;Drug degradation rates previosly determined
$THETA  0.020 FIX     ; KD_PCG
$THETA  0 FIX         ; KD_ERYTRO
$THETA  0.026 FIX     ; KD_CEFUR
$THETA  0 FIX         ; KD_VANCO
$THETA  0 FIX         ; KD_MOXI

$THETA  (0,1,100)     ; KEO_PCG
$THETA  100 FIX       ; KEO_ERYTRO
$THETA  (0,1,100)     ; KEO_CEFUR
$THETA  100 FIX       ; KEO_VANCO
$THETA  (0,1,100)     ; KEO_MOXI

;No IIV variabillity assumed
$OMEGA  0 FIX         ; ETA_EMAX

$SIGMA  0.98           ; RES
$SIGMA  BLOCK(1) 0.47  ; RRES1
$SIGMA  BLOCK(1) SAME  ; RRES2
$SIGMA  BLOCK(1) SAME  ; RRES3

$EST METH=1 INTER SIGDIG=3  MAXEVAL=0 POSTHOC  PRINT=1
;$SIM (12345) ONLYSIM
;$COV

$TABLE ID TYPE TIME C IPRED IWRES NOPRINT ONEHEADER FILE=sdtab_semimech_PKPD_antibiotics
$TABLE ATOT ID TYPE TIME C A3 A4 KGS SR RS BMAX EST KK EMAX EC50 GAM FEED IPRED NOPRINT ONEHEADER FILE=patab_semimech_PKPD_antibiotics

;; 1. Based on: 43
;;Date of run 2003-03-13
;;1.Rerun of run38 but eta changed to additive
;;2.Structural model
;;Two-compartment pk-model and effect in central compartment
;;3.Covariate model
;;Not applicable
;;4.Interindividual variability 
;;Additative IIV on Baseline
;;Exponential IIV on EC50
;;5.Interoccasion variability
;;Not applicable
;;6.Residual variability 
;;Not applicable 
;;7.Estimation
;;LAPLACE NUM
;;8.Other
;;Probabalistic PD model
$PROBLEM    ATM-2 PD modelling of receptor expression high or low
$INPUT      ID CENT PANO DAT2=DROP TIME AMT DV FLAG SEX AGE WT HT BMI
            ICL IV1 IQ IV2 EVID CMT=DROP MDV
$DATA      Simulated_PKPD_Tcell_mAb_multiplesclerosis.csv IGNORE=#
$SUBROUTINE ADVAN7 TRANS1
$MODEL      COMP=(CENT DEFDOSE DEFOBS) COMP=PERI
$PK
       CL = ICL
       V1 = IV1
       V2 = IV2
       Q = IQ

       K10 = CL/V1
       K12 = Q/V1
       K21 = Q/V2

       S1 = V1

$ERROR
CP    =A(1)/V1
; ------------------ PD -------------
; Baseline values
 B1  =THETA(1)+ETA(1)
; Dose-effect relationship
 EMAX =THETA(2)
 EC50 =THETA(3)*EXP(ETA(2))
 CP   =F
 DRUG =F*EMAX/(F+EC50)
; Logits for Y>=4, Y>=5
 A1   =B1+DRUG
 C1   =EXP(A1)
; Probabilities for Y>=4
 P1   =C1/(1+C1)
; Probabilities for Y=4, Y=5
PA   =1-P1
PB   =P1
; Select appropriate P(Y=m)
IF(DV.EQ.4) Y=PB
IF(DV.EQ.5) Y=PA
$THETA  -5.84 FIX           ; BASE
$THETA (0,14.2323023216728) ; EMAX
$THETA (0,487.552840587381) ; EC50
$OMEGA  0  FIX              ; BASE
$OMEGA 0.29386601136254     ; EC50
$ESTIMATION NOABORT NUMERICAL MAX=0 PRINT=5 METHOD=COND LAPLACE LIKE
; $COVARIANCE
; $TABLE      ID TIME DRUG AMT CP NOPRINT ONEHEADER FILE=sdtab43
; $TABLE      ID EMAX EC50 ETA1 NOPRINT ONEHEADER FILE=patab43
; $TABLE      ID SEX NOPRINT ONEHEADER FILE=catab43
; $TABLE      ID AGE WT HT BMI AMT NOPRINT ONEHEADER FILE=cotab43


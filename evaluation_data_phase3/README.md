# Quantum computing innovation challenge

Follow-up case with new data set

Challenge subject

A phase 1 placebo-controlled clinical trial testing three dose levels (0.1 mg, 0.3 mg, and 1 mg) of a new compound dosed once daily for three weeks has just been completed. Data from 48 subjects is available (36 on active treatment and 12 on placebo). The subjects weigh between 50 and 130 kg. Pharmacokinetic (PK) data (measurements of compound concentration) is available for the subjects randomized to active treatment, and pharmacodynamic (PD) data (measurement of a biomarker) is available for all subjects. Prior knowledge about the structure of a PK-PD model that can describe the observed data is limited, but the data shows that the new compound has a long half-life and that it can increase the level of the biomarker.

Tasks to solve

Your task is to develop a quantum-enhanced model to answer the following questions:

What is the daily dose level (in whole multiples of 0.1 mg) that ensures that 90% of all subjects in a population similar to the one studied in the phase 1 trial achieve increase of the biomarker above a clinically relevant threshold (10 ng/mL) throughout a 24-hour dosing interval at steady-state?
    
Which weekly dose level (in whole multiples of 1 mg) has the same effect over a 168-hour dosing interval at steady-state, if the compound was dosed once-weekly?

Dataset explanations and metadata:

Dataset shape: 2820 rows x 10 columns
    
Dataset columns:
    
1. ID: Subject identifier.
	
2. BW: Body weight of the subject in kg.
	
3. DOSE: Dose level in mg. This column represents the amount of drug administered to the subject. It is typically measured in units like milligrams (mg) or micrograms (µg). The dose is a critical factor in determining the drug’s pharmacokinetic and pharmacodynamic properties.
	
4. TIME: Time in hours. Indicates the time elapsed since the start of the first drug administration. Time is typically measured in hours or minutes and is essential for plotting concentration-time profiles.
	
5. DV (Dependent Variable): Compound concentration (mg/L) for DVID=1. Biomarker level (ng/mL) for DVID=2. This column usually represents observed data, such as drug concentration in plasma or another biological matrix. It can also refer to the measurement of a biomarker or response variable affected by the drug, such as blood pressure or heart rate.
	
6. EVID (Event ID): This is an event identifier used in NONMEM (a common software in pharmacometrics). It signifies the type of event occurring:

   EVID = 0 for observation events (e.g., a concentration measurement),

   EVID = 1 for dosing events.
	    
7. MDV (Missing Dependent Variable): This value indicates whether the dependent variable (DV) is missing. An MDV of 1 means the DV value is missing, while 0 means it is present.
	
8. AMT (Amount): Stands for the actual dose amount administered, especially applicable during infusion dosing. AMT will be zero for observation records.
	
9. CMT (Compartment): This denotes the compartment where the event (e.g., dosing or sampling) occurs. In PK models, different compartments (like central and peripheral) help to describe drug kinetics.
	
10. DVID (Dependent Variable ID): This identifier helps distinguish between different types of DVs in the dataset. For example, you might have multiple measurement types such as concentrations and biomarkers.
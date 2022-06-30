# Case Study 2: Real-time Process Monitoring for Fed-Batch Mammalian Cell Culture Bioreactors Using an In-Line Raman Spectroscopy Prob

## Introduction
Raman spectroscopy measures the amount of light scattered inelastically at different frequencies by molecular vibrations, which results in a linear combination of individual molecular fingerprints (aka single - component spectra) with high chemical specificity.
An in-situ Raman spectroscopy probe or instrument can be used to online monitor cell culture process parameters simultaneously and continuously, including cell nutrients and metabolites (e.g., glucose, lactate, ammonium, viable cell density, and whole cell density).
In this case study, students will learn how to build multivariate statistical models for relating Raman spectra to glucose and lactate concentrations and the use of such models in troubleshooting of operational problems. The case study is derived from real industrial problems in which Raman spectroscopy is used for real-time monitoring of multiple process parameters in bioreactors[2,3,5],and in control of glucose and lactate concentrations [1,4]

## Background
One of the objectives in the production of therapeutic proteins is to maximize cell culture productivity. A key limitation in the cell culture processes is the accumulation of lactate. Since excess glucose can lead to increased lactate production through increased glycolytic activity, online Raman spectroscopic measurements of glucose and lactate concentrations enable both monitoring operational problems and tighter control of glucose and lactate concentrations [1].

The Raman spectra are correlated with glucose, glutamine, glutamate, lactate, ammonium, viable cell density (VCD), and total cell density (TCD) using multivariate statistical analysis. This case study will educate the trainees on how to use in-situ Raman spectroscopy data and build cell culture models to support real-time bioprocess monitoring, prediction, and control.

## Objective
Apply Raman spectroscopy for in-line monitoring for mammalian cell culture bioreactorsto ensure the specifications on CQAs are satisfied.

## PAT Methodology
1. **Raman spectroscopic analysis and cell culture process modeling**. Mechanistic and hybrid models can be used to characterize interactions between cells and extracellular environments, and improve understanding, support reliable prediction, and guide process control. The specific steps are
    1. Calibrate the models for nutrients, metabolites, and cell densities, with a focus on glucose and lactate concentrations and viable cell density[2,3,5].
    2. Preprocess Raman spectra by applying standard normalization and the second derivative to distinguish peaks and eliminate baseline shift [2,3,5].
    3. Process Raman spectroscopy by applying multivariate statistical (aka chemometric) models to provide in-situ process state measurement and quantify the measurement errors.
    4. Construct a cell culture process mechanistic model characterizing the dynamic evolution of state and inherent stochasticity.
2. Model-based bioprocess prediction: Use the bioprocess mechanistic model to provide insightful predictions on how the effect of inputs (e.g., glucose and lactate concentrations and feeding strategy) affect cell culture process dynamics.
3. Model Evaluation: The model performance is evaluated by using cross-validation method and prediction on an independent test set. The statistical metrics will include root mean square errors and square correlation coefficients (R2) for calibration, cross-validation, and prediction.


## Example Experimental Design
The bioreactors are operated in fed-batch mode using proprietary basal and feed medium formulations. Commercially available chemically defined basal medium (CD) is used in both cell expansion and production stages. A fed-batch operation has nutrient feeds entering on predetermined schedule. Medium, feeds, and stock solutions are known, and all cultures are harvested at the same time. The bioreactors were operated at 37°C. The pH and DO are maintained at 6.94±0.04 and 50% of air saturation, respectively. The bioreactor pH
is controlled by addition of CO2 gas to decrease pH and addition of 3N NaOH base to increase pH as needed. Inoculum conditions and feeding rate of substrate are controlled in order to provide conditions in the bioreactor to achieve a wider range in nutrients, metabolites, and cell densities than would be observed under normal operating conditions, which in turn allows for generation of more robust calibration models. Experiments are run to generate the at-line measurements of glucose, glutamine, glutamate, lactate, ammonium, viable cell density (VCD), and total cell density (TCD). Cell densities and viability are measured at-line. Culture samples are also analyzed at-line to monitor pH, pCO2, pO2, glucose, glutamine, glutamate, lactate, and/or ammonium.


## References
1. Thomas E. Matthews, Brandon N. Berry, John Smelko, Justin Moretto, Brandon Moore, Kelly Wiltberger (Biogen). Closed loop control of lactate concentration in Mammalian cell culture by Raman spectroscopy leads to improve cell density, viability, and biopharmaceutical protein production. Biotechnology and Bioengineering, Vol. 113, No. 11, Pages 2416-2424, November 2016. https://pubmed.ncbi.nlm.nih.gov/27215441/
2. Hamidreza Mehdizadeh, David Lauri, Krizia M. Karry, and Mojgan Moshgbar, Renee Procopio-Melino, and Denis Drapeau (Pfizer). Generic Raman-based calibration models enabling real-timemonitoring of cell culture bioreactors. Biotechnology Progress, Vol. 31, No. 4, Pages 1004-1013, July-August 2015. https://pubmed.ncbi.nlm.nih.gov/25825868/
3. Jessica Whelan, Stephen Craven, and Brian Glennon. In situ Raman spectroscopy for simultaneous monitoring of multiple process parameters in Mammalian cell culture bioreactors. Biotechnology Progress, Vol. 28, No. 5, Pages 1355-1362, September-October 2012.https://pubmed.ncbi.nlm.nih.gov/22740438/
4. Stephen Craven, Jessica Whelan, Brian Glennon (Applied Process Company). Glucose concentration control of a fed-batch mammalian cell bioprocess using a nonlinear model predictive controller. Journal of Process Control, Vol.24,344-357, 2014.
5. Stephen Goldrick, Alexandra Umprecht, Alison Tang, Roman Zakrzewski, Matthew Cheeks, Richard Turner, Aled Charles, Karolina Les, Martyn Hulley, Chris Spencer, and Suzanne S. Farid(nearly all employees of AstraZeneca). High-throughput Raman spectroscopy combined with innovate data analysis workflow to enhance biopharmaceutical process development. Processes, 8, 1179, 2020.
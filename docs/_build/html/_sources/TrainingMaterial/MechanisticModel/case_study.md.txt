# Case Study 1: Troubleshooting a Change in the Distribution of Galactosylated Glycoforms

## Introduction
This case study is a troubleshooting scenario in which the distribution of galactosylated glycoforms for a Chinese Hamster Ovary (CHO)-derived recombinant protein changes considerably after switching to a new chemically defined (CD) media. The case study is adapted from a real industrial problem reported by Eli Lilly and Company (McCracken et al., 2014), and the task is to track down which changes in the CD media, nutrients and metabolite resulted in the change of the distribution of galactosylated glycoforms. With this information, the composition of the new CD media can be modified to shift back to being near the original glycosylation profile. This case study requires constructing the relationship between a change in the raw materials (media)/metabolites and a critical quality attribute (the distribution of galactosylated glycoforms (G0F, G1F, and G2F)). This case study considers a fed-batch process and may consider two different length scales (flasks and bioreactors).
The objective of this case study is to improve the trainees’ understanding of process variation in critical quality attributes (CQAs) and how to identify the critical process parameters (CPPs), infer the causal relationships between CPPs and CQAs, and use process analytical technology (PAT) tools to find the root cause for the CQA variation. The trainees are expected to quantitatively analyze the process variation and its connection to critical process parameters, process variables, and other potential factors. In addition, they can conduct predictive analysis to guide decisions towards achieving specifications on the CQAs.

## Background
Cell culture process conditions, including media components and bioreactor operation conditions, have a profound impact on cellular metabolic state and recombinant protein quality attributes. Glycosylation is a complex biological process that can impact activity, immunogenicity, and efficacy of therapies. Considerable changes in the distribution of galactosylated glycoforms (G0F, G1F, and G2F) are observed for a CHO-derived monoclonal antibody when switching to new chemically defined (CD) media, namely, significantly lower G0F and higher G1F and G2F percentages are observed. These changes are of interest as a change in glycosylation can impact the safety and effectiveness of therapeutic proteins.

## Objective
Understand the root cause of the change and design a strategy for correcting the distribution of galactosylated glycoforms.


## Example Experimental Design at the Flask Scale to Generate Data for Troubleshooting
Flasks are seeded with approximately the same cell densities. A fed-batch operation has nutrient feeds entering on predetermined schedule. Culture cell density and viability and metabolites are monitored. Medium, feeds, and stock solutions are known and all cultures are harvested at the same time. Some experimental conditions are repeated to assess reproducibility, whereas others are not, to save time. The experiments are designed with the goal of determining which specific components in the media (such as Sialic Acid, Galactose, etc.), nutrients (glucose) and metabolite (such as ammonia) resulted in the change in the glycosylation profile. 

## Potential Root Causes
- Feed rate
- Media and metabolite
- Mn concentration
- Ammonia concentration
- Galactose feed

## Tasks
- Task 1: Use the simulator to generate a data set of glycan profiles with new CD media versus old media.
- Task 2: Illustrate the distribution difference in the digital twin interface (see figures adapted from McCracken et al., 2014)
- Task 3: Ask students to troubleshoot the problem.


# Case Study 3: Integrated Process Control of the Bioreactor and Chromatography

## Introduction
In this case, students will study two key unit operations, namely bioreactor and chromatography, in an typical biomanufacturing process. The unit operations will be highly integrated, with automatic material handling and transfer between unit operations. Students are expected to learn various control methods and trouble-shooting analysis in the integrated fed-batch bioreactor and chromatography process and then analyze its impact on the process yield and product purity. The case study is derived from real industrial problems for real-time monitoring of multiple process parameters in bioreactors[2,3,5], and in control of glucose and lactate concentrations[1,4] and chromatographic process[9,10]

## Background
One of the objectives in the production of therapeutic proteins is to maximize productivity while respecting required product quality. A key limitation in the cell culture processes is the accumulation of lactate. Since excess glucose can lead to increased lactate production through increased glycolytic activity, online measurements of glucose and lactate concentrations enable both monitoring operational problems and tighter control of glucose and lactate concentrations. In addition, the perturbation in the upstream (caused by sensor noise, raw material, or operation errors) will inevitably impact on the product quality in the downstream. Thus, it is of great interest to quantify the impact and analyze probable root causes.

In the downstream, chromatography is a critical separation step that is widely used in research and industry as an analytical or preparative tool. Separation can be achieved by differences in affinity, hydrophobicity, surface charge and molecular weight, and other factors. It involves a column which is charged with pulses of the feed solution. While travelling through the column, the more adsorptive species is retained longer by the adsorbent, thus leaving the column after the less adsorptive specie.

## Objective
Apply process control for mammalian cell culture bioreactors and downstream chromatography to improve productivity and product quality while ensuring that specifications on CQAs are satisfied. 

## Example Experimental Design
The bioreactors are operated in fed-batch mode using proprietary basal and feed medium formulations. Commercially available chemically defined basal medium (CD) is used in both cell expansion and production stages. A fed-batch operation has nutrient feeds entering on predetermined schedule. Medium, feeds, and stock solutions are known, and all cultures are harvested at the same time. The bioreactors were operated at 37°C. The pH and DO are maintained at 6.94±0.04 and 50% of air saturation, respectively. The bioreactor pH is controlled by addition of CO2 gas to decrease pH and addition of 3 N NaOH base to increase pH as needed. Inoculum conditions and feeding rate of substrate are controlled in order to provide conditions in the bioreactor to achieve a wider range in nutrients, metabolites, and cell densities than would be observed under normal operating conditions, which in turn allows for generation of more robust calibration models. Experiments are run to generate the at-line measurements of the viable cell, glucose, glutamine, lactate, ammonium and monoclonal antibody concentration. The chromatography will be simulated in the digital twin using GRM and SMA. The results of the SMA are in the form of chromatograms. Yields and purities of the desired product can be computed at various points along the simulated model chromatogram.


## Tasks
1. Sensitivity/predictive analysis of fermentation to the perturbation of CPPs. The perturbation of critical process parameters will affect the productivity and product quality. The importance of maintaining a proper level of CPPs (e.g. VCD/nitrogen concentration) has been well recognized. The sensitivity of the fermentation to these perturbations has considerable practical importance and is therefore the subject of the case study.

2. Chromatography Optimization. The inlet protein and salt concentration, entering the chromatography columns, are set from the bioreactor upstream to the separation unit. In a chromatography setup, sensors can only be placed at the entrance and exit of the column–not within its interior. This limits the applicability of classical feedback control. For example, if a load solution has an unusually high protein concentration, this will only be detected by premature breakthrough during the loading phase. Even if loading was terminated upon the detection of breakthrough, additional protein would be lost during the wash phase, compromising recovery. To mitigate these effects, open loop optimal control can be used as a control strategy. In this case, measures of the load solution concentration, obtained by UV280 absorbance or Raman spectroscopy, can be used as parameters for the chromatography model. Optimization can then be used to select the operating conditions to be used to achieve optimal process performance.

In the chronography, the column was assumed to operate in a 4-stage bind and elute mode. A total of 6 operating variables were relevant to the problem, listed in table below,

Table 1: Operating Variable of Chronography

| #   | Stage  | Time | Protein Conc. | Salt Conc.   |
|-----|--------|------|---------------|--------------|
| 1   | Load   | t_l  | $c_{p,i}^\#$  | $c_{p,l}^\#$ |
| 2   | Wash   | t_w  | 0             | $c_{s,w}$    | 
| 3   | Purge  | t_p  | 0             | $c_{s,e}$    |
| 4   | Elute  | t_e  | 0             | $c_{s,e}$    |


For a chromatographic process with given design parameters, the determination of the optimal operating regime can be posed as a combination of four competing objectives: purity, recovery, concentration, and productivity; see details in [11].


**Reference**:
1. Thomas E. Matthews, Brandon N. Berry, John Smelko, Justin Moretto, Brandon Moore, Kelly Wiltberger (Biogen). Closed loop control of lactate concentration in Mammalian cell culture by Raman spectroscopy leads to improve cell density, viability, and biopharmaceutical protein production. Biotechnology and Bioengineering, Vol. 113, No. 11, Pages 2416-2424, November 2016. https://pubmed.ncbi.nlm.nih.gov/27215441/
2. Hamidreza Mehdizadeh, David Lauri, Krizia M. Karry, and Mojgan Moshgbar, Renee Procopio-Melino, and Denis Drapeau (Pfizer). Generic Raman-based calibration models enabling real-time monitoring of cell culture bioreactors. Biotechnology Progress, Vol. 31, No. 4, Pages 1004-1013, July-August 2015. https://pubmed.ncbi.nlm.nih.gov/25825868/
3. Jessica Whelan, Stephen Craven, and Brian Glennon. In situ Raman spectroscopy for simultaneous monitoring of multiple process parameters in Mammalian cell culture bioreactors. Biotechnology Progress, Vol. 28, No. 5, Pages 1355-1362, September-October 2012. https://pubmed.ncbi.nlm.nih.gov/22740438/
4. Stephen Craven, Jessica Whelan, Brian Glennon (Applied Process Company). Glucose concentration control of a fed-batch mammalian cell bioprocess using a nonlinear model predictive controller. Journal of Process Control, Vol. 24, 344-357, 2014.
5. Stephen Goldrick, Alexandra Umprecht, Alison Tang, Roman Zakrzewski, Matthew Cheeks, Richard Turner, Aled Charles, Karolina Les, Martyn Hulley, Chris Spencer, and Suzanne S. Farid (nearly all employees of AstraZeneca). High-throughput Raman spectroscopy combined with innovate data analysis workflow to enhance biopharmaceutical process development. Processes, 8, 1179, 2020.
6. Hua Zheng, Wei Xie, Ilya O. Ryzhov, Dongming Xie (2021). Policy Optimization in Bayesian Network Hybrid Models of Biomanufacturing Processes. https://arxiv.org/pdf/2105.06543.pdf
7. Kornecki, M. and J. Strube. "Process Analytical Technology for Advanced Process Control in Biologics Manufacturing with the Aid of Macroscopic Kinetic Modeling." Bioengineering (Basel), Vol. 5, No. 1, 2018, doi:10.3390/bioengineering5010025.
8. Shekhawat, L. K., & Rathore, A. S. (2019). An overview of mechanistic modeling of liquid chromatography. Preparative Biochemistry and Biotechnology, 49(6), 623-638.
9. Iyer, H., Tapper, S., Lester, P., Wolk, B., & van Reis, R. (1999). Use of the steric mass action model in ion-exchange chromatographic process development. Journal of Chromatography A, 832(1-2), 1-9.
10. Osberghaus, A., Hepbildikler, S., Nath, S., Haindl, M., Von Lieres, E., & Hubbuch, J. (2012). Optimizing a chromatographic three component separation: A comparison of mechanistic and empiric modeling approaches. Journal of chromatography A, 1237, 86-95.
11. Lu, A. E., Paulson, J. A., Mozdzierz, N. J., Stockdale, A., Versypt, A. N. F., Love, K. R., ... & Braatz, R. D. (2015, September). Control systems technology in the advanced manufacturing of biologic drugs. In 2015 IEEE Conference on Control Applications (CCA) (pp. 1505-1515). IEEE.

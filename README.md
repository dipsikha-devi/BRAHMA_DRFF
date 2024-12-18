# BRAHMA_DRFF
This repository contains the source code for a coupled reservoir operation and one-dimensional hydraulic model integrated with an IoT-based water level acquisition system. The model was developed by the Indian Institute of Technology Guwahati (IITG) for the North Eastern Electric Power Corporation Limited (NEEPCO). BRAHMA_DRFF stands for Braided River Aid-Hydro Morphological Analyzer_Dam Release Flow Forecasting with an IoT based water level acquisition system to include pre-release status of downstream water level. The forecasted dam release, as issued by the dam operators, is given as input into the model. The model estimates flood levels in downstream flood-prone areas and provides early warnings if these levels exceed the warning level. The reservoir operating policy adopted is the Standard Operating Policy (SOP) with an objective of maximum power generation. For the hydraulic mode, the numerical scheme employed is the McCormak Predictor Corrector with Total Variation Diminishing finite difference Scheme. There is a provision of user defined pre-release strategy to release some amount before the arrival of the flood peak thereby managing the downstream floods. This inclusion of pre-release will change the existing SOP and is termed as the Adjusted Operating Policy (AOP).This model can be applied to any river system with a hydropower dam situated upstream. The BRAHMA_DRFF model has been successfully implemented in the Ranganadi Hydroelectric Power Project (RHEP) in Arunachal Pradesh, India. Three sensor locations were identified based on field surveys and reconnaissance.

<div align="center">
  <img src="https://github.com/user-attachments/assets/59708e57-9dc6-4379-9e22-1558f9305d22" alt="BRAHMA_DRFF Diagram" width="500">
</div>

## Repository Content
This repository contains the application of BRAHMA_DRFF on RHEP.This includes a bathymetry file in Excel format, which consists of three columns: the first column represents the chainage of the river, the second column indicates the river's width, and the third column contains the bathymetry data. Additionally, the repository provides a 24-hour inflow dataset for the reservoir. It also includes two other Excel files detailing the Reduced Levels of the sensors and hypothetical sensor readings.

## Outputs
The output of the model will generate the results of excel sheets for SOP and AOP including downstream flood levels, discharge, velocity,power generated and reliability of the system.

## Citation

If you use this repository or its components in your work, please cite it as follows:

Dipsikha Devi, Amit Kalita, Arup Kumar Sarma, Sanasam Ranbir Singh. "BRAHMA_DRFF" GitHub, 2024, https://github.com/dipsikha-devi/BRAHMA_DRFF.

## References
If this repository helped you, you can cite:
1. Devi, D., & Sarma, A. K. (2024). Optimal advanced release scheme based on effective forecast horizon to minimize flood at downstream of a hydroelectric project. Journal of Hydrology, 631, 130822.
2. Devi, D., & Sarma, A. K. (2023). Flow assessment downstream of a hydroelectric project in an ungauged area. Journal of Hydrologic Engineering, 28(11), 05023023.



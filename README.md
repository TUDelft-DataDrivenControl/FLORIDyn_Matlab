# FLORIDyn framework
The Flow Redirection and Induction Dynamics Model (FLORIDyn) is a dynamic wake model designed for model-based wind farm flow control.

The code is written for Matlab; for Python, see the FLORIDyn implementation in [OFF](https://github.com/TUDelft-DataDrivenControl/OFF).

## Model features
- Simulate wind farms dynamically at a low computational cost
- Estimate the power generated, added turbulence, and wake-induced losses.
- Apply heterogeneous and time-varying wind speeds and directions
- Test different modeling approaches

## Control features
- Three families of yaw angle trajectory derivation methods with each different implementations
- Different cost functions (e.g., max energy, power, shifted energy)
- Switch flow field predictions on / off or choose a transition between full and no knowledge.

## State estimation features
- Ensemble Kalman Filter flow field state estimation based on turbine measurements
- Correction of Wind speed, direction, and ambient turbulence intensity
- Different ways to calculate the Kalman Gain Matrix

## Closed-loop application examples
-  Demo files to apply closed-loop model-based wind farm flow control based on FLORIDyn in tandem with the LES [SOWFA](https://github.com/TUDelft-DataDrivenControl/SOWFA)

# Getting started
This code is folder-based. This means that a folder is referenced in the main.m, and the code draws all information from that folder. Look under „Simulations“ to find example cases. Each case consists of data (e.g., wind speeds), parameters for FLORIS and FLORIDyn, a setup.mlx, and information about the turbine placement in turbineArrayProperties.m. The closed-loop cases also have a clc_settings.mlx file and the Ensemble Kalman Filter settings.

To change settings, investigate setup.mlx, the file contains all major settings and explanations. The same holds true for the clc_stettings.mlx and EnKF_settings.m

To set up your own case, copy past an existing one and replace what is needed. If you select data input sources that don’t have a related .csv file, the code will generate one for you in the fitting format.

# Requirements
The code is based on MATLAB and tested in R2023a. The code uses the optimization toolbox, as well as the parallelization toolbox by MathWorks.

# References
Citation of the FLORIDyn model:
FLORIDyn - A dynamic and flexible framework for real-time wind farm control, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, http://doi.org/10.1088/1742-6596/2265/3/032103

Used FLORIS model:
Experimental and theoretical study of wind turbine wakes in yawed conditions, M. Bastankhah, F. Porté-Agel, 2020, http://doi.org/10.1017/jfm.2016.595

Additional references for smaller subcomponents can be found in the code or in the related publications.

## Publications that use FLORIDyn
- Ensemble-Based Flow Field Estimation Using the Dynamic Wind Farm Model FLORIDyn, M. Becker, D. Allaerts, J.W. van Wingerden, 2022, http://doi.org/10.3390/en15228589
- Wind pattern clustering of high frequent field measurements for dynamic wind farm flow control, M. Becker, D. Allaerts, J.W. van Wingerden, 2024, http://doi.org/10.1088/1742-6596/2767/3/032028 
- Sensitivity analysis and Bayesian calibration of a dynamic wind farm control model: FLORIDyn, V.V. Dighe, M. Becker, T. Göçmen, B. Sanderse, J.W. van Wingerden, 2022, http://doi.org/10.1088/1742-6596/2265/2/022062
- Time-shifted cost function design for more efficient dynamic wind farm flow control, M. Becker, D. Allaerts, J.W. van Wingerden, 2024, http://doi.org/10.1109/CCTA60707.2024.10666535
- Suitability of Dynamic Wake Models for AEP Estimation: A Wind Farm-Scale Validation Study, M. Van der Straeten, http://resolver.tudelft.nl/uuid:f35617a2-2409-439b-8bc2-6334b807ce1f 
- Scaling DMD modes for modeling Dynamic Induction Control wakes in various wind speeds, J. Gutknecht, M. Becker, C. Muscari, T. Lutz, J.W. van Wingerden, 2023, http://doi.org/10.1109/CCTA54093.2023.10252400
- Model predictive control of wakes for wind farm power tracking, A. Sterle, C.A. Hans, J. Raisch, 2024, http://doi.org/10.1088/1742-6596/2767/3/032005

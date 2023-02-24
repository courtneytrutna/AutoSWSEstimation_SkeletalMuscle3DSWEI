# AutoSWSEstimation_SkeletalMuscle3DSWEI
Automatic algorithms for shear wave speed (SWS) estimation from rotational 3D shear wave elasticity imaging (SWEI) data in skeletal muscle.

Information about the performance and validation of these algorithms, along with lateral range optimization, will be available in a manuscript currently in preparation.

For further information, please contact Courtney Trutna Paley (courtney.trutna.paley@duke.edu)  

# Algorithm Overview 
Several SWS estimation algorithms are included in this repository. All algorithms take in SWEI particle motion data as a function of time, lateral/radial position, and acquisition rotation angle. Primary outputs are the estimate of SWS parallel and perpendicular to the muscle fibers. The rotation of the muscle fibers relative to nominal 0 degrees is also determined, and for some algorithms, a data quality metric is also available.    

The algorithms included are organized around three factors:
* SWS Estimation Method: 
  * Radon Sum
  * Time-To-Peak
  * Cross Correlation
* Lateral Range Selection: 
  * Fixed: Single lateral range for all rotation angles
  * Angle Informed: lateral range changes based on rotational distance from muscle fiber direction
  * Dynamic: lateral range is dynamically set using quality metrics
* Methods for Improving Robustness
  * 3D-SWS Estimation: Detecting wave trajectory and SWS from all angles simultaneously, forcing SWS ellipse expected in TI materials
  * Multiwave Detection: Detecting multiple waves when present in data, and using a priori information to determine which one should be used for further analysis

Most iterations of these factors are included.

# Running SWS Estimation
All implemented SWS estimation methods are controlled through `RunAuto3DSWSFromPlanes.m`. Instructions for fitting are set using a structure called `setupdataprocessing`. Specifically, `setupdataprocessing.SWSEstimationMethod` is a string that determines which method is used, and `setupdataprocessing.SWSEstimationParams` is a string that determines additional, method-dependent parameters for SWS estimation. 

Parameter settings for each method are shown in `runexamples.m`. 

# Data Setup
Input data is assumed to be organized in a certain way and have a certain format. A small dataset is included in this repository as an example. The script `runexamples.m` can be used to run SWS estimation over this sample data.

The expected organization and format is as follows: 
Each individual 3D-SWEI acquisition (set of all rotation angles) should have it's own directory. 
Within this directory should be a folder called `PlanesForAnalysis`. 
Within this folder, data should be saved in a .mat file labeled `dataplanes***.mat`, where ***'s can be used to designated different preprocessing methods
This file should have a variable `out`, a 1xn structure, where n is the number of acquisition angles. 
The fields of this structure are `latmm` (for lateral position in mm), `tms` (time in ms), and `plane` (particle motion data, with size = length(latmm) x length(tms))

For Angle Informed lateral ranges, the acquisition directory should have a second folder named `FiberAngleEstimate`, containing a .mat file. 
The name of this .mat file can be set in `setupdataprocessing`, but the file should include a variable `rot_angles`, which is a single value indicating the estimated difference from the fiber rotation angle to the nominal zero degrees.
Note this is fiber rotational angle is an a priori estimate, and not the output of SWEI processing.  
  
Output data will be saved within this directory in a new folder named `SWS3DEstimation`. 
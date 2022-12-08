# TCO-estimator-current-future-HCVs
**Total ownership-cost estimation of operating current and future high-capacity vehicles in current infrastructure**

This repository contains:
1. A simulink model of an Internal combustion high-capacity vehicle with the environment.
2. A simulink model of an Electric high-capacity vehicle with the environment.
3. A simulink model of a Hybrid combination with an Internal combustion tractor and an Electric semitrailer with the environment.
4. The Matlab scripts required to process the GPS data, run the appropriate model and plot the relevant results.
5. The TCO calculator in excel.


To simulate each of the models and get a TCO estimate, follow the respective instructions:

**Step 1: Parameterizing a new vehicle**

(a) Internal combustion
- Copy one of the existing Internal combustion vehicle definition folder and rename it according to your vehicle.
- Utilize the data available in the report appendix to choose and parameterize the files for your choice of vehicle.
- The files that need to be parameterized are named hereby : 
- Save all the parameterization files in your vehicle folder.

(b) Electric
- Copy the existing Electric vehicle definition folder and rename it according to your vehicle.
- Utilize the data available from your vehicle manufacturer to parameterize your choice of vehicle.
- The files that need to be parameterized are named hereby : 
- Save all the parameterization files in your vehicle folder.

(c) Hybrid combination



**Step 2: Getting the specific GPS route**

- Open google maps from 'https://www.google.co.in/maps'.
- Select the route between your start and end destinations; and get the directions in the driving mode.
- Manually identify and note the slow sections and stop locations on the route.
- Copy the complete webpage url.
- Open GPS visualizer from 'https://www.gpsvisualizer.com/convert_input'.
- Enter the copied url in the blank space.
- Paste your GOOGLE API key in the space provided.
- Select slope and distance along with best available elevation data.
- Save file as comma separated values.
- Open and edit obtained csv file to read data in Matlab by removing the first 4 rows of data.
- Rename file and save it in this Matlab directory.


**Step 3: Preparing an Operating Cycle for the route**

- Copy one of the existing '...OC_description.xlsx' file and rename it according to your route.
- Parameterize the 'Slow' sheet with all the slow sections representing Urban areas, traffic and intersections.
- Parameterize the 'Stop' sheet with specifically planned stop points on the route.
- Save it in the same Matlab directory.


**Step 4: Running the simulation to find the energy consumption on the GPS route**

(a) Internal Combustion
- Open the Matlab script file 'ICE_Script_GPS.m'.
- Specify your vehicle parameterization folder location in all the 'xlsread' commands, thus updating to your parameterization files.
- Replace the route and operating cycle filenames in the script with your saved route.
- The files that need to be read in Matlab are named hereby :
- Run the script, automatically simulating 'Refined_Model_45_Autofactor_GPS.slx'
- Get the total energy consumption from the resulting plots.

(b) Electric

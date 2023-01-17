# TCO-estimator-current-future-HCVs
**Total ownership-cost estimation of operating current and future high-capacity vehicles in current infrastructure**

**Major Project - Kshitij Soni (659097)**

This repository contains:
1. A simulink model of an Internal combustion high-capacity vehicle with the environment.
2. A simulink model of an Electric high-capacity vehicle with the environment.
3. A simulink model of a Hybrid combination with an Internal combustion tractor and an Electric semitrailer with the environment.
4. The Matlab scripts required to process the GPS data, run the appropriate model and plot the relevant results.
5. The TCO calculator in excel.


To simulate each of the models and get a TCO estimate, follow the respective instructions:

**Step 1: Parameterizing a new vehicle**

(a) Internal combustion
- Copy the existing Internal combustion vehicle definition folder 'Vehicle 2 - Class 5 TnT' and rename it according to your vehicle.
- Utilize the data available in the **report appendix** to choose and parameterize the files for your choice of vehicle.
- The files that need to be parameterized are named hereby : 

      325kW (2).xlsx / 175kW (2).xlsx             Engine BSFC Map
      325kW.xlsx / 175kW.xlsx                     Engine Torque Map
      All Gears.xlsx                              All Gears Loss Maps
      AxlLoss.xlsx                                Axle Loss Maps
      Basic_GBXDef.xlsx                           Basic gearbox ratio definition
      FC_Map.xlsx                                 Fuel Consumption Map
      Input_Definition_File.xlsx                  Vehicle and environment Parameterization
      Throttle Torque.xlsx                        Throttle Map
      
- Save all the parameterization files in your vehicle folder.

(b) Electric
- Copy the existing Electric vehicle definition folder 'Vehicle 5 - Class 3 BEV loaded' and rename it according to your vehicle.
- Utilize the data available from your vehicle manufacturer to parameterize the file for your choice of vehicle.
- The files that need to be parameterized are named hereby : 

      Input_Definition_File.xlsx                  Vehicle and environment Parameterization
      
- Save all the parameterization files in your vehicle folder.

(c) Hybrid combination
- Copy the existing ICE tractor and E-trailer combination definition folder 'Vehicle 6 - Class 5 tractorEtrailer' and rename it according to your vehicle.
- Utilize the data available from your **tractor manufacturer and the report appendix** to parameterize the file for your choice of ICE tractor.
- Utilize the data available from your E-trailer manufacturer to parameterize the file for your choice of e-trailer.
- The files that need to be parameterized are named hereby : 

      325kW (2).xlsx / 175kW (2).xlsx             Engine BSFC Map
      325kW.xlsx / 175kW.xlsx                     Engine Torque Map
      All Gears.xlsx                              All Gears Loss Maps
      AxlLoss.xlsx                                Axle Loss Maps
      Basic_GBXDef.xlsx                           Basic gearbox ratio definition
      FC_Map.xlsx                                 Fuel Consumption Map
      Input_Definition_File.xlsx                  Vehicle and environment Parameterization
      Throttle Torque.xlsx                        Throttle Map
      Input_Definition_File_Etrailer.xlsx         Trailer drivetrain and environment Parameterization
      
- Save all the parameterization files in your vehicle combination folder.


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


**Step 4: Set the dirctory in the first line of code to the location of the folder - Total Simulink Files**


**Step 5: Running the simulation to find the energy consumption on the GPS route**

(a) Internal Combustion
- Open the Matlab script file 'ICE_Script_GPS.m'.
- Specify your vehicle parameterization folder location in all the 'xlsread' commands, thus updating to your parameterization files.
- Replace the route and operating cycle filenames in the script with your saved route.
- Run the script section by section, automatically simulating 'Refined_Model_45_Autofactor_GPS.slx'
- Get the total fuel consumption in litres from the resulting plots.

(b) Electric
- Open the Matlab script file 'BEV_Script_GPS.m'.
- Specify your vehicle parameterization folder location in all the 'xlsread' commands, thus updating to your parameterization files.
- Replace the route and operating cycle filenames in the script with your saved route.
- Run the scriptsection by section, automatically simulating 'Refined_Model_43_Elec5_ch.slx'
- Get the end battery status of charge, thus determining the total battery capacity used in kWh from the resulting plots.

(c) Hybrid Combination
- Open the Matlab script file 'ICEBEV_Script_GPS.m'.
- Specify your combination parameterization folder location in all the 'xlsread' commands for both tractor and e-trailer parameterization.
- Replace the route and operating cycle filenames in the script with your saved route.
- Run the script section by section, automatically simulating 'Refined_Model_45_Autofactor_GPS.slx' and 'Refined_Model_45_Etrailer.slx'.
- Get the end fuel consumption and battery status of charge, thus determining the total energy used by both drivetrains.

**Step 6: TCO estimation**

- Open the excel document 'TCO calculator.xlsx' and choose the cost sheet for your respective vehicle type.
- The blocks highlighted as Inputs need to be parameterized according to the economic factors and variable costs affecting the operation of the vehicle on the simulated route.
- More details about the acquired cost data are available in the report.
- From the simulation results, fill the trip length in kms and the energy consumed in Litres or kWh, according to the vehicle.
- The calculator will provide a Total ownership-cost estimate in terms of Euros/tonne. km i.e. Payload and distance specific costs.



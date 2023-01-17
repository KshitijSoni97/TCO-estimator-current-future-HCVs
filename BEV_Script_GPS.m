clc
clear 

%% Set the directory as the folder containing the 'Total Simulink  Model files'
cd ('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files');

%% Vehicle Definition Combustion
VehicleDefine = xlsread('Vehicle 5 - Class 3 BEV loaded\Input_Definition_File_old.xlsx');
Whlbs = VehicleDefine(1,3);                                               % Total Wheelbase in m
NoWhlFA = VehicleDefine(2,3);                                                % Number of wheels on the front axle
NoRA = VehicleDefine(3,3);                                                   % Number of Rear Axles
NoWhlRA = VehicleDefine(4,3);                                                % Number of wheels on the rear axle
DynRolRad = VehicleDefine(5,3);                                          % Dynamic rolling radius in m
rho = VehicleDefine(6,3);                                                  % Air density in kg/m^3 
Cx = VehicleDefine(7,3);                                                   % Air drag coefficient
Sx = VehicleDefine(8,3);                                                   % Vehicle frontal area in m^2
Mv = VehicleDefine(9,3);                                                 % Mass of the Vehicle in kgs 
Me = VehicleDefine(10,3);                                                     % Expendable mass/ Extra mass in kgs
Mratio = VehicleDefine(11,3);                                                % Mass split on front axle in percentage
Mt = VehicleDefine(12,3);                                                     % Mass of trailer in kgs
Mtot = Mv+Me+Mt;                                            % Total moving mass
Mvtot = Mv+Me;                                              % Mass on Truck excluding trailer
Fr = VehicleDefine(13,3);                                                  % Rolling resistance coefficient for smooth tarmac roads
g = VehicleDefine(14,3);                                                   % Gravitational acceleration in m/s^2
Jk = VehicleDefine(15,3);                                                    % Rotational inertia of one wheel in kgm^2 
InpSpeed = VehicleDefine(16,3);         % Maximum input vehicle speed in kmph
InpSpeedmps = InpSpeed/3.6;                         % Input Speed in m/s 
RxnTime = VehicleDefine(17,3);                   % Driver Reaction Time in seconds  
mu = VehicleDefine(18,3);                       % Coefficient of friction between road and tire
aMaxLat = VehicleDefine(19,3);                   % Maximum permitted lateral acceleration
FDR = VehicleDefine(20,3);
InitVel = VehicleDefine(21,3);                  % Define vehicle velocity at simulation start
B_road = VehicleDefine(22,3);                   % Define road width for speed profile smoothening
B_vehicle = VehicleDefine(23,3);                % Define width of vehicle
V = VehicleDefine(24,3);                        % Battery Nominal Voltage
C_cap = VehicleDefine(25,3);                    % Battery Charge Capacity in Ampere-hours
W_mmax = VehicleDefine(26,3);               	% Motor maximum speed in rotations per minute
P_mmax = VehicleDefine(27,3);                   % Motor peak power in kilo Watts
MaxDec = 1;                                     % Defining maximum deceleration according to literature in m/s^2
MaxAcc = 1;                                     % Defining maximum acceleration according to literature in m/s^2


% %% GPS Route Definition
% % Open maps and select route:
% maps = 'https://www.google.co.in/maps';
% web(maps)
% % Select origin and destination and get directions
% % Select mode as driving
% % copy webpage url
% 
% %% Acquiring route data
% gpsvisual = 'https://www.gpsvisualizer.com/convert_input';
% web(gpsvisual)
% % Enter the copied url in the blank space
% % Paste the GOOGLE API key in the space provided
% % Select slope and distance along with elevation data
% % Save file as comma separated values
% % Open and edit obtained csv file to read data in Matlab
% % Rename file and save it in this Matlab directory
% %____________________________________________________________________________________________________________________________________

%% Read saved file for GPS to get long-lat and elevation data and filter
R = readmatrix('GPS Routes and OC\Brno - Vienna.csv');                                 % Replace and read the saved csv file as a table
Slowtimes = readmatrix('GPS Routes and OC\Brno - Vienna - OC_description.xlsx','Sheet','Slow');     % Manually entering the urban sections and slow sections
Stoptimes = readmatrix('GPS Routes and OC\Brno - Vienna - OC_description.xlsx','Sheet','Stop');     % Manually entering the stop points
Rrep = R(:,7)==0;                                                   % Deleting the repeated waypoints
R(Rrep,:) = [];
R(isnan(R)) = 0;                                                    % Replacing NaN values with 0

%% Store values as arrays
Latitude = R(:,2);                                            
Longitude = R(:,3);
Alt = R(:,4);                                               % Altitude in Meters
Slope = R(:,5);                                             % Slope in Percentage
Dist = [1:(R(end,6)*1000)];                                 % Distance from Origin in meters
DistInt = R(:,7);

%% Convert Lat-Long to local coordinates using WGS84 Datum
[X, Y] = LatLonToMeters(Latitude,Longitude);

%% Calibrating according to the Origin Location
[Xcor, Ycor] = calibr(X,Y);

%% Interpolating to get waypoints 1m intervals
NoSteps = (round(Dist(end))+3);                             % Final covered distance with loop value correction
XYZintp = interparc(NoSteps,Xcor,Ycor,Alt);

%% Getting the slope between the interpolated points       
[IntpSlope] = slope(XYZintp(:,3));
IntpSlope = round(IntpSlope,1);

%% Calculating the heading angle at every m
HeadAng = HeadingAngleFcn(XYZintp);                      
HeadAngRnd = round(HeadAng);
for i = 1:length(HeadAngRnd)
    if HeadAngRnd(i) > 90
        HeadAngRnd(i) = 180 - HeadAngRnd(i);
    elseif HeadAngRnd(i) < -90
        HeadAngRnd(i) = 180 + HeadAngRnd(i);
    end
end

%% Generating route segments based on heading angle
[SegNo, SegLeninm, CurRadinm, CovDist] = SegRoute(HeadAngRnd);     % Creating route segments from the heading angle
RoadProfile = [SegNo, SegLeninm, CurRadinm, CovDist];

%% Generating a target speed in m/s for the 1m route segments based on steering inputs
SegSPLatinmps = (sqrt(aMaxLat*g*CurRadinm));                       % Using defined max lateral acceleration to determine speed on the route
SegSPLatinmps(SegSPLatinmps>InpSpeedmps) = InpSpeedmps;            % Setting the user defined maximum speed in the speed profile
MaxSpeedProfile = repelem(SegSPLatinmps,SegLeninm);                % Getting the complete speed profile according to lateral acceleration at every 1m interval for Max Lat Acc.                                                           
MSPO = MaxSpeedProfile;

%% Filtering the TSP according to a factor between 30 and 80 meters 
% This is done to reduce the speed fluctuations caused by rogue values    
s = divisors(length(MaxSpeedProfile));                             % Finding a factor to smoothen the TSP
for i = 1: length(s)
    if s(i) >= 10 && s(i) <= 150
        sect = max(s(i));
    end
end
for i = 2:length(SegSPLatinmps)-2                                  % Removing extreme reductions as rogue data if length of reduction between stable speeds is less than 250 m (Take-off distance average)
    if SegSPLatinmps(i-1) >= 13 && SegSPLatinmps(i+1) >= 13 && SegLeninm(i) <= 3*sect && SegSPLatinmps(i) < SegSPLatinmps(i+1)
        SegSPLatinmps(i) = SegSPLatinmps(i+1);
    end
end
FMSP = repelem(SegSPLatinmps,SegLeninm);
FMSP1 = FMSP;
for j = 1:sect:length(FMSP)
    i = j;
    FMSP1(j:j+(sect-1)) = min(FMSP(i:i+(sect-1)));
end
IntpSlope(length(FMSP1)+1:end) = []; 

%% Inserting stop points and slow sections
Stoptimes(isnan(Stoptimes)) = length(FMSP1);                    % Automatically entering the end distance
Slowtimes(isnan(Slowtimes)) = length(FMSP1);

for i = 1:length(Stoptimes(:,1))                                % Adding the stop times
    for j = 1:length(FMSP1)
        if j == Stoptimes(i,1)
            FMSP1(i) = 0;
        end
    end
end
for i = 1:length(Slowtimes(:,1))                                % Adding the slow sections
    j = Slowtimes(i,1);
    k = Slowtimes(i,2);
    SP = FMSP1(j:k);
    SP(SP > (Slowtimes(i,3)/3.6)) = (Slowtimes(i,3)/3.6);
    FMSP1(j:k) = SP;
end  
Dist = [1:length(FMSP1)];
%____________________________________________________________________________________________________________________________________

%% Running the Simulation
out = sim('Refined_Model_43_Elec5_ch.slx',30000);           
%____________________________________________________________________________________________________________________________________
 
% %% Plot Route only for GPS maps
% figure('Name','Route and Interpolated Map')                 % Route in X-Y coordinates
% title('Interpolated Map')
% plot(XYZintp(:,1),XYZintp(:,2),'k.','Linewidth',1.5)
% grid on
% axis equal
% xlabel('Interpolated X (m)')
% ylabel('Interpolated Y (m)')
% 
% figure('Name','Gradient Profile')                           % Gradient Profile
% title('Interpolated Altitude')
% plot(XYZintp(:,3),'k.','Linewidth',1.5)
% grid on
% xlabel('Interpolated Waypoints (m)')
% ylabel('Interpolated Altitude (m)')
% 
% figure('Name','Heading change/steering requirement')         % Route heading angle change
% plot(HeadAngRnd,'b','Linewidth',1.5)
% grid on
% xlabel('meters')
% ylabel('degrees')
% 
% figure('Name','Operating Cycle')                           % Speed profile plot
% hold on
% title('Generated Operating Cycle')
% plot(MaxSpeedProfile,'k','Linewidth',1.5)
% hold on
% ylim([0 17])
% xlabel('Displacement [m]','FontSize', 15)
% ylabel('Target Speed [m/s]','FontSize', 15)
% hold off
% %____________________________________________________________________________________________________________________________________

%% Behaviour Plots

figure('Name','Operating Cycle')
hold on
title('Operating Cycle','FontSize', 20)
plot(FMSP1,'b','Linewidth',1.5)
grid on
xlabel('Displacement [m]','FontSize', 15)
ylabel('Target Speed [m/s]','FontSize', 15) 


%% Vehicle Behaviour Plots

figure('Name','Vehicle Velocity')
hold on
title('Vehicle Velocity [m/s]','FontSize', 20)
plot(out.dist.data, out.vel.data,'b','Linewidth',1.5)
hold on
plot(out.dist.data, out.targ.data,'k','Linewidth',1.5)
grid on
set(gca,'FontSize',20)
legend('Vehicle Speed [m/s]','Target Speed [m/s]','FontSize', 20)
xlabel('Distance [m]','FontSize', 20)
ylabel('Speed [m/s]','FontSize', 20)

figure('Name','Energy Consumption [kWh]')
hold on
title('Battery SOC [%]','FontSize', 20)
plot(out.dist.data,out.SOC.data ,'r','Linewidth',1.5)
grid on
set(gca,'FontSize',20)
xlabel('Distance [m]','FontSize', 20)
ylabel('SOC %','FontSize', 20)


% ____________________________________________________________________________________________________________________________________


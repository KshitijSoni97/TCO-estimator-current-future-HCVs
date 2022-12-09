clc
clear 

%% Vehicle Definition Combustion
VehicleDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 5 - Class 3 BEV loaded\Input_Definition_File.xlsx');
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


%% VECTO Routes 1.1:  Target and Gradient Profiles from VECTO
% Comment section out for using BASIC or GPS test route
clc
ValDat = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\VECTO Results\Class 5 Engg\Class 5 all routes_LongHaul_1Hz.csv');
ValTime = ValDat(:,1);
ValDist = ValDat(:,3);
ValVel = ValDat(:,4);
ValTarg = ValDat(:,5);
ValGrad = ValDat(:,7);
ValWeng = ValDat(:,9);
ValTqeng = ValDat(:,10);
ValFC = ValDat(:,47);
ValG = ValDat(:,8);
ValFCcor = ValDat(:,49);

% Profile = [SegNo CovDist MaxSpeedProfile];
[Dist,MaxSpeedProfile,NoSteps,IntpSlope] = VectoSPpm(ValDist,ValTarg,ValGrad);
SegmentNo = groupConsec(MaxSpeedProfile);
j = unique(SegmentNo);
G = groupcounts(SegmentNo);
for i = 1:max(SegmentNo)
    E(i) = sum(G(1:i));
    SegNo(i,1) = i;
    idx = (E(i));
    CovDist(i,1) = Dist(idx,1);
    SegSPLatinmps(i,1) = MaxSpeedProfile(idx,1);
end
Profile = [SegNo CovDist SegSPLatinmps];
l = find(SegSPLatinmps == 0);
Stoptimes(:,1) = CovDist(l);
Stoptimes(:,2) = [ 21 ];                   %% Input all stop durations
% 21, 19, 8, 37, 19, 7, 260, 313, 11, 31, 5
%____________________________________________________________________________________________________________________________________

% %% BASIC Route 1:  Define Straight and no slope Test Route, with start stop at 300000
% % Comment section out for using VECTO or GPS test route
% clc
% NoSteps = 300000;                                                 % Enter total distance in m
% A = zeros(NoSteps,1);                                         % Straight road X = 0
% B = [0:1:NoSteps-1];                                              % Straight road Y = 0:Nosteps
% B = B';
% IntpSlope = A;                                                  % No Slope
% XYZIntp = [A B IntpSlope];
% HeadAngRnd = A;
% % RoadProfile = [SegNo, SegLeninm, CurRadinm, CovDist];
% RoadProfile = [1 1 inf 1; 2 NoSteps-2 inf NoSteps-1; 5  1 0 NoSteps];
% SegNo = RoadProfile(:,1);
% SegLeninm = RoadProfile(:,2);
% CurRadinm = RoadProfile(:,3);                                   % Target speed here
% CovDist = RoadProfile(:,4);
% l = find(CurRadinm == 0);
% Stoptimes(:,1) = CovDist(l);
% Stoptimes(:,2) = [10];
% Dist = [1:CovDist(end)];
% SegSPLatinmps = (sqrt(aMaxLat*g*CurRadinm));
% SegSPLatinmps(SegSPLatinmps>InpSpeedmps) = InpSpeedmps;                  
% MaxSpeedProfile = repelem(SegSPLatinmps,SegLeninm);      
% % %%
out = sim('Refined_Model_43_Elec4_regen.slx',30000);           
% %____________________________________________________________________________________________________________________________________
 

% %% GPS Route Definition
% % Comment Section out to use BASIC or VECTO test routes
% % Open maps and select route:
% maps = 'https://www.google.co.in/maps';
% web(maps)
% % Select origin and destination and get directions
% % Select mode as driving
% % copy webpage url
% 
% %% Acquiring route data
% % Comment Section out to use BASIC or VECTO test routes
% gpsvisual = 'https://www.gpsvisualizer.com/convert_input';
% web(gpsvisual)
% % Enter the copied url in the blank space
% % Paste the GOOGLE API key in the space provided
% % Select slope and ditance along with elevation data
% % Save file as comma separated values
% % Open and edit obtained csv file to read data in Matlab
% % Rename file and save it in this Matlab directory
% %____________________________________________________________________________________________________________________________________
% 
% %% Read saved file for GPS to get long-lat and elevation data and filter
% % % Comment following sections out to use BASIC or VECTO test routes
% R = readmatrix('Arnhem to Gildemeestersplein.csv');                   % Replace and read the saved csv file as a table
% Rrep = R(:,7)==0;                                            % Deleting the repeated waypoints
% R(Rrep,:) = [];
% R(isnan(R)) = 0;                                            % Replacing NaN values with 0
% 
% %% Store values as arrays
% Latitude = R(:,2);                                            
% Longitude = R(:,3);
% Alt = R(:,4);                                                % Altitude in Meters
% Slope = R(:,5);                                              % Slope in Percentage
% Dist = [1:(R(end,6)*1000)];                                               % Distance from Origin in meters
% DistInt = R(:,7);
% 
% %% Convert Lat-Long to local coordinates
% [X, Y] = LatLonToMeters(Latitude,Longitude);
% 
% %% Calibrating according to the Origin Location
% [Xcor, Ycor] = calibr(X,Y);
% 
% %% Interpolating to get waypoints at 1m intervals
% NoSteps = round(Dist(end));
% XYZintp = interparc(NoSteps,Xcor,Ycor,Alt);
% 
% %% Getting the slope between 1m interpolated points       
% [IntpSlope] = slope(XYZintp(:,3));
% IntpSlope = round(IntpSlope,1);
% 
% %% Calculating the heading angle at every m
% HeadAng = HeadingAngleFcn(XYZintp);                      
% HeadAngRnd = round(HeadAng);
% for i = 1:length(HeadAngRnd)
%     if HeadAngRnd(i) > 90
%         HeadAngRnd(i) = 180 - HeadAngRnd(i);
%     elseif HeadAngRnd(i) < -90
%         HeadAngRnd(i) = 180 + HeadAngRnd(i);
%     end
% end
% 
% %% Generating route segments based on heading angle
% [SegNo, SegLeninm, CurRadinm, CovDist] = SegRoute(HeadAngRnd);           % Creating route segments from the heading angle
% RoadProfile = [SegNo, SegLeninm, CurRadinm, CovDist];
% 
% %% Generating a target speed in m/s for the 1m route segments based on steering inputs
% SegSPLatinmps = (sqrt(aMaxLat*g*CurRadinm));                       % Using defined max lateral acceleration to determine speed on the route
% SegSPLatinmps(SegSPLatinmps>InpSpeedmps) = InpSpeedmps;                  % Setting the user defined maximum speed in the speed profile
% MaxSpeedProfile = repelem(SegSPLatinmps,SegLeninm);                    % Getting the complete speed profile according to lateral acceleration at every 1m interval for Max Lat Acc.                                                           
% 
% %% Positioning the vehicle on the path
% Dist = [floor(Whlbs):(R(end,6)*1000)];          % Start point for front axle
% cnta = length(Dist);
% cntb = length(IntpSlope);
% cntc = length(MaxSpeedProfile);
% IntpSlope = IntpSlope(cntb-cnta+1:end);                 % Offsetting data according to the start point
% MaxSpeedProfile = MaxSpeedProfile(cntc-cnta+1:end);     % Offsetting data according to the start point
% 
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

figure('Name','Operating Cycle: VECTO')
hold on
title('Operating Cycle: VECTO','FontSize', 20)
plot(MaxSpeedProfile,'b','Linewidth',1.5)
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


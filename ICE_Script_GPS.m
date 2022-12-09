clc
clear 
%%
VehicleDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\Input_Definition_File.xlsx');
Whlbs = VehicleDefine(1,3);                                 % Total Wheelbase in m
NoWhl = VehicleDefine(2,3);                                 % Number of wheels on the front axle
NoAxl = VehicleDefine(3,3);                                 % Number of Rear Axles
NoWhlDri = VehicleDefine(4,3);                              % Number of wheels on the rear axle
DynRolRad = VehicleDefine(5,3);                             % Dynamic rolling radius in m
rho = VehicleDefine(6,3);                                   % Air density in kg/m^3 
Cx = VehicleDefine(7,3);                                    % Air drag coefficient
Sx = VehicleDefine(8,3);                                    % Vehicle frontal area in m^2
Mv = VehicleDefine(9,3);                                    % Mass of the Vehicle in kgs 
Me = VehicleDefine(10,3);                                   % Expendable mass/ Extra mass in kgs
Mratio = VehicleDefine(11,3);                               % Mass split on front axle in percentage
Mt = VehicleDefine(12,3);                                   % Mass of trailer in kgs
Mtot = Mv+Me+Mt;                                            % Total moving mass
Mvtot = Mv+Me;                                              % Mass on Truck excluding trailer
Fr = VehicleDefine(13,3);                                   % Rolling resistance coefficient for smooth tarmac roads
g = VehicleDefine(14,3);                                    % Gravitational acceleration in m/s^2
Jk = VehicleDefine(15,3);                                   % Rotational inertia of one wheel in kgm^2 
Jm = VehicleDefine(16,3);                                   % Total interia of Driveline components in kgm^2
InpSpeed = VehicleDefine(17,3);                             % Maximum input vehicle speed in kmph
InpSpeedmps = InpSpeed/3.6;                                 % Input Speed in m/s 
RxnTime = VehicleDefine(18,3);                              % Driver Reaction Time in seconds  
mu = VehicleDefine(19,3);                                   % Coefficient of friction between road and tire
DrLayout = VehicleDefine(20,3);                             % Driveline layout
aMaxLat = VehicleDefine(21,3);                              % Maximum permitted lateral acceleration
FDR = VehicleDefine(22,3);                                  % Final Drive Ratio
etaClutch = VehicleDefine(23,3);                            % Clutch Efficiency
InitVel = VehicleDefine(24,3);                              % Define vehicle velocity at simulation start
Je = VehicleDefine(25,3);                                   % Define engine inertia
B_road = VehicleDefine(26,3);                               % Define road width for speed profile smoothening
B_vehicle = VehicleDefine(27,3);                            % Define width of vehicle
PTOF = VehicleDefine(28,3);                                 % Define high FC factor
CF = VehicleDefine(29,3);                                   % Define high FC factor
dt = 0.1;                                                   % Simulation sample time
MaxDec = 1;                                                 % Defining maximum deceleration according to literature in m/s^2
MaxAcc = 1;                                                 % Defining maximum acceleration according to literature in m/s^2

%% ICE  - Gearbox Definition / Loss Maps
GBXDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\Basic_GBXDef.xlsx');
GNo = GBXDefine(:,1);
GRa = GBXDefine(:,2);
etaGRa = GBXDefine(:,3);

%% ICE  - All speeds loss map
% Format the gearbox loss input file into a single excel as the example [GearNo,GInprpm,GInpTor,GTorLoss]
GLoss = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\All Gears.csv');
[GearNo,GInprpm,GInpTor,GTorLoss] = GearLossMap(GLoss);

%% ICE - BSFC and Throttle Map
clc
ThrottleMap = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\Throttle Torque.xlsx');              % Upload throttle map
FCMap = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\FC_Map.xlsx');                             % Upload Fuel Consumption Map with same breakpoints
FCTableVal = FCMap(2:end,2:end);
ThrotTableVal = ThrottleMap(2:end,2:end);
ThrotLev = ThrottleMap(:,1);
ThrotLev(isnan(ThrotLev)) = [];
RpmBrkPnt = ThrottleMap(1,:);
RpmBrkPnt(isnan(RpmBrkPnt)) = [];
ThrotPos = linspace(0,1,length(ThrotLev));                  % Defining the available throttle position as an indice from 0-1
ThrotPos = ThrotPos';
MaxTorValThrot = max(ThrotTableVal,[],2);

%% ICE - Vehicle Powertrain modelling
clc
TorCurRaw= xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\325kW.csv');
NoLev = round(length(TorCurRaw),-2);                                        % No of throttle increments
EngineRpmRaw = TorCurRaw(:,1);                                              % Engine speed in rpm
FLoadTorRaw = TorCurRaw(:,2);                                               % Engine Torque in Nm
MotoringTorRaw = TorCurRaw(:,3);                                            % Engine Motoring torque in Nm
IntpTorCur = interparc(NoLev,EngineRpmRaw,FLoadTorRaw,MotoringTorRaw);      % Interpolating the torque curve to fill all data for black box model
EngineRpm = round(IntpTorCur(:,1));                                         
FLoadTor = IntpTorCur(:,2);
MotoringTor = IntpTorCur(:,3);
h = find(FLoadTor == max(FLoadTor));
EngineRpmUse1 = EngineRpm(1:h-1);
s = round(FLoadTor);
k = find(s == max(s));
FLoadTorUse2 = unique(FLoadTor(k));
idx = length(FLoadTorUse2);                                                 % Splitting the engine rpm into utilizable and too high for fuel consumption
EngineRpmUse2 = EngineRpm(h+1:h+idx);
EngineRpmUse = [EngineRpmUse1; EngineRpmUse2];
for i = 1:length(GNo)
    MaxSpinG(i) = (0.104722*max(EngineRpmUse)*DynRolRad)./(GRa(i).*FDR);    % Max possible speed at max torque in first 3 gears in m/s
    MaxWhlTorinG(i) = max(MaxTorValThrot).*GRa(i).*FDR.*(etaGRa(i)./100);   % Max available driving torque at wheels in each gear in Nm
end

MaxAdhForce = mu*(Mvtot-(Mvtot*Mratio));                                    % Driver Traction Limits for Rear wheel drive
%____________________________________________________________________________________________________________________________________

%% Axle loss definition (Using loss tables from VECTO)
AxlLoss = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 4 - Class 5 Super Eco Combi\AxlLoss.xlsx');
ALTableVal = AxlLoss(2:end,2:end);
ALInpTor = AxlLoss(:,1);
ALInpSpeed = AxlLoss (1,:);
ALInpTor(isnan(ALInpTor)) = [];
ALInpSpeed(isnan(ALInpSpeed)) = [];

%____________________________________________________________________________________________________________________________________

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
R = readmatrix('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\GPS Routes and OC\Turin - Lyon.csv');                                 % Replace and read the saved csv file as a table
Slowtimes = readmatrix('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\GPS Routes and OC\Turin - Lyon - OC_description.xlsx','Sheet','Slow');     % Manually entering the urban sections and slow sections
Stoptimes = readmatrix('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\GPS Routes and OC\Turin - Lyon - OC_description.xlsx','Sheet','Stop');     % Manually entering the stop points
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
sim('Refined_Model_45_Autofactor_GPS.slx',30000);    

%% Plot Route only for GPS maps
figure('Name','Route and Interpolated Map')                 % Route in X-Y coordinates
title('Interpolated Map')
plot(XYZintp(:,1),XYZintp(:,2),'k.','Linewidth',1.5)
grid on
axis equal
xlabel('Interpolated X (m)','FontSize', 20)
ylabel('Interpolated Y (m)','FontSize', 20)
% 
% figure('Name','Gradient Profile')                           % Gradient Profile
% title('Interpolated Altitude')
% plot(XYZintp(:,3),'k.','Linewidth',1.5)
% grid on
% xlabel('Interpolated Waypoints (m)')
% ylabel('Interpolated Altitude (m)')
%%
figure('Name','Operating Cycle')                           % Speed profile plot
hold on
title('Generated Operating Cycle','FontSize',20)
hold on
grid on
plot(FMSP1,'b','Linewidth',1.5)
xlabel('Distance [m]','FontSize', 20)
ylabel('Target Speed [m/s]','FontSize', 20)
hold off
%____________________________________________________________________________________________________________________________________

%% Vehicle Behaviour Plots
figure('Name','Vehicle Speed')
title('Vehicle Speed [m/s]')
plot(out.distval.data,out.Vmpsval.data,'b','Linewidth',1.5)
hold on
plot(out.distval.data,out.TargVval.data,'k','Linewidth',1.5)
grid on
legend('Vehicle Speed [m/s]','Target Speed [m/s]')
xlabel('Distance [m]','FontSize', 15)
ylabel('Speed [m/s]','FontSize', 15)

figure('Name','Fuel Consumption [L]')
title('Fuel Consumption Simulation model [L]','FontSize', 20)
plot(out.distval.data,out.FCtot.data,'b','Linewidth',1.5)
hold on
grid on
xlabel('Distance [m]','FontSize', 15)
ylabel('Fuel Consumption [L]','FontSize', 15)

% ____________________________________________________________________________________________________________________________________


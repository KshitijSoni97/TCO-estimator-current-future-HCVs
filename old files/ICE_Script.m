clc
clear 

%% Vehicle Definition Combustion
VehicleDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\Input_Definition_File.xlsx');
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
Jm = VehicleDefine(16,3);                                                   % Total interia of Driveline components in kgm^2
InpSpeed = VehicleDefine(17,3);         % Maximum input vehicle speed in kmph
InpSpeedmps = InpSpeed/3.6;                         % Input Speed in m/s 
RxnTime = VehicleDefine(18,3);                   % Driver Reaction Time in seconds  
mu = VehicleDefine(19,3);                       % Coefficient of friction between road and tire
DrLayout = VehicleDefine(20,3);                 % Driveline layout
aMaxLat = VehicleDefine(21,3);                   % Maximum permitted lateral acceleration
FDR = VehicleDefine(22,3);
etaClutch = VehicleDefine(23,3);
InitVel = VehicleDefine(24,3);                  % Define vehicle velocity at simulation start
Je = VehicleDefine(25,3);                       % Define engine inertia
B_road = VehicleDefine(26,3);                   % Define road width for speed profile smoothening
B_vehicle = VehicleDefine(27,3);                % Define width of vehicle
dt = 0.1;
MaxDec = 1;                                     % Defining maximum deceleration according to literature in m/s^2
MaxAcc = 1;                                     % Defining maximum acceleration according to literature in m/s^2

%% ICE  - Gearbox Definition / Loss Maps
GBXDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\Basic_GBXDef.xlsx');
GNo = GBXDefine(:,1);
GRa = GBXDefine(:,2);
etaGRa = GBXDefine(:,3);

%% ICE  - All speeds loss map
% Format the gearbox loss input file into a single excel as the example [GearNo,GInprpm,GInpTor,GTorLoss]
GLoss = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\All Gears.csv');
[GearNo,GInprpm,GInpTor,GTorLoss] = GearLossMap(GLoss);

%% ICE - BSFC and Throttle Map
clc
ThrottleMap = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\Throttle Torque.xlsx');              % Upload throttle map
FCMap = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\FC_Map.xlsx');                             % Upload Fuel Consumption Map with same breakpoints
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
TorCurRaw= xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\325kW.csv');
NoLev = round(length(TorCurRaw),-2);
EngineRpmRaw = TorCurRaw(:,1);
FLoadTorRaw = TorCurRaw(:,2);
MotoringTorRaw = TorCurRaw(:,3);
IntpTorCur = interparc(NoLev,EngineRpmRaw,FLoadTorRaw,MotoringTorRaw);
EngineRpm = round(IntpTorCur(:,1));
FLoadTor = IntpTorCur(:,2);
MotoringTor = IntpTorCur(:,3);
h = find(FLoadTor == max(FLoadTor));
EngineRpmUse1 = EngineRpm(1:h-1);
s = round(FLoadTor);
k = find(s == max(s));
FLoadTorUse2 = unique(FLoadTor(k));
idx = length(FLoadTorUse2);
EngineRpmUse2 = EngineRpm(h+1:h+idx);
EngineRpmUse = [EngineRpmUse1; EngineRpmUse2];
for i = 1:length(GNo)
    MaxSpinG(i) = (0.104722*max(EngineRpmRaw*0.8)*DynRolRad)./(GRa(i).*FDR);   % Max possible speed at max torque in first 3 gears in m/s
    MaxWhlTorinG(i) = max(MaxTorValThrot).*GRa(i).*FDR.*(etaGRa(i)./100);     % Max available driving torque at wheels in each gear in Nm
end

MaxAdhForce = mu*(Mvtot-(Mvtot*Mratio));                                    % Driver Traction Limits for Rear wheel drive
%____________________________________________________________________________________________________________________________________

%% Axle loss definition
AxlLoss = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 2 - Class 5 TnT\AxlLoss.xlsx');
ALTableVal = AxlLoss(2:end,2:end);
ALInpTor = AxlLoss(:,1);
ALInpSpeed = AxlLoss (1,:);
ALInpTor(isnan(ALInpTor)) = [];
ALInpSpeed(isnan(ALInpSpeed)) = [];

%____________________________________________________________________________________________________________________________________

%% GPS Route Definition
% Comment Section out to use BASIC or VECTO test routes
% Open maps and select route:
maps = 'https://www.google.co.in/maps';
web(maps)
% Select origin and destination and get directions
% Select mode as driving
% copy webpage url

%% Acquiring route data
% Comment Section out to use BASIC or VECTO test routes
gpsvisual = 'https://www.gpsvisualizer.com/convert_input';
web(gpsvisual)
% Enter the copied url in the blank space
% Paste the GOOGLE API key in the space provided
% Select slope and ditance along with elevation data
% Save file as comma separated values
% Open and edit obtained csv file to read data in Matlab
% Rename file and save it in this Matlab directory
%____________________________________________________________________________________________________________________________________

%% Read saved file for GPS to get long-lat and elevation data and filter
% % Comment following sections out to use BASIC or VECTO test routes
R = readmatrix('Rotterdam - arnhem.csv');                   % Replace and read the saved csv file as a table
Rrep = R(:,7)==0;                                            % Deleting the repeated waypoints
R(Rrep,:) = [];
R(isnan(R)) = 0;                                            % Replacing NaN values with 0

Stoptimes(:,1) = [ 5000, 10000, 13000, 18000, 19724 ];
Stoptimes(:,2) = [ 10, 10, 10, 10, 10 ];                   %% Input all stop durations

%% Store values as arrays
Latitude = R(:,2);                                            
Longitude = R(:,3);
Alt = R(:,4);                                                % Altitude in Meters
Slope = R(:,5);                                              % Slope in Percentage
Dist = [1:(R(end,6)*1000)];                                               % Distance from Origin in meters
DistInt = R(:,7);

%% Convert Lat-Long to local coordinates
[X, Y] = LatLonToMeters(Latitude,Longitude);

%% Calibrating according to the Origin Location
[Xcor, Ycor] = calibr(X,Y);

%% Interpolating to get waypoints at 50m intervals
NoSteps = round(Dist(end));
XYZintp = interparc(NoSteps,Xcor,Ycor,Alt);

%% Getting the slope between the 50m interpolated points       
[IntpSlope] = slope(XYZintp(:,3));
IntpSlope = round(IntpSlope,1);

%% Calculating the heading angle at every 50 m
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
[SegNo, SegLeninm, CurRadinm, CovDist] = SegRoute(HeadAngRnd);           % Creating route segments from the heading angle
RoadProfile = [SegNo, SegLeninm, CurRadinm, CovDist];

%% Generating a target speed in m/s for the 1m route segments based on steering inputs
SegSPLatinmps = (sqrt(aMaxLat*g*CurRadinm));                       % Using defined max lateral acceleration to determine speed on the route
SegSPLatinmps(SegSPLatinmps>InpSpeedmps) = InpSpeedmps;                  % Setting the user defined maximum speed in the speed profile
MaxSpeedProfile = repelem(SegSPLatinmps,SegLeninm);                    % Getting the complete speed profile according to lateral acceleration at every 1m interval for Max Lat Acc.                                                           
Dist = 1:length(MaxSpeedProfile);


%% Filtering the TSP
for i = 1:50:length(MaxSpeedProfile)
    for j = i
        FMSP(j:j+49) = min(MaxSpeedProfile(i:i+49));
    end
end
MaxSpeedProfile = FMSP;

IntpSlope(length(MaxSpeedProfile)+1:end) = [];
% %% Positioning the vehicle on the path
% Dist = [floor(Whlbs):(R(end,6)*1000)];          % Start point for front axle
% cnta = length(Dist);
% cntb = length(IntpSlope);
% cntc = length(MaxSpeedProfile);
% IntpSlope = IntpSlope(cntb-cnta+2:end);                 % Offsetting data according to the start point
% MaxSpeedProfile = MaxSpeedProfile(cntc-cnta+2:end);     % Offsetting data according to the start point

%% Inserting stop points
for i = 1:length(FMSP)
    for j = 1:length(Stoptimes)
        if i == Stoptimes(j,1)
            FMSP(i) = 0;
        end
    end
end
%____________________________________________________________________________________________________________________________________

%% Plot Route only for GPS maps
figure('Name','Route and Interpolated Map')                 % Route in X-Y coordinates
title('Interpolated Map')
plot(XYZintp(:,1),XYZintp(:,2),'k.','Linewidth',1.5)
grid on
axis equal
xlabel('Interpolated X (m)')
ylabel('Interpolated Y (m)')

figure('Name','Gradient Profile')                           % Gradient Profile
title('Interpolated Altitude')
plot(XYZintp(:,3),'k.','Linewidth',1.5)
grid on
xlabel('Interpolated Waypoints (m)')
ylabel('Interpolated Altitude (m)')

figure('Name','Operating Cycle')                           % Speed profile plot
hold on
title('Generated Operating Cycle')
plot(FMSP,'k','Linewidth',1.5)
hold on
ylim([0 17])
xlabel('Displacement [m]','FontSize', 15)
ylabel('Target Speed [m/s]','FontSize', 15)
hold off
%____________________________________________________________________________________________________________________________________

%% Vehicle Behaviour Plots
figure('Name','Vehicle position and Speed')
subplot(2,1,1)
hold on
title('Vehicle Position [m]')
plot(out.distval,'b','Linewidth',1.5)
hold on
grid on
xlabel('Time [s]','FontSize', 15)
ylabel('Displacement [m]','FontSize', 15)

subplot(2,1,2)
hold on
title('Vehicle Speed [m/s]')
plot(out.Vmpsval,'b','Linewidth',1.5)
hold on
plot(out.TargVval,'k','Linewidth',1.5)
grid on
legend('Vehicle Speed [m/s]','Target Speed [m/s]')
xlabel('Time [s]','FontSize', 15)
ylabel('Speed [m/s]','FontSize', 15)

figure('Name','Fuel Consumption [L]')
title('Fuel Consumption Simulation model [L]','FontSize', 20)
plot(out.FCtot,'b','Linewidth',1.5)
hold on
grid on
xlabel('Time [s]','FontSize', 15)
ylabel('Fuel Consumption [L]','FontSize', 15)

% ____________________________________________________________________________________________________________________________________


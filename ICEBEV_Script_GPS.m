clc
clear 

%% Vehicle Definition Combustion
VehicleDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\Input_Definition_File.xlsx');
Whlbs = VehicleDefine(1,3);                                               % Total Wheelbase in m
NoWhl = VehicleDefine(2,3);                                                % Number of wheels on the front axle
NoAxl = VehicleDefine(3,3);                                                   % Number of Rear Axles
NoWhlDri = VehicleDefine(4,3);                                                % Number of wheels on the rear axle
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
% FCH2 = VehicleDefine(28,3);                % Define high FC factor
% FCL2 = VehicleDefine(29,3);                % Define high FC factor
PTOF = 0.003255;
CF = 0.002408;
dt = 0.1;
MaxDec = 1;                                     % Defining maximum deceleration according to literature in m/s^2
MaxAcc = 1;                                     % Defining maximum acceleration according to literature in m/s^2

%% ICE  - Gearbox Definition / Loss Maps
GBXDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\Basic_GBXDef.xlsx');
GNo = GBXDefine(:,1);
GRa = GBXDefine(:,2);
etaGRa = GBXDefine(:,3);

%% ICE  - All speeds loss map
% Format the gearbox loss input file into a single excel as the example [GearNo,GInprpm,GInpTor,GTorLoss]
GLoss = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\All Gears.csv');
[GearNo,GInprpm,GInpTor,GTorLoss] = GearLossMap(GLoss);

%% ICE - BSFC and Throttle Map
clc
ThrottleMap = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\Throttle Torque.xlsx');              % Upload throttle map
FCMap = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\FC_Map.xlsx');                             % Upload Fuel Consumption Map with same breakpoints
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
TorCurRaw= xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\325kW.csv');
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

MaxAdhForce = mu*(Mvtot-(Mvtot*Mratio));                                    % Driver Traction Limits for Rear wheel drive
%____________________________________________________________________________________________________________________________________

%% Axle loss definition
AxlLoss = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\AxlLoss.xlsx');
ALTableVal = AxlLoss(2:end,2:end);
ALInpTor = AxlLoss(:,1);
ALInpSpeed = AxlLoss (1,:);
ALInpTor(isnan(ALInpTor)) = [];
ALInpSpeed(isnan(ALInpSpeed)) = [];

%____________________________________________________________________________________________________________________________________

%% Shift map
MinSp1 = ((0.377*(min(EngineRpmRaw))*DynRolRad)/GRa(1)*FDR)/3.6;
for i = 1:length(GNo)
    MaxSpinG(i) = ((0.377*max(EngineRpmRaw*0.8)*DynRolRad)./(GRa(i).*FDR))/3.6;   % Max possible speed at max torque in first 3 gears in m/s
    MaxWhlTorinG(i) = max(MaxTorValThrot).*GRa(i).*FDR.*(etaGRa(i)./100);     % Max available driving torque at wheels in each gear in Nm
    MinEngRpminG(1) = 0;
    MinEngRpminG(i+1) = (MaxSpinG(i).*FDR.*GRa(i))/(0.377.*DynRolRad);
end

%____________________________________________________________________________________________________________________________________

%% E-Trailer Definition

Etr = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 tractorEtrailer\Input_Definition_File_Etrailer.xlsx');
DynRolRadE = Etr(1,3);                          % E-trailer tyre rolling radius
Mtr = Etr(2,3);                                 % Mass of the Trailer
Mload = Etr(3,3);                               % Expendable mass/ Extra mass on vehicle
Mtrailer = Mtr + Mload;
Jk = Etr(1,3);                                  % Rotational inertia of one wheel
FDRE = Etr(1,3);                                % Final Drive Ratio of E-trailer
PrimRed = Etr(1,3);                             % Motor primary reduction
V = Etr(1,3);                                   % Battery nominal voltage
C_cap = Etr(1,3);                               % Battery charge capacity
W_mmax = Etr(1,3);                              % Motor maximum speed
P_mmax = Etr(1,3);                              % Motor peak power

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

%% Run the ICE tractor first and the E-trailer next
out1 = sim('Refined_Model_45_Autofactor_GPS.slx',20000);
TSP = out1.Vmpsval;
out2 = sim('Refined_Model_45_Etrailer.slx',20000);

%% VECTO Validation Plots

figure('Name','Operating Cycle: VECTO')
hold on
title('Operating Cycle: VECTO','FontSize', 20)
plot(FMSP1,'b','Linewidth',1.5)
grid on
set(gca,'FontSize',20)
xlabel('Displacement [m]','FontSize', 20)
ylabel('Target Speed [m/s]','FontSize', 20) 
%%
figure('Name','Vehicle velocity')
title('Velocity Comparison','FontSize', 20)
plot(out1.distval.data,out1.Vmpsval.data,'b','Linewidth',1.5)
hold on
plot(out2.dist.data,out2.vel.data,'k','Linewidth',1.5)
grid on
set(gca,'FontSize',20)
legend('Simulink model [m/s]','VECTO [m/s]','Model target [m/s]','FontSize', 12)
xlabel('Displacement [m]','FontSize', 20)
ylabel('Velocity [m/s]','FontSize', 20)

%____________________________________________________________________________________________________________________________________


%% Fuel Consumption Plots

figure('Name','Fuel Consumption [L] and Battery SOC status [%]')
title('Fuel Consumption [L] and Battery SOC status [%]')
yyaxis left
plot(out1.distval.data,out1.FCtot.data,'b','Linewidth',1.5)
xlabel('Distance [m]','FontSize', 20)
ylabel('Fuel Consumption [L]','FontSize', 20)
hold on
grid on
yyaxis right
plot(out2.dist.data,out2.SOC.data ,'r','Linewidth',1.5)
ylabel('SOC %','FontSize', 20)
grid on
set(gca,'FontSize',20)
%%
figure('Name','Driving Force Comparison')
title('Driving Force Comparison [N]')
hold on
plot(out1.distval.data,out1.T_e.data,'b','Linewidth',1.5)
hold on
plot(out2.dist.data,out2.T_m.data,'r','Linewidth',1.5)
grid on
set(gca,'FontSize',20)
legend('From engine [N]','From motor [N]','FontSize', 20)
xlabel('Distance [m]','FontSize', 20)
ylabel('Driving force on wheels [N]','FontSize', 20)
%____________________________________________________________________________________________________________________________________

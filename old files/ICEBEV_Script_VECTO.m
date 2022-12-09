clc
clear 

%% Vehicle Definition Combustion
VehicleDefine = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Vehicle 6 - Class 5 + BEV trailer loaded\Input_Definition_File.xlsx');
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
PTOF = VehicleDefine(28,3);                % Define high FC factor
CF = VehicleDefine(29,3);                % Define high FC factor
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

%% Shift map
MinSp1 = ((0.377*(min(EngineRpmRaw))*DynRolRad)/GRa(1)*FDR)/3.6;
for i = 1:length(GNo)
    MaxSpinG(i) = ((0.377*max(EngineRpmRaw*0.8)*DynRolRad)./(GRa(i).*FDR))/3.6;   % Max possible speed at max torque in first 3 gears in m/s
    MaxWhlTorinG(i) = max(MaxTorValThrot).*GRa(i).*FDR.*(etaGRa(i)./100);     % Max available driving torque at wheels in each gear in Nm
    MinEngRpminG(1) = 0;
    MinEngRpminG(i+1) = (MaxSpinG(i).*FDR.*GRa(i))/(0.377.*DynRolRad);
end


%% VECTO Route Input:  Target and Gradient Profiles from VECTO
clc
ValDat = xlsread('C:\Users\Kshitij\Desktop\Thesis docs\Code\Total Simulink  Model files\Class 5 Engg\Class 5 all routes_RegionalDelivery_1Hz.csv');
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

% Profile = [SegNo CovDist MaxSpeedProfile]
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
Stoptimes(:,2) = [ 25 ];                   %% Input all stop durations
% 21, 19, 8, 37, 19, 7, 260, 313, 11, 31, 5
%___________________________________________________________________________________________________________________________________

% %% BASIC Route 1:  Define Straight and no slope Test Route, with start stop at 300000
% % Comment section out for using VECTO or GPS test route
% clc
% InpSpeed = 54;         % Maximum input vehicle speed in kmph
% InpSpeedmps = InpSpeed/3.6;                         % Input Speed in m/s 
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
% 
% % %%
% % out = sim('Refined_Model_43_Elec.slx',300);           
% %____________________________________________________________________________________________________________________________________

%% Compare models
out = sim('Refined_Model_44_ICEtractor_BEVtrailer_2.slx',300);

%____________________________________________________________________________________________________________________________________

%% Vehicle Behaviour Plots
figure('Name','Vehicle Speed')
title('Vehicle Speed [m/s]')
plot(out.distval.data,out.Vmpsval.data,'b','Linewidth',1.5)
hold on
plot(out.distval.data,out.TargVval.data,'k','Linewidth',1.5)
grid on
xlim ([0 2950])
legend('Vehicle Speed [m/s]','Target Speed [m/s]')
xlabel('Distance [m]','FontSize', 15)
ylabel('Speed [m/s]','FontSize', 15)

%% Fuel Consumption Plots

figure('Name','Fuel Consumption [L] and Battery SOC status [%]')
title('Fuel Consumption [L] and Battery SOC status [%]')
yyaxis left
plot(out.distval.data,out.FCtot.data,'b','Linewidth',1.5)
xlabel('Distance [m]','FontSize', 20)
ylabel('Fuel Consumption [L]','FontSize', 20)
hold on
grid on
yyaxis right
plot(out.distval.data,out.SOC.data ,'r','Linewidth',1.5)
ylabel('SOC %','FontSize', 20)
grid on
xlim ([0 2950])
set(gca,'FontSize',20)

figure('Name','Driving Force Comparison')
title('Driving Force Comparison [N]')
hold on
plot(out.distval.data,out.T_e.data,'b','Linewidth',1.5)
hold on
plot(out.distval.data,out.T_m.data,'r','Linewidth',1.5)
grid on
xlim ([0 2950])
set(gca,'FontSize',20)
legend('From engine [N]','From motor [N]','FontSize', 20)
xlabel('Distance [m]','FontSize', 20)
ylabel('Driving force on wheels [N]','FontSize', 20)





%____________________________________________________________________________________________________________________________________


%% Creating the VECTO profile per meter
function [Dist,MSP,NoSteps,IntpSlope] = VectoSPpm(ValDist,ValTarg,ValGrad)
D = round(ValDist);
D(1) = [];
DP = unique(D);
NoSteps = max(D);
SP = (round(ValTarg))/3.6;
SP(1) = [];
S = round(ValGrad);
S(1) = [];
for i = 1:length(DP)                                               % Selecting minimum value at every displacement
    Dis(i,1) = DP(i);
    E = find(D == Dis(i,1));                                               % Finding all locations with a stop in the 
    if length(E)==1
        TSP(i,1) = SP(E);
        Slope(i,1) = S(E);
    else
        TSP(i,1) = min(SP(E));
        Slope(i,1) = min(S(E));
    end
end
for i = 1:NoSteps
    Dist(i,1) = i;
    if ismember(i,Dis)
        k = find(Dis(:) == Dist(i,1));
        MSP(i,1) = TSP(k,1);
        IntpSlope(i,1) = Slope(k,1);
    else
       MSP(i,1) = MSP(i-1,1);
       IntpSlope(i,1) = IntpSlope(i-1,1);
    end     
end



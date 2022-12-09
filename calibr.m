%% Calibration function
function [Xcal, Ycal] = calibr(X,Y)
    Xor = X(1,:);
    Yor = Y(1,:);
    Xcal = X-Xor;
    Ycal = Y-Yor;
end

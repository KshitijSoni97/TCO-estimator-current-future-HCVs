%% Slope Function
 function [slpm] = slope(Altitude)
    i = 1:(length(Altitude)-1);
    DelApm = Altitude((i+1),:) - Altitude(i,:);
    slpm = tan(asin(DelApm));
    slpm = rad2deg(slpm);
 end

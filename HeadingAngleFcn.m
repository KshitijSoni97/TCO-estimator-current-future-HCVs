%% Heading angle function in degrees
function [HeadAng] = HeadingAngleFcn(XYZintps)
for chgst = 1:length(XYZintps)-2                            
    j = chgst+1;
    k = j+1;
    X1 = XYZintps(chgst,1);
    X2 = XYZintps(j,1);
    X3 = XYZintps(k,1);
    Y1 = XYZintps(chgst,2);
    Y2 = XYZintps(j,2);
    Y3 = XYZintps(k,2);
    theta1 = atand([Y2 - Y1]/[X2 - X1]);                       % Change in heading direction of road in degrees
    theta2 = atand([Y3 - Y2]/[X3 - X2]);
    HeadAng(chgst) = theta1-theta2;
end
for i = 2:size(HeadAng)-1
  if HeadAng(i) >= 90 || HeadAng(i) <= -90                          % Replacing rogue values with average points
    HeadAng(i) = (HeadAng(i+1)+HeadAng(i-1))/2;
  end
end
HeadAng(isnan(HeadAng)) = 0;
HeadAng = HeadAng';
end

%% Utilizing the Gear Loss Maps
function [GearNo,Inprpm,InpTor,TorLoss] = GearLossMap(GLoss)
GearNo = unique(GLoss(:,1));
Inprpm = unique(GLoss(:,2));
InpTor = unique(GLoss(:,3));
NoTorVal = length(InpTor);
NorpmVal = length(Inprpm);
Full = GLoss(:,4);
TorLoss = reshape(Full,NoTorVal,NorpmVal,length(GearNo));
end

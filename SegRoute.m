%% Grouping the road into segments from the 1m points
function [SegNo, SegLeninm, CurRadinm, CovDist] = SegRoute(HeadAngRnd) 
SegmentNo = groupConsec(HeadAngRnd);                % Numbering the sections of constant curves
TotSeg = unique(SegmentNo);
for i = 1:length(TotSeg)-1
    A = SegmentNo == i;                                  % Gives index location of required heading angle
    B = SegmentNo == i;
    SegCurAng(i) = max(HeadAngRnd(A));                              % Angle change between 1m sections in the segment
    CurRadinm(i)  = RadiiCalc(SegCurAng(i));                     % Defining the curve radius to determine max velocity and type of environment
    SegLeninm(i) = sum(B);                                  % Segment length in meters  
    CovDist(i) = sum(SegLeninm(1:i));                       % Total covered distance
    SegNo(i) = i;                                         % Segment number
end

function labelmap = groupConsec(x)
%Takes a logical input vector and returns a label map of all the consecutive 
%ones.
% 
%       labelmap = groupConsec(x)       
%                            
%IN:                         
%                            
%    x: an input vector of numeric values.               
%               
%OUT:                        
%                            
%    labelmap: a vector of labels for the consecutive elements in x
%
%EXAMPLE:
% 
%     >> x=[20 20 7 7 7 -1 37 37 37 37 37];
%     >> labelmap = groupConsec(x)
% 
%      labelmap =
% 
%          1     1     2     2     2     3     4     4     4     4     4
 assert(isvector(x),'Input must be a vector')
 
 y=reshape(diff([0,x(:).'])~=0,size(x));
 
 labelmap=cumsum(abs(y));
end

end

function [CurRadinm] = RadiiCalc(SegCurAng)
S = (180-abs(SegCurAng))/2;
S = deg2rad(S);
CurRadinm = (5/2)./(cos(S));        % Calculating for 5m segments
end

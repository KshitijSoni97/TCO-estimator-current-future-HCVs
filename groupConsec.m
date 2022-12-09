%% Group Arrays
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

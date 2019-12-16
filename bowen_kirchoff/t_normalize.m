% Trace Normalize
% author @Bowen Guo 2014-11-19

function [out]=t_normalize(in)


[~,ng]=size(in);

out=zeros(size(in));
for ig=1:ng
    out(:,ig)=Normalize(in(:,ig));
end

return
% Kirchhoff Modeling
% @version 1 2014-05-20
% @author Bowen Guo

% ttt_s-source time table
% ttt_r-receive time table
function [d]=k_forward(refl,source,sx,sz,gx,gz,ttt_s,ttt_r,nt,dt,dx)

[nz,nx]=size(refl);

ng=numel(gx);

d=zeros(nt,ng);

dist_s=cal_dist(nz,nx,sx,sz,dx);

for ig=1:ng
    trace=zeros(nt,1);
    
    t=round((ttt_s+ttt_r(:,:,ig))/dt)+1;
    
    dist_r=cal_dist(nz,nx,gx(ig),gz(ig),dx);
    
    dist=sqrt(dist_s.*dist_r);
    
    for ix=1:nx
        for iz=1:nz
            if (t(iz,ix)<=nt)
                trace(t(iz,ix))=trace(t(iz,ix))+refl(iz,ix)/dist(iz,ix);
            end
        end
    end
    
    trace=conv2(trace,source);
    
    d(:,ig)=trace(1:nt);
end

end







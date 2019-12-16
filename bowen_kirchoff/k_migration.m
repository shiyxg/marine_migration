% Kirchhoff Migration
% @version 1 2014-05-20
% @author Bowen Guo

function [mig]=k_migration(data,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,nz,nx)


mig=zeros(nz,nx);

[nt,ng]=size(data);

dist_s=cal_dist(nz,nx,sx,sz,dx);

for ig=1:ng
    
    trace2=xcorr(data(:,ig),source);
    
    trace=trace2(nt:2*nt-1);
    
    dist_r=cal_dist(nz,nx,gx(ig),gz(ig),dx);
    
    t=round((ttt_s+ttt_r(:,:,ig))/dt)+1;
    
    dist=sqrt(dist_s.*dist_r);
    
    for iz=1:nz
        
        for ix=1:nx
            
            if (t(iz,ix)<=nt)
                
                mig(iz,ix)=mig(iz,ix)+trace(t(iz,ix))/dist(iz,ix);
            end
        end
    end
end
    
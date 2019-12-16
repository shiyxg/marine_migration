% Conjugate Gradient Least Squares Migration 
% @version 1 2014-05-24
% @author Bowen Guo


function [mig,res]=lsm_under(data,mig0,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,iter)


[nt,~]=size(data);


res=zeros(iter,1);
[nz,nx]=size(mig0);


mig=mig0;
for i=1:iter
    
    if (i==1)
       
        d=k_forward(mig,source,sx,sz,gx,gz,ttt_s,ttt_r,nt,dt,dx);
        r=d-data;
        g=k_migration(r,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,nz,nx);
        dir=-g;
    end
    
    % Calculate Step length            
    d=k_forward(dir,source,sx,sz,gx,gz,ttt_s,ttt_r,nt,dt,dx);
    alpha=sum(dir(:).*g(:))/sum(d(:).*d(:));
    
    % Update model and gradient
    mig_new=mig-alpha*dir;
    d=k_forward(mig_new,source,sx,sz,gx,gz,ttt_s,ttt_r,nt,dt,dx);
    
    r=d-data;
    res(i)=sum(r(:).^2);
    % Calculate new gradient
    
    g_new=k_migration(r,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,nz,nx);
    
    % Calculate new searching direction
    beta=(sum(g_new(:).^2)-sum(g_new(:).*g(:)))/sum(g(:).*g(:));
    dir_new=-g_new+beta*dir;
    
    % Update image, gradient and searching direction
    dir=dir_new;
    g=g_new;
    mig=mig_new;
end

end
    


    











        
        
        
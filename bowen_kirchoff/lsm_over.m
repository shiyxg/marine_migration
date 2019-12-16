function [mig,res]=lsm_over(data,mig0,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,iter)

[nz,nx]=size(mig0);

[nt,ng,ns]=size(data);

res=zeros(iter,1);

mig=mig0;

d=zeros(nt,ng,ns);
for i=1:iter
    display(['Iteration No. = ',num2str(i)]);
    if (i==1)
        parfor is=1:ns
            d(:,:,is)=k_forward(mig,source,sx(is),sz(is),gx,gz,ttt_s(:,:,is),ttt_r,nt,dt,dx);
        end
       r=d-data;
       g=zeros(nz,nx);
       parfor is=1:ns
           g=g+k_migration(r(:,:,is),source,sx(is),sz(is),gx,gz,ttt_s(:,:,is),ttt_r,dt,dx,nz,nx);
       end
       dir=-g;
    end
    
    % Calculate step length
    parfor is=1:ns
        d(:,:,is)=k_forward(dir,source,sx(is),sz(is),gx,gz,ttt_s(:,:,is),ttt_r,nt,dt,dx);
    end
    alpha=sum(dir(:).*g(:))/sum(d(:).*d(:));
    
    % Update migration image and gradient
    mig_new=mig-alpha*dir;
    parfor is=1:ns
        d(:,:,is)=k_forward(mig_new,source,sx(is),sz(is),gx,gz,ttt_s(:,:,is),ttt_r,nt,dt,dx);
    end
    
    r=d-data;
    res(i)=sum(r(:).^2);
    
    g_new=zeros(nz,nx);
    parfor is=1:ns
        g_new=g_new+k_migration(r(:,:,is),source,sx(is),sz(is),gx,gz,ttt_s(:,:,is),ttt_r,dt,dx,nz,nx);
    end
    
    % Calculate beta and new searching direction 
    
    beta=(sum(g_new(:).^2)-sum(g(:).*g_new(:)))/sum(g(:).*g(:));
    dir_new=-g_new+beta*dir;
    
    % Update image, gradient and searching direction
    dir=dir_new;
    g=g_new;
    mig=mig_new;
    
end

end
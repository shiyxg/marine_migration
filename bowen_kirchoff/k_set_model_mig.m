% building model

v0=1000;
refl=zeros(500,500);refl(248:252,:)=10;
[nz,nx]=size(refl); dx=1; x = (0:nx-1)*dx; z = (0:nz-1)*dx; 
nt=2500; dt=0.00046; t=(0:nt-1)*dt;


% gx=(0:2:(nx-1))*dx; gz=zeros(size(gx)); ng=numel(gx);
% sx=(0:20:(nx-1))*dx; 
sx=nx/2*dx;sz=0;
gx=(0:2:nx-1)*dx;gz=zeros(size(gx)); ng=numel(gx);

freq=25; 
source=ricker(freq,dt);
source=expand_source(source,nt);

ttt_s=zeros(nz,nx);
ttt_r=zeros(nz,nx,ng);
ttt_s=time_table(sx,sz,nx,nz,dx)/v0;

for i=1:ng
    ttt_r(:,:,i)=time_table(gx(i),gz(i),nx,nz,dx)/v0;
end


d=k_forward(refl,source,sx,sz,gx,gz,ttt_s,ttt_r,nt,dt,dx);
figure(3);
subplot(211);imagesc(refl);colormap(gray);
subplot(212);imagesc(d);colormap(gray);

igg=1;
%mig=k_migration(d(:,igg),source,sx,sz,gx(:,igg),gz(:,igg),ttt_s,ttt_r(:,:,igg),dt,dx,nz,nx);
mig=k_migration(d,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,nz,nx);
figure(4);
imagesc(mig);colormap(gray);

% iter=5;
% [mig_lsm,res]=lsm_over(data,mig0,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,iter);
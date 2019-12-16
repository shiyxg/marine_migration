% building model
clc
clear
v0=1e8*ones(160,200); 
refl=zeros(160,200);refl(80,50)=10; refl(80,80)=10; 
[nz,nx]=size(refl); dx=0.025; x = (0:nx-1)*dx; z = (0:nz-1)*dx; 
nt=4000; dt=5e-11; t=(0:nt-1)*dt;


% gx=(0:2:(nx-1))*dx; gz=zeros(size(gx)); ng=numel(gx);
% sx=(0:20:(nx-1))*dx; 

freq=100e6; 
source=ricker(freq,dt);
source=expand_source(source,nt);

is0 = 0.75;
group_dx = 0.05;
geophone_dx = 0.25;


sx=nx/2*dx;sz=0;
gx=(0:2:nx-1)*dx;gz=zeros(size(gx)); ng=numel(gx);


for ixsrc=1:1:51;
    ixsrc
      parfor idsrc = 1:6 %
            isx = fix(is0/dx) + fix(group_dx/dx)*(ixsrc-1) + (idsrc-1)*fix(geophone_dx/dx)
            SP = [1 isx];   % Shot Coordinates
            [Ttable]=Mray(SW,SP,DX);
            traveltimesrc(:,:,idsrc,ixsrc)=Ttable(:,:);
%           ixsrc
            imagesc(x,y,Ttable);title('Traveltime Field');colorbar;text(nx*DX+5*DX,-DX,'Time (s)');
            xlabel('X (m)');ylabel('Y (m)');            
            pause(.015)
      end
end






% sx=nx/2*dx;sz=0;
% gx=(0:2:nx-1)*dx;gz=zeros(size(gx)); ng=numel(gx);

% freq=100e6; 
% source=ricker(freq,dt);
% source=expand_source(source,nt);

ttt_r = zeros(nz,nx,ng);
ttt_s=eikonal2d(1./v0,sx,sz,dx);
  
for i=1:ng
    ttt_r(:,:,i)=eikonal2d(1./v0,gx(i),gz(i),dx);
%     imagesc(ttt_r(:,:,i))
end
 
% ttt_s=zeros(nz,nx);
% ttt_r=zeros(nz,nx,ng);
% ttt_s=time_table(sx,sz,nx,nz,dx)/v0;

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
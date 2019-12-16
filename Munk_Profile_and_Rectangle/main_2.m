clear all
close all
clc

zhat = 2*([1:50:5000]-1300)./1300;
vel = 1500.*(1+0.00737.*(zhat - 1 + exp(-zhat)));

vel=repmat(vel',[1 201]);
vel(20:40,50:100)= 1500;

[nz,nx]=size(vel);

dx=5;
x = (0:nx-1)*dx;
z = (0:nz-1)*dx;

sx = (nx-1)/2*dx; 
sz = 0;

gx=(0:2:(nx-1))*dx;
gz=zeros(size(gx));
nbc=40; 
nt=2001; 
dt=0.0005;
isFS=false;

freq=25; 
s=ricker(freq,dt);

figure(1);
set(gcf,'position',[0 0 800 400]);
subplot(221);
imagesc(x,z,vel);
colorbar;
xlabel('X (m)'); 
ylabel('Z (m)'); title('velocity');
figure(1);subplot(222);plot((0:numel(s)-1)*dt,s);
xlabel('Time (s)'); ylabel('Amplitude');title('wavelet');
tic;
seis=a2d_mod_abc28_snapshot(vel,nbc,dx,nt,dt,s,sx,sz,gx,gz,isFS);
toc;
% figure
% contourf(vel)
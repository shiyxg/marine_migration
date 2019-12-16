% for synthetic test  Dec.3rd
vel=[repmat(1000,[1,30]), repmat(1200,[1,30]), repmat(1500,[1,21])];
vel=repmat(vel',[1 201]);

for i= 71:85
        vel(31:30+(i-70),i)=1000;
end
for i= 86:114
        vel(31:45,i)=1000;
end
for i=115:129
        vel(31:30+(130-i),i)=1000;    
end
figure;imagesc(vel);xlabel('X (m)','FontSize',14); ylabel('Z (m)','FontSize',14);
h=colorbar; set(get(h,'title'),'string','velocity(m/s)','FontSize',14);%colorbar('velocity(m/s)');
title('Sythetic Velocity Model','FontSize',18)
[vel_ss,~]=vel_smooth(vel,11,11,5); 
figure;imagesc(vel_ss);xlabel('X (m)','FontSize',14); ylabel('Z (m)','FontSize',14);
h=colorbar; set(get(h,'title'),'string','velocity(m/s)','FontSize',14);
title('Smooth Velocity Model','FontSize',20)%colorbar;
s=1./vel_ss;


v0=1000;
nx=201;nz=81;
dx=5; x = (0:nx-1)*dx; z = (0:nz-1)*dx; 
nt=2001; dt=0.0005; t=(0:nt-1)*dt;


% gx=(0:2:(nx-1))*dx; gz=zeros(size(gx)); ng=numel(gx);
% sx=(0:20:(nx-1))*dx; 
% sx=(0:100:nx-1)*dx;sz=zeros(size(sx));ns=numel(sx);
% gx=(0:2:nx-1)*dx;gz=zeros(size(gx)); ng=numel(gx);

gx=(0:2:(nx-1))*dx; gz=zeros(size(gx)); ng=numel(gx);
sx=(0:20:(nx-1))*dx; ns=numel(sx); sz=zeros(size(sx));

freq=25; 
source=ricker(freq,dt);
source=expand_source(source,nt);
 
ttt_s=zeros(nz,nx,ns);
ttt_r=zeros(nz,nx,ng);

for i=1:ns
    ttt_s(:,:,i)=eikonal2d(s,sx(i),sz(i),dx);
end

for i=1:ng
    ttt_r(:,:,i)=eikonal2d(s,gx(i),gz(i),dx);
end

load seis_refl;
d=seis_refl;

iter=10;
mig0=zeros(nz,nx);
[mig_lsm,res]=lsm_over(d,mig0,source,sx,sz,gx,gz,ttt_s,ttt_r,dt,dx,iter);
figure(6);
imagesc(mig_lsm);colormap(gray);title('LSM-Kirchhoff Image','FontSize',18);
xlabel('X (m)','FontSize',14); ylabel('Z (m)','FontSize',14); 
figure(7);
plot(1:1:iter,res);title('LSM-Kirchhoff Residual Curve','FontSize',18);
xlabel('Iteration','FontSize',14); ylabel('residual','FontSize',14);




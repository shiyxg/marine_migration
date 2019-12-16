function ttt=eikonal3d(sss,sx,sy,sz,dx)

%EIKONAL: 3D Eikonal solver of direct arrival.
%
%  t_out=eikonal3d(sss,sx,sy,sz,dx);
%
%  IN   sss   : slowness model
%       sx    : source X
%       sy    : source Y
%       sz    : source Z
%       dx    : x integral for calculation
%
%  OUT  ttt(:,:,:) : traveltime table of this source 
%
%  Example
%    nx = 7; ny=7; nz = 7; dx = 5; dy = 5; dz = 5; 
%    vel = 1000.0 * ones(nz,ny,nx); %vel((round(nz+1)/2):end,:)=2500.0;
%    sss=1.0/vel;
%    sx=(nx-1)/2*dx;sy=(ny-1)/2*dy;sz=(nz-1)/2*dz;
%    ttt=eikonal3d(sss,sx,sy,sz,dx);
%    figure;subplot(121);surface(ttt);
%    subplot(122);surface(ttt);
%
[nz,ny,nx]=size(sss);

h=dx; 

ttt=zeros(nz,ny,nx); % claim the matrix of t for cacultation

% decide the closest node point to the source position
isx=round(sx/h)+1;isy=round(sy/h)+1;isz=round(sz/h)+1; 

% calculate the traveltime around the source position up to 5X5
ttt=eik3d_around_source(ttt,sss,sx,sy,sz,isx,isy,isz,h);
% expand the wave front
ttt=eik3d_square_expand(ttt,sss,isx,isy,isz,h);

end

function ttt=eik3d_around_source(ttt,sss,sx,sy,sz,isx,isy,isz,h)

[nz,ny,nx]=size(sss);

ix1=max(isx-2,1);ix2=min(isx+2,nx);
iy1=max(isy-2,1);iy2=min(isy+2,ny);
iz1=max(isz-2,1);iz2=min(isz+2,nz);

for i3=ix1:ix2
    iix1=min(i3,isx);iix2=max(i3,isx);
    for i2=iy1:iy2
        iiy1=min(i2,isy);iiy2=max(i2,isy);
        for i1=iz1:iz2
            iiz1=min(i1,isz);iiz2=max(i1,isz);
            ns=(iix2-iix1+1)*(iiy2-iiy1+1)*(iiz2-iiz1+1);
            s_ave=sum(sum(sum(sss(iiz1:iiz2,iiy1:iiy2,iix1:iix2))))/ns;
            ttt(i1,i2,i3)=( ((i3-1)*h-sx)^2 ...
                        + ((i2-1)*h-sy)^2 ...
                        + ((i1-1)*h-sz)^2 )^0.5*s_ave;
        end
    end
end

end

% function ttt=eik3d_around_source(ttt,sss,sx,sy,sz,isx,isy,isz,h)
% 
% [nz,ny,nx]=size(sss);
% 
% ix1=max(isx-3,1);ix2=min(isx+3,nx);
% iy1=max(isy-3,1);iy2=min(isy+3,ny);
% iz1=max(isz-3,1);iz2=min(isz+3,nz);
% 
% for i3=ix1:ix2
%     iix1=min(i3,isx);iix2=max(i3,isx);
%     for i2=iy1:iy2
%         iiy1=min(i2,isy);iiy2=max(i2,isy);
%         for i1=iz1:iz2
%             iiz1=min(i1,isz);iiz2=max(i1,isz);
%             ns=(iix2-iix1+1)*(iiy2-iiy1+1)*(iiz2-iiz1+1);
%             s_ave=sum(sum(sum(sss(iiz1:iiz2,iiy1:iiy2,iix1:iix2))))/ns;
%             ttt(i1,i2,i3)=( ((i3-1)*h-sx)^2 ...
%                         + ((i2-1)*h-sy)^2 ...
%                         + ((i1-1)*h-sz)^2 )^0.5*s_ave;
%         end
%     end
% end
% 
% end



function ttt=eik3d_square_expand(ttt,sss,isx,isy,isz,h)

[nz,ny,nx]=size(sss);

for i=3:max([isx-1,isy-1,isz-1,nx-isx,ny-isy,nz-isz])
    % top expand
    if (isz-i)>=1
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        ss0=sss(isz-i,iy1:iy2,ix1:ix2);
        ss1=sss(isz-(i-1),iy1:iy2,ix1:ix2);
%         ss2=sss(isz-(i-2),iy1:iy2,ix1:ix2);
        tt1=ttt(isz-(i-1),iy1:iy2,ix1:ix2);
%         tt2=ttt(isz-(i-2),iy1:iy2,ix1:ix2);
%         tt0=eik3d_surface_expand(tt1,tt2,ss0,ss1,ss2,h);   
        tt0=eik3d_surface_expand(tt1,ss0,ss1,h);        
        ttt(isz-i,iy1:iy2,ix1:ix2)=tt0;
    end
    % bottom expand
    if (isz+i)<=nz
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        ss0=sss(isz+i,iy1:iy2,ix1:ix2);
        ss1=sss(isz+(i-1),iy1:iy2,ix1:ix2);
%         ss2=sss(isz+(i-2),iy1:iy2,ix1:ix2);
        tt1=ttt(isz+(i-1),iy1:iy2,ix1:ix2);
%         tt2=ttt(isz+(i-2),iy1:iy2,ix1:ix2);
%         tt0=eik3d_surface_expand(tt1,tt2,ss0,ss1,ss2,h);
        tt0=eik3d_surface_expand(tt1,ss0,ss1,h);
        ttt(isz+i,iy1:iy2,ix1:ix2)=tt0;
    end    
    % left expand
    if (isx-i)>=1
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        ss0=sss(iz1:iz2,iy1:iy2,isx-i);
        ss1=sss(iz1:iz2,iy1:iy2,isx-(i-1));
%         ss2=sss(iz1:iz2,iy1:iy2,isx-(i-2));
        tt1=ttt(iz1:iz2,iy1:iy2,isx-(i-1));
%         tt2=ttt(iz1:iz2,iy1:iy2,isx-(i-2));
%         tt0=eik3d_surface_expand(tt1,tt2,ss0,ss1,ss2,h);
        tt0=eik3d_surface_expand(tt1,ss0,ss1,h);
        ttt(iz1:iz2,iy1:iy2,isx-i)=tt0;
    end     
    % right expand
    if (isx+i)<=nx
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        ss0=sss(iz1:iz2,iy1:iy2,isx+i);
        ss1=sss(iz1:iz2,iy1:iy2,isx+(i-1));
%         ss2=sss(iz1:iz2,iy1:iy2,isx+(i-2));
        tt1=ttt(iz1:iz2,iy1:iy2,isx+(i-1));
%         tt2=ttt(iz1:iz2,iy1:iy2,isx+(i-2));
%         tt0=eik3d_surface_expand(tt1,tt2,ss0,ss1,ss2,h);
        tt0=eik3d_surface_expand(tt1,ss0,ss1,h);
        ttt(iz1:iz2,iy1:iy2,isx+i)=tt0;
    end       
    % back expand
    if (isy-i)>=1
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        ss0=sss(iz1:iz2,isy-i,ix1:ix2);
        ss1=sss(iz1:iz2,isy-(i-1),ix1:ix2);
%         ss2=sss(iz1:iz2,isy-(i-2),ix1:ix2);
        tt1=ttt(iz1:iz2,isy-(i-1),ix1:ix2);
%         tt2=ttt(iz1:iz2,isy-(i-2),ix1:ix2);
%         tt0=eik3d_surface_expand(tt1,tt2,ss0,ss1,ss2,h);
        tt0=eik3d_surface_expand(tt1,ss0,ss1,h);
        ttt(iz1:iz2,isy-i,ix1:ix2)=tt0;
    end     
    % front expand
    if (isy+i)<=ny
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        ss0=sss(iz1:iz2,isy+i,ix1:ix2);
        ss1=sss(iz1:iz2,isy+(i-1),ix1:ix2);
%         ss2=sss(iz1:iz2,isy+(i-2),ix1:ix2);
        tt1=ttt(iz1:iz2,isy+(i-1),ix1:ix2);
%         tt2=ttt(iz1:iz2,isy+(i-2),ix1:ix2);
%         tt0=eik3d_surface_expand(tt1,tt2,ss0,ss1,ss2,h);
        tt0=eik3d_surface_expand(tt1,ss0,ss1,h);
        ttt(iz1:iz2,isy+i,ix1:ix2)=tt0;
    end       
    
    % eadge expand describe by the share face
    % top back face
    if (isz-i) >= 1 && (isy-i) >=1
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        s0=sss(isz-i,isy-i,ix1:ix2);
        s1=sss(isz-(i-1),isy-i,ix1:ix2);
%         s11=sss(isz-(i-2),isy-i,ix1:ix2);
        s2=sss(isz-i,isy-(i-1),ix1:ix2);
%         s22=sss(isz-i,isy-(i-2),ix1:ix2);
        s3=sss(isz-(i-1),isy-(i-1),ix1:ix2);
        t1=ttt(isz-(i-1),isy-i,ix1:ix2);
%         t11=ttt(isz-(i-2),isy-i,ix1:ix2);
        t2=ttt(isz-i,isy-(i-1),ix1:ix2);
%         t22=ttt(isz-i,isy-(i-2),ix1:ix2);
        t3=ttt(isz-(i-1),isy-(i-1),ix1:ix2);
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz-i,isy-i,ix1:ix2)=t0;
    end
     % top front face
    if (isz-i) >= 1 && (isy+i) <=ny
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        s0=sss(isz-i,isy+i,ix1:ix2);
        s1=sss(isz-(i-1),isy+i,ix1:ix2);
%         s11=sss(isz-(i-2),isy+i,ix1:ix2);
        s2=sss(isz-i,isy+(i-1),ix1:ix2);
%         s22=sss(isz-i,isy+(i-2),ix1:ix2);
        s3=sss(isz-(i-1),isy+(i-1),ix1:ix2);
        t1=ttt(isz-(i-1),isy+i,ix1:ix2);
%         t11=ttt(isz-(i-2),isy+i,ix1:ix2);
        t2=ttt(isz-i,isy+(i-1),ix1:ix2);
%         t22=ttt(isz-i,isy+(i-2),ix1:ix2);
        t3=ttt(isz-(i-1),isy+(i-1),ix1:ix2);
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz-i,isy+i,ix1:ix2)=t0;
    end   
    % top left face
    if (isz-i) >= 1 && (isx-i) >=1
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        s0=sss(isz-i,iy1:iy2,isx-i);
        s1=sss(isz-(i-1),iy1:iy2,isx-i);
%         s11=sss(isz-(i-2),iy1:iy2,isx-i);
        s2=sss(isz-i,iy1:iy2,isx-(i-1));
%         s22=sss(isz-i,iy1:iy2,isx-(i-2));
        s3=sss(isz-(i-1),iy1:iy2,isx-(i-1));
        t1=ttt(isz-(i-1),iy1:iy2,isx-i);
%         t11=ttt(isz-(i-2),iy1:iy2,isx-i);
        t2=ttt(isz-i,iy1:iy2,isx-(i-1));
%         t22=ttt(isz-i,iy1:iy2,isx-(i-2));
        t3=ttt(isz-(i-1),iy1:iy2,isx-(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz-i,iy1:iy2,isx-i)=t0;
    end  
    % top right face
    if (isz-i) >= 1 && (isx+i) <= ny
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        s0=sss(isz-i,iy1:iy2,isx+i);
        s1=sss(isz-(i-1),iy1:iy2,isx+i);
%         s11=sss(isz-(i-2),iy1:iy2,isx+i);
        s2=sss(isz-i,iy1:iy2,isx+(i-1));
%         s22=sss(isz-i,iy1:iy2,isx+(i-2));
        s3=sss(isz-(i-1),iy1:iy2,isx+(i-1));
        t1=ttt(isz-(i-1),iy1:iy2,isx+i);
%         t11=ttt(isz-(i-2),iy1:iy2,isx+i);
        t2=ttt(isz-i,iy1:iy2,isx+(i-1));
%         t22=ttt(isz-i,iy1:iy2,isx+(i-2));
        t3=ttt(isz-(i-1),iy1:iy2,isx+(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz-i,iy1:iy2,isx+i)=t0;
    end  
    % front left face
    if (isy-i) >= 1 && (isx-i) >=1
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        s0=sss(iz1:iz2,isy-i,isx-i);
        s1=sss(iz1:iz2,isy-(i-1),isx-i);
%         s11=sss(iz1:iz2,isy-(i-2),isx-i);
        s2=sss(iz1:iz2,isy-i,isx-(i-1));
%         s22=sss(iz1:iz2,isy-i,isx-(i-2));
        s3=sss(iz1:iz2,isy-(i-1),isx-(i-1));
        t1=ttt(iz1:iz2,isy-(i-1),isx-i);
%         t11=ttt(iz1:iz2,isy-(i-2),isx-i);
        t2=ttt(iz1:iz2,isy-i,isx-(i-1));
%         t22=ttt(iz1:iz2,isy-i,isx-(i-2));
        t3=ttt(iz1:iz2,isy-(i-1),isx-(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(iz1:iz2,isy-i,isx-i)=t0;
    end      
    % front right face
    if (isy-i) >= 1 && (isx+i) <= nx
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        s0=sss(iz1:iz2,isy-i,isx+i);
        s1=sss(iz1:iz2,isy-(i-1),isx+i);
%         s11=sss(iz1:iz2,isy-(i-2),isx+i);
        s2=sss(iz1:iz2,isy-i,isx+(i-1));
%         s22=sss(iz1:iz2,isy-i,isx+(i-2));
        s3=sss(iz1:iz2,isy-(i-1),isx+(i-1));
        t1=ttt(iz1:iz2,isy-(i-1),isx+i);
%         t11=ttt(iz1:iz2,isy-(i-2),isx+i);
        t2=ttt(iz1:iz2,isy-i,isx+(i-1));
%         t22=ttt(iz1:iz2,isy-i,isx+(i-2));
        t3=ttt(iz1:iz2,isy-(i-1),isx+(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(iz1:iz2,isy-i,isx+i)=t0;
    end      
    % back left face
    if (isy+i) <= ny && (isx-i) >=1
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        s0=sss(iz1:iz2,isy+i,isx-i);
        s1=sss(iz1:iz2,isy+(i-1),isx-i);
%         s11=sss(iz1:iz2,isy+(i-2),isx-i);
        s2=sss(iz1:iz2,isy+i,isx-(i-1));
%         s22=sss(iz1:iz2,isy+i,isx-(i-2));
        s3=sss(iz1:iz2,isy+(i-1),isx-(i-1));
        t1=ttt(iz1:iz2,isy+(i-1),isx-i);
%         t11=ttt(iz1:iz2,isy+(i-2),isx-i);
        t2=ttt(iz1:iz2,isy+i,isx-(i-1));
%         t22=ttt(iz1:iz2,isy+i,isx-(i-2));
        t3=ttt(iz1:iz2,isy+(i-1),isx-(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(iz1:iz2,isy+i,isx-i)=t0;
    end      
    % back right face
    if (isy+i) <= ny && (isx+i) <= nx
        iz1=max(isz-(i-1),1);iz2=min(isz+(i-1),nz);
        s0=sss(iz1:iz2,isy+i,isx+i);
        s1=sss(iz1:iz2,isy+(i-1),isx+i);
%         s11=sss(iz1:iz2,isy+(i-2),isx+i);
        s2=sss(iz1:iz2,isy+i,isx+(i-1));
%         s22=sss(iz1:iz2,isy+i,isx+(i-2));
        s3=sss(iz1:iz2,isy+(i-1),isx+(i-1));
        t1=ttt(iz1:iz2,isy+(i-1),isx+i);
%         t11=ttt(iz1:iz2,isy+(i-2),isx+i);
        t2=ttt(iz1:iz2,isy+i,isx+(i-1));
%         t22=ttt(iz1:iz2,isy+i,isx+(i-2));
        t3=ttt(iz1:iz2,isy+(i-1),isx+(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(iz1:iz2,isy+i,isx+i)=t0;
    end      
    % bottom back face
    if (isz+i) <= nz && (isy-i) >=1
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        s0=sss(isz+i,isy-i,ix1:ix2);
        s1=sss(isz+(i-1),isy-i,ix1:ix2);
%         s11=sss(isz+(i-2),isy-i,ix1:ix2);
        s2=sss(isz+i,isy-(i-1),ix1:ix2);
%         s22=sss(isz+i,isy-(i-2),ix1:ix2);
        s3=sss(isz+(i-1),isy-(i-1),ix1:ix2);
        t1=ttt(isz+(i-1),isy-i,ix1:ix2);
%         t11=ttt(isz+(i-2),isy-i,ix1:ix2);
        t2=ttt(isz+i,isy-(i-1),ix1:ix2);
%         t22=ttt(isz+i,isy-(i-2),ix1:ix2);
        t3=ttt(isz+(i-1),isy-(i-1),ix1:ix2);
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz+i,isy-i,ix1:ix2)=t0;
    end
    % bottom front face
    if (isz+i) <= nz && (isy+i) <=ny
        ix1=max(isx-(i-1),1);ix2=min(isx+(i-1),nx);
        s0=sss(isz+i,isy+i,ix1:ix2);
        s1=sss(isz+(i-1),isy+i,ix1:ix2);
%         s11=sss(isz+(i-2),isy+i,ix1:ix2);
        s2=sss(isz+i,isy+(i-1),ix1:ix2);
%         s22=sss(isz+i,isy+(i-2),ix1:ix2);
        s3=sss(isz+(i-1),isy+(i-1),ix1:ix2);
        t1=ttt(isz+(i-1),isy+i,ix1:ix2);
%         t11=ttt(isz+(i-2),isy+i,ix1:ix2);
        t2=ttt(isz+i,isy+(i-1),ix1:ix2);
%         t22=ttt(isz+i,isy+(i-2),ix1:ix2);
        t3=ttt(isz+(i-1),isy+(i-1),ix1:ix2);
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz+i,isy+i,ix1:ix2)=t0;
    end   
    % bottom left face
    if (isz+i) <= nz && (isx-i) >=1
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        s0=sss(isz+i,iy1:iy2,isx-i);
        s1=sss(isz+(i-1),iy1:iy2,isx-i);
%         s11=sss(isz+(i-2),iy1:iy2,isx-i);
        s2=sss(isz+i,iy1:iy2,isx-(i-1));
%         s22=sss(isz+i,iy1:iy2,isx-(i-2));
        s3=sss(isz+(i-1),iy1:iy2,isx-(i-1));
        t1=ttt(isz+(i-1),iy1:iy2,isx-i);
%         t11=ttt(isz+(i-2),iy1:iy2,isx-i);
        t2=ttt(isz+i,iy1:iy2,isx-(i-1));
%         t22=ttt(isz+i,iy1:iy2,isx-(i-2));
        t3=ttt(isz+(i-1),iy1:iy2,isx-(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz+i,iy1:iy2,isx-i)=t0;
    end  
    % bottom right face
    if (isz+i) <= nz && (isx+i) <= ny
        iy1=max(isy-(i-1),1);iy2=min(isy+(i-1),ny);
        s0=sss(isz+i,iy1:iy2,isx+i);
        s1=sss(isz+(i-1),iy1:iy2,isx+i);
%         s11=sss(isz+(i-2),iy1:iy2,isx+i);
        s2=sss(isz+i,iy1:iy2,isx+(i-1));
%         s22=sss(isz+i,iy1:iy2,isx+(i-2));
        s3=sss(isz+(i-1),iy1:iy2,isx+(i-1));
        t1=ttt(isz+(i-1),iy1:iy2,isx+i);
%         t11=ttt(isz+(i-2),iy1:iy2,isx+i);
        t2=ttt(isz+i,iy1:iy2,isx+(i-1));
%         t22=ttt(isz+i,iy1:iy2,isx+(i-2));
        t3=ttt(isz+(i-1),iy1:iy2,isx+(i-1));
%         t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h);
        t0=eik3d_eadge_expand(s0,s1,s2,s3,t1,t2,t3,h);
        ttt(isz+i,iy1:iy2,isx+i)=t0;
    end  
    
    % conor expand, use the face around it to represent
    % top back left
    if ((isx-i) >= 1 && (isy-i) >= 1 && (isz-i) >= 1)
        ttt(isz-i,isy-i,isx-i)=eik3d_coner_expand(sss(isz-i,isy-i,isx-i),...
                                            sss(isz-i+1,isy-i,isx-i),sss(isz-i+2,isy-i,isx-i),...
                                            sss(isz-i,isy-i+1,isx-i),sss(isz-i,isy-i+2,isx-i),...
                                            sss(isz-i,isy-i,isx-i+1),sss(isz-i,isy-i,isx-i+2),...
                                            sss(isz-i+1,isy-i,isx-i+1),sss(isz-i,isy-i+1,isx-i+1),...
                                            sss(isz-i+1,isy-i+1,isx-i),...
                                            ttt(isz-i+1,isy-i,isx-i),ttt(isz-i+2,isy-i,isx-i),...
                                            ttt(isz-i,isy-i+1,isx-i),ttt(isz-i,isy-i+2,isx-i),...
                                            ttt(isz-i,isy-i,isx-i+1),ttt(isz-i,isy-i,isx-i+2),...
                                            ttt(isz-i+1,isy-i,isx-i+1),ttt(isz-i,isy-i+1,isx-i+1),...
                                            ttt(isz-i+1,isy-i+1,isx-i),h);
    end
    % top back right
    if ((isx+i) <= nx && (isy-i) >= 1 && (isz-i) >= 1) 
        ttt(isz-i,isy-i,isx+i)=eik3d_coner_expand(sss(isz-i,isy-i,isx+i),...
                                            sss(isz-i+1,isy-i,isx+i),sss(isz-i+2,isy-i,isx+i),...
                                            sss(isz-i,isy-i+1,isx+i),sss(isz-i,isy-i+2,isx+i),...
                                            sss(isz-i,isy-i,isx+i-1),sss(isz-i,isy-i,isx+i-2),...
                                            sss(isz-i+1,isy-i,isx+i-1),sss(isz-i,isy-i+1,isx+i-1),...
                                            sss(isz-i+1,isy-i+1,isx+i),...
                                            ttt(isz-i+1,isy-i,isx+i),ttt(isz-i+2,isy-i,isx+i),...
                                            ttt(isz-i,isy-i+1,isx+i),ttt(isz-i,isy-i+2,isx+i),...
                                            ttt(isz-i,isy-i,isx+i-1),ttt(isz-i,isy-i,isx+i-2),...
                                            ttt(isz-i+1,isy-i,isx+i-1),ttt(isz-i,isy-i+1,isx+i-1),...
                                            ttt(isz-i+1,isy-i+1,isx+i),h);        
                                        
    end
    % top front left
    if ((isx-i) >= 1 && (isy+i) <= ny && (isz-i) >= 1)
        ttt(isz-i,isy+i,isx-i)=eik3d_coner_expand(sss(isz-i,isy+i,isx-i),...
                                            sss(isz-i+1,isy+i,isx-i),sss(isz-i+2,isy+i,isx-i),...
                                            sss(isz-i,isy+i-1,isx-i),sss(isz-i,isy+i-2,isx-i),...
                                            sss(isz-i,isy+i,isx-i+1),sss(isz-i,isy+i,isx-i+2),...
                                            sss(isz-i+1,isy+i,isx-i+1),sss(isz-i,isy+i-1,isx-i+1),...
                                            sss(isz-i+1,isy+i-1,isx-i),...
                                            ttt(isz-i+1,isy+i,isx-i),ttt(isz-i+2,isy+i,isx-i),...
                                            ttt(isz-i,isy+i-1,isx-i),ttt(isz-i,isy+i-2,isx-i),...
                                            ttt(isz-i,isy+i,isx-i+1),ttt(isz-i,isy+i,isx-i+2),...
                                            ttt(isz-i+1,isy+i,isx-i+1),ttt(isz-i,isy+i-1,isx-i+1),...
                                            ttt(isz-i+1,isy+i-1,isx-i),h);
    end
    % top front right
    if ((isx+i) >= 1 && (isy+i) <= ny && (isz-i) >= 1)
        ttt(isz-i,isy+i,isx+i)=eik3d_coner_expand(sss(isz-i,isy+i,isx+i),...
                                            sss(isz-i+1,isy+i,isx+i),sss(isz-i+2,isy+i,isx+i),...
                                            sss(isz-i,isy+i-1,isx+i),sss(isz-i,isy+i-2,isx+i),...
                                            sss(isz-i,isy+i,isx-i+1),sss(isz-i,isy+i,isx-i+2),...
                                            sss(isz-i+1,isy+i,isx+i-1),sss(isz-i,isy+i-1,isx+i-1),...
                                            sss(isz-i+1,isy+i-1,isx+i),...
                                            ttt(isz-i+1,isy+i,isx+i),ttt(isz-i+2,isy+i,isx+i),...
                                            ttt(isz-i,isy+i-1,isx+i),ttt(isz-i,isy+i-2,isx+i),...
                                            ttt(isz-i,isy+i,isx-i+1),ttt(isz-i,isy+i,isx-i+2),...
                                            ttt(isz-i+1,isy+i,isx+i-1),ttt(isz-i,isy+i-1,isx+i-1),...
                                            ttt(isz-i+1,isy+i-1,isx+i),h);
    end
    % bottom back left
    if ((isx-i) >= 1 && (isy-i) >= 1 && (isz+i) <= nz)
        ttt(isz+i,isy-i,isx-i)=eik3d_coner_expand(sss(isz+i,isy-i,isx-i),...
                                            sss(isz+i-1,isy-i,isx-i),sss(isz+i-2,isy-i,isx-i),...
                                            sss(isz+i,isy-i+1,isx-i),sss(isz+i,isy-i+2,isx-i),...
                                            sss(isz+i,isy-i,isx-i+1),sss(isz+i,isy-i,isx-i+2),...
                                            sss(isz+i-1,isy-i,isx-i+1),sss(isz+i,isy-i+1,isx-i+1),...
                                            sss(isz+i-1,isy-i+1,isx-i),...
                                            ttt(isz+i-1,isy-i,isx-i),ttt(isz+i-2,isy-i,isx-i),...
                                            ttt(isz+i,isy-i+1,isx-i),ttt(isz+i,isy-i+2,isx-i),...
                                            ttt(isz+i,isy-i,isx-i+1),ttt(isz+i,isy-i,isx-i+2),...
                                            ttt(isz+i-1,isy-i,isx-i+1),ttt(isz+i,isy-i+1,isx-i+1),...
                                            ttt(isz+i-1,isy-i+1,isx-i),h);
    end
    % bottom back right
    if ((isx+i) <= nx && (isy-i) >= 1 && (isz+i) <= nz) 
        ttt(isz+i,isy-i,isx+i)=eik3d_coner_expand(sss(isz+i,isy-i,isx+i),...
                                            sss(isz+i-1,isy-i,isx+i),sss(isz+i-2,isy-i,isx+i),...
                                            sss(isz+i,isy-i+1,isx+i),sss(isz+i,isy-i+2,isx+i),...
                                            sss(isz+i,isy-i,isx+i-1),sss(isz+i,isy-i,isx+i-2),...
                                            sss(isz+i-1,isy-i,isx+i-1),sss(isz+i,isy-i+1,isx+i-1),...
                                            sss(isz+i-1,isy-i+1,isx+i),...
                                            ttt(isz+i-1,isy-i,isx+i),ttt(isz+i-2,isy-i,isx+i),...
                                            ttt(isz+i,isy-i+1,isx+i),ttt(isz+i,isy-i+2,isx+i),...
                                            ttt(isz+i,isy-i,isx+i-1),ttt(isz+i,isy-i,isx+i-2),...
                                            ttt(isz+i-1,isy-i,isx+i-1),ttt(isz+i,isy-i+1,isx+i-1),...
                                            ttt(isz+i-1,isy-i+1,isx+i),h);        
    end
    % bottom front left
    if ((isx-i) >= 1 && (isy+i) <= ny && (isz+i) <= nz)
        ttt(isz+i,isy+i,isx-i)=eik3d_coner_expand(sss(isz+i,isy+i,isx-i),...
                                            sss(isz+i-1,isy+i,isx-i),sss(isz+i-2,isy+i,isx-i),...
                                            sss(isz+i,isy+i-1,isx-i),sss(isz+i,isy+i-2,isx-i),...
                                            sss(isz+i,isy+i,isx-i+1),sss(isz+i,isy+i,isx-i+2),...
                                            sss(isz+i-1,isy+i,isx-i+1),sss(isz+i,isy+i-1,isx-i+1),...
                                            sss(isz+i-1,isy+i-1,isx-i),...
                                            ttt(isz+i-1,isy+i,isx-i),ttt(isz+i-2,isy+i,isx-i),...
                                            ttt(isz+i,isy+i-1,isx-i),ttt(isz+i,isy+i-2,isx-i),...
                                            ttt(isz+i,isy+i,isx-i+1),ttt(isz+i,isy+i,isx-i+2),...
                                            ttt(isz+i-1,isy+i,isx-i+1),ttt(isz+i,isy+i-1,isx-i+1),...
                                            ttt(isz+i-1,isy+i-1,isx-i),h);
    end
    % bottom front right
    if ((isx+i) >= 1 && (isy+i) <= ny && (isz+i) <= nz)
        ttt(isz+i,isy+i,isx+i)=eik3d_coner_expand(sss(isz+i,isy+i,isx+i),...
                                            sss(isz+i-1,isy+i,isx+i),sss(isz+i-2,isy+i,isx+i),...
                                            sss(isz+i,isy+i-1,isx+i),sss(isz+i,isy+i-2,isx+i),...
                                            sss(isz+i,isy+i,isx-i+1),sss(isz+i,isy+i,isx-i+2),...
                                            sss(isz+i-1,isy+i,isx+i-1),sss(isz+i,isy+i-1,isx+i-1),...
                                            sss(isz+i-1,isy+i-1,isx+i),...
                                            ttt(isz+i-1,isy+i,isx+i),ttt(isz+i-2,isy+i,isx+i),...
                                            ttt(isz+i,isy+i-1,isx+i),ttt(isz+i,isy+i-2,isx+i),...
                                            ttt(isz+i,isy+i,isx-i+1),ttt(isz+i,isy+i,isx-i+2),...
                                            ttt(isz+i-1,isy+i,isx+i-1),ttt(isz+i,isy+i-1,isx+i-1),...
                                            ttt(isz+i-1,isy+i-1,isx+i),h);
    end    
end

end

function t0=eik3d_surface_expand(t1,s0,s1,h)
% function t0=eik3d_surface_expand(t1,t2,s0,s1,s2,h)

t1=squeeze(t1);%t2=squeeze(t2);
s0=squeeze(s0);s1=squeeze(s1);%s2=squeeze(s2);

t0=zeros(size(t1));
[n1,n2]=size(t1);
h2=h*h;
hsq=2^0.5 * h;

% find the minimal points of t1
[min_loc_fast,min_num_fast,min_close_fast]=eik3d_local_min(t1,is1,is2);
[max_loc_fast,max_num_fast,max_close_fast]=eik3d_local_max(t1,is1,is2);

% calculate the traveltime just outside the local minimum points
for i2=1:n2
    for imin=1:min_num_fast(i2)
        t0(min_loc_fast(imin,1),min_loc(imin,2))=eik3d_expand_outside(min_loc(imin,1:2),t1,s0,s1,h);
    end
end

% find the minimal point the most closest to the source point
i_close=eik3d_closest(min_loc,min_num,is1,is2);

% use this point to decide 


% use the i_close to divide the 2D surface into 2 groups  

% sweep left and right



for i1=1:n1
    for i2=1:n2
        if (i1==1)
            if (i2==1)
%                 s=(s0(i1,i2)+s1(i1,i2)+s1(i1+1,i2)+s1(i1,i2+1)+s2(i1,i2))/5.0;
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1+1,i2)+s1(i1,i2+1))/4.0;
                delta=4*h2*s^2-4*(t1(i1,i2)-t1(i1+1,i2))^2 - 4*(t1(i1,i2)-t1(i1,i2+1))^2;
                if (delta<0.0)
                    t0(i1,i2)=t1(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1+1,i2)+hsq*s,t1(i1,i2+1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end
            elseif (i2>1 && i2<n2)
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1+1,i2)+s1(i1,i2+1)+s1(i1,i2-1)+s2(i1,i2))/6.0;
                delta=4*h2*s^2-4*(t1(i1,i2)-t1(i1+1,i2))^2- (t1(i1,i2-1)-t1(i1,i2+1))^2; 
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1+1,i2)+hsq*s,t1(i1,i2+1)+hsq*s, t1(i1,i2-1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end                
            else
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1+1,i2)+s1(i1,i2-1)+s2(i1,i2))/5.0;
                delta=4*h2*s^2-4*(t1(i1,i2)-t1(i1+1,i2))^2 - 4*(t1(i1,i2-1)-t1(i1,i2))^2;
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1+1,i2)+hsq*s,t1(i1,i2-1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end     
            end
        elseif (i1>1 && i1<n1)
             if (i2==1)
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1-1,i2)+s1(i1+1,i2)+s1(i1,i2+1)+s2(i1,i2))/6.0;
                delta=4*h2*s^2-(t1(i1-1,i2)-t1(i1+1,i2))^2 - 4*(t1(i1,i2)-t1(i1,i2+1))^2;
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1-1,i2)+hsq*s, t1(i1+1,i2)+hsq*s,t1(i1,i2+1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end
            elseif (i2>1 && i2<n2)
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1-1,i2)+s1(i1+1,i2)+s1(i1,i2+1)+s1(i1,i2-1)+s2(i1,i2))/7.0;
                delta=4*h2*s^2-(t1(i1-1,i2)-t1(i1+1,i2))^2- (t1(i1,i2-1)-t1(i1,i2+1))^2; 
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1-1,i2)+hsq*s,t1(i1+1,i2)+hsq*s,t1(i1,i2+1)+hsq*s, t1(i1,i2-1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end                
            else
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1-1,i2)+s1(i1+1,i2)+s1(i1,i2-1)+s2(i1,i2))/6.0;
                delta=4*h2*s^2-(t1(i1-1,i2)-t1(i1+1,i2))^2 - 4*(t1(i1,i2-1)-t1(i1,i2))^2;
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1-1,i2)+hsq*s,t1(i1+1,i2)+hsq*s,t1(i1,i2-1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end    
            end          
        else
            if (i2==1)
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1-1,i2)+s1(i1,i2+1)+s2(i1,i2))/5.0;
                delta=4*h2*s^2-4*(t1(i1-1,i2)-t1(i1,i2))^2 - 4*(t1(i1,i2)-t1(i1,i2+1))^2;
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1-1,i2)+hsq*s, t1(i1,i2+1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end
            elseif (i2>1 && i2<n2)
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1-1,i2)+s1(i1,i2+1)+s1(i1,i2-1)+s2(i1,i2))/6.0;
                delta=4*h2*s^2-4*(t1(i1-1,i2)-t1(i1,i2))^2- (t1(i1,i2-1)-t1(i1,i2+1))^2; 
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1-1,i2)+hsq*s,t1(i1,i2+1)+hsq*s, t1(i1,i2-1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end                
            else
                s=(s0(i1,i2)+s1(i1,i2)+s1(i1-1,i2)+s1(i1,i2-1)+s2(i1,i2))/5.0;
                delta=4*h2*s^2-4*(t1(i1-1,i2)-t1(i1,i2))^2 - 4*(t1(i1,i2-1)-t1(i1,i2))^2;
                if (delta<0.0)
                    t0(i1,i2)=t2(i1,i2) + min ([t1(i1,i2)+h*s, t1(i1-1,i2)+hsq*s,t1(i1,i2-1)+hsq*s]);
                else
                    t0(i1,i2)=t2(i1,i2) + delta^0.5;
                end     
            end          
        end
    end
end

end
            
        
function t0=eik3d_eadge_expand(s0,s1,s11,s2,s22,s3,t1,t11,t2,t22,t3,h)

n=length(s0);
hsq=2^0.5*h;
h2=h^2;
t0=zeros(size(s0));

t0_temp1=eik3d_eadge_expand_internal(s0,s1,s11,s3,t1,t11,t3,n,h,hsq,h2);
t0_temp2=eik3d_eadge_expand_internal(s0,s2,s22,s3,t2,t22,t3,n,h,hsq,h2);

for i=1:n
    t0(i)=min(t0_temp1(i),t0_temp2(i));
end

end

function t0=eik3d_eadge_expand_internal(s0,s1,s11,s3,t1,t11,t3,n,h,hsq,h2)
    
t0=zeros(size(s0));

for i=1:n
    if i==1
        s=(s0(i)+s1(i)+s1(i+1)+s11(i)+s3(i))/5.0;
        delta=4*h2*s^2 - 4*(t1(i)-t1(i+1))^2 - 4*(t3(i)-t1(i))^2;
        if (delta<0)
            t0(i)=t11(i)+min([t1(i)+h*s,t1(i+1)+hsq*s,t3(i)+hsq*s]);
        else
            t0(i)=t11(i)+delta^0.5;i
        end
    elseif i>1 && i<n
        s=(s0(i)+s1(i-1)+s1(i)+s1(i+1)+s11(i)+s3(i))/6.0;
        delta=4*h2*s^2 - (t1(i-1)-t1(i+1))^2 - 4*(t3(i)-t1(i))^2;
        if (delta<0)
            t0(i)=t11(i)+min([t1(i)+h*s,t1(i-1)+hsq*s,t1(i+1)+hsq*s,t3(i)+hsq*s]);
        else
            t0(i)=t11(i)+delta^0.5;i
        end 
    else
        s=(s0(i)+s1(i-1)+s1(i)+s11(i)+s3(i))/5.0;
        delta=4*h2*s^2 - 4*(t1(i-1)-t1(i))^2 - 4*(t3(i)-t1(i))^2;
        if (delta<0)
            t0(i)=t11(i)+min([t1(i)+h*s,t1(i-1)+hsq*s,t3(i)+hsq*s]);
        else
            t0(i)=t11(i)+delta^0.5;i
        end         
    end
end

end

function t0=eik3d_coner_expand(s0,s1,s11,s2,s22,s3,s33,s13,s23,s12,t1,t11,t2,t22,t3,t33,t13,t23,t12,h)

hsq=2^0.5*h;
h2=h^2;
t0_temp1=eik3d_coner_expand_internal(s0,s1,s11,s12,s13,t1,t11,t12,t13,h,hsq,h2);
t0_temp2=eik3d_coner_expand_internal(s0,s2,s22,s12,s23,t2,t22,t12,t23,h,hsq,h2);
t0_temp3=eik3d_coner_expand_internal(s0,s3,s33,s13,s23,t3,t33,t13,t23,h,hsq,h2);
t0=min([t0_temp1,t0_temp2,t0_temp3]);

end

function t0=eik3d_coner_expand_internal(s0,s1,s11,s12,s13,t1,t11,t12,t13,h,hsq,h2)

s=(s0+s1+s11+s12+s13)/5.0;
delta= 4 *h2 * s^2 - 4*(t1-t12)^2 - 4*(t1-t13)^2;
if (delta <0)
    t0=t11+min([t1+h*s,t12+hsq*s,t13+hsq*s]);
else
    t0=t11+delta^0.5;
end

end




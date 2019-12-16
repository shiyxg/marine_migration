function [dist]=cal_dist(nz,nx,sx,sz,dx)


x=(0:nx-1)*dx; z=(0:nz-1)*dx;
[xx, zz]=meshgrid(x,z);
dist=sqrt((sx-xx).^2+(sz-zz).^2);
dist=max(dist,dx);
end
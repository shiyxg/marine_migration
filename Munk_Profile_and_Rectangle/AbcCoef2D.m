function damp=AbcCoef2D(vel,nbc,dx)
[nzbc,nxbc]=size(vel);
velmin=min(vel(:));
nz=nzbc-2*nbc;nx=nxbc-2*nbc;
a=(nbc-1)*dx;
kappa = 3.0 * velmin * log(10000000.0) / (2.0 * a);
% setup 1D BC damping array
damp1d=kappa*((0:(nbc-1))*dx/a).^2;
% setup 2D BC damping array
damp=zeros(nzbc,nxbc);
% divide the whole area to 9 zones, and 5th is the target zone
%  1   |   2   |   3
%  ------------------
%  4   |   5   |   6
%  ------------------
%  7   |   8   |   9
% fulltill zone 1, 4, 7 and 3, 6, 9
for iz=1:nzbc
    damp(iz,1:nbc)=damp1d(nbc:-1:1);
    damp(iz,nx+nbc+1:nx+2*nbc)=damp1d(:);
end
% full fill zone 2 and 8
for ix=nbc+1:nbc+nx
    damp(1:nbc,ix)=damp1d(nbc:-1:1);
    damp(nbc+nz+1:nz+2*nbc,ix)=damp1d(:);
end
end
function [ time_table ] = time_table( x,z,nx,nz,dx )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

time_table=zeros(nz,nx);
for i=1:nx
    for j=1:nz
         time_table(j,i)=sqrt((x-dx*i)^2+(z-dx*j)^2);  
    end
end


end


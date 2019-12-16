function t=eikonal2d(s,sx,sz,h)
%  Copyright (C) 2010 Center for Subsurface Imaging and Fluid Modeling (CSIM),
%  King Abdullah University of Science and Technology, All rights reserved.
%
%  author:   Xin Wang
%  email:    xin.wang@kaust.edu.sa
%  date:     Feb 15 , 2013
%
%  EIKONAL: Eikonal solver of direct arrival.
%
%  t_out=eikonal2d(vel,sx,sz,dx,dz,dx_out,dz_out);
%
%  IN   s     : slowness model
%       sx    : source X 
%       sz    : source Z
%       dx    : x integral for calculation
%       dz    : z integral for calculation
%       dx_out: x intergral for output
%       dz_out: z intergral for output
%  OUT  t(:,:) : traveltime table of this source 
%
%  Example
%    nx = 201; nz = 101; dx = 5; 
%    vel = 1500.0 * ones(nz,nx); vel((round(nz+1)/2):end,:)=2500.0;
%    sx=(nx-1)/2*dx;sz=(nz-11)/2*dx; s=1.0./vel;
%    t=eikonal2d(s,sx,sz,dx);
%    figure;subplot(121);imagesc(vel);
%    subplot(122);imagesc(t);
%

[nz,nx]=size(s);

t=zeros(nz,nx); % claim the matrix of t for cacultation

% decide the closest node point to the source position
isx=round(sx/h)+1;isz=round(sz/h)+1; 

% calculate the traveltime around the source position up to 5X5
t=eik2d_around_source(t,s,sx,sz,nx,nz,isx,isz,h);
% expand the wave front
t=eik2d_square_expand(t,s,nx,nz,isx,isz,h);

% down sample the traveltime table
if nargin==7
    t=t(1:ds:end,1:ds:end);
end

end

function t=eik2d_around_source(t,s,sx,sz,nx,nz,isx,isz,h)

ix1=max(isx-2,1);ix2=min(isx+2,nx);
iz1=max(isz-2,1);iz2=min(isz+2,nz);

for i2=ix1:ix2
    for i1=iz1:iz2
        iix1=min(i2,isx);iix2=max(i2,isx);
        iiz1=min(i1,isz);iiz2=max(i1,isz);
        ns=(iix2-iix1+1)*(iiz2-iiz1+1);
        s_ave=sum(sum(s(iiz1:iiz2,iix1:iix2)))/ns;
        t(i1,i2)=(((i2-1)*h-sx)^2+((i1-1)*h-sz)^2)^0.5*s_ave;
    end
end

end

function t=eik2d_square_expand(t,s,nx,nz,isx,isz,h)

%  The order will go through "Top" to "Right" to "Bottom" to "Left"        
%
%                        L---T---T---T---T---T---R
%                        |                       |
%                        L   *---*---*---*---*   R
%                        |   |               |   |
%                        L   *               *   R
%                        |   |               |   |
%                        L   *               *   R
%                        |   |               |   |
%                        L   *               *   R
%                        |   |               |   |
%                        L   *---*---*---*---*   R
%                        |                       |
%                        L---B---B---B---B---B---B

for i=3:max(max(isx-1,isz-1),max(nx-isx,nz-isz))
    % deal with the top edge
    if (isz-i>=1)
       i_outer=isz-i;i_inner=i_outer+1;
       j1=max(1,isx-i+1);j2=min(isx+i-1,nx);
       s_index=isx-j1+1;
       t_inner=t(i_inner,j1:j2);
       s_inner=s(i_inner,j1:j2);
       s_outer=s(i_outer,j1:j2);
       t_outer=eik2d_trat(t_inner,s_inner,s_outer,s_index,h);
       t(i_outer,j1:j2)=t_outer;
    end
    % deal with the right edge
    if (isx+i <= nx)
       i_outer=isx+i;i_inner=i_outer-1;
       j1=max(1,isz-i);j2=min(isz+i-1,nz);
       s_index=isz-j1+1;
       t_inner=t(j1:j2,i_inner);
       s_inner=s(j1:j2,i_inner);
       s_outer=s(j1:j2,i_outer);
       t_outer=eik2d_trat(t_inner,s_inner,s_outer,s_index,h);
       t(j1:j2,i_outer)=t_outer;
    end
    % deal with the bottom edge
    if (isz+i<=nz)
       i_outer=isz+i;i_inner=i_outer-1;
       j1=max(1,isx-i+1);j2=min(isx+i,nx);
       s_index=isx-j1+1;
       t_inner=t(i_inner,j1:j2);
       s_inner=s(i_inner,j1:j2);
       s_outer=s(i_outer,j1:j2);
       t_outer=eik2d_trat(t_inner,s_inner,s_outer,s_index,h);
       t(i_outer,j1:j2)=t_outer;
    end
    % deal with the left edge
    if (isx-i>=1)
       i_outer=isx-i;i_inner=i_outer+1;
       j1=max(1,isz-i);j2=min(isz+i,nz);
       s_index=isz-j1+1;
       t_inner=t(j1:j2,i_inner);
       s_inner=s(j1:j2,i_inner);
       s_outer=s(j1:j2,i_outer);
       t_outer=eik2d_trat(t_inner,s_inner,s_outer,s_index,h);
       t(j1:j2,i_outer)=t_outer;
    end
end

end

function t_outer=eik2d_trat(t_inner,s_inner,s_outer,s_index,h)

n_array=length(t_inner);
t_outer=zeros(n_array,1);
[min_number,min_loc]=eik2d_local_min(t_inner);
[max_number,max_loc]=eik2d_local_max(t_inner);

% calculate the traveltime just outside the local minimum points
for i=1:min_number
    i_min=min_loc(i);
    if (i_min==1 || i_min==n_array)
        t_outer(i_min)=t_inner(i_min)+h*s_outer(i_min); 
        % question: solver it with 3 points instead of 1 point?
    else
        sq=(h*s_outer(i_min))^2-0.25*(t_inner(i_min-1)-t_inner(i_min+1))^2;
        if (sq < 0.0) 
            t_outer(i_min)=t_inner(i_min)+h*s_outer(i_min);
        else
            t_outer(i_min)=t_inner(i_min)+sqrt(sq);
        end
    end
end

% find the local minimum point closest to the source (this will divide the 
% whole eadge into two segments
if min_number==1
    i_close_min=1;
elseif min_number==2  
    % debug on Oct 1, 2012:
    % the matlab built-in function local_max cannnot give the local_max for
    % array [2,2], but here needs at least one.
    if abs(min_loc(1)-s_index)==abs(min_loc(2)-s_index)
        i_close_min=1;
    else
        i_close_min=local_max(-abs(min_loc(:)-s_index));
    end
else
    i_close_min=local_max(-abs(min_loc(:)-s_index));
end

% Left segment
% sweep from left to right, incoming wave
if (i_close_min==1) % no minimum left of the close point, no incoming wave
else
   for i=1:i_close_min-1
       i_min=min_loc(i);
       i_max=0;
       for j=1:max_number
          if (max_loc(j) > i_min ) 
            i_max=max_loc(j);
            break;
          end
       end
       for j=i_min+1:i_max
          s_ave=0.25*(s_inner(j)+s_inner(j-1)+s_outer(j)+s_outer(j-1));
          sq=2*(h*s_ave)^2-(t_inner(j)-t_outer(j-1))^2;
         if (sq<0)
            t_outer(j)=min(t_inner(j),t_outer(j-1))+h*s_ave;
         else
            t_outer(j)=t_inner(j-1)+sqrt(sq);
         end
       end
   end
end

% sweep from right to left, outgoint wave
if (min_loc(i_close_min)==1) % the close point is on left end, no outgoing)
else
   for i=i_close_min:-1:1
       i_min=min_loc(i);
       i_max=2*length(t_inner);
       for j=max_number:-1:1
           if ( max_loc(j) < i_min ) 
              i_max=max_loc(j);
              break;
           end
       end
       if (i_max ~= 2*length(t_inner)) 
          for j=i_min-1:-1:i_max
              s_ave=0.25*(s_inner(j)+s_inner(j+1)+s_outer(j)+s_outer(j+1));
              sq=2*(h*s_ave)^2-(t_inner(j)-t_outer(j+1))^2;
              if sq < 0.0 
                 t_outer(j)=min(t_inner(j),t_outer(j+1))+h*s_ave;
              else
                 t_outer(j)=t_inner(j+1)+sqrt(sq);
              end
          end
       end
   end
end

% Right segment
% Sweep from right to left,incoming wave
if (i_close_min==min_number) % no local min right to the close point,no incoming wave
else
   for i=min_number:-1:i_close_min+1
       i_min=min_loc(i);
       i_max=2*length(t_inner);
       for j=max_number:-1:1
           if (max_loc(j) < i_min) 
              i_max=max_loc(j);
              break;
           end
       end
       for j=i_min-1:-1:i_max
           s_ave=0.25*(s_inner(j)+s_inner(j+1)+s_outer(j)+s_outer(j+1));
           sq=2*(h*s_ave)^2-(t_inner(j)-t_outer(j+1))^2;
           if (sq < 0.0 )
              t_outer(j)=min(t_inner(j),t_outer(j+1))+h*s_ave;
           else
              t_outer(j)=t_inner(j+1)+sqrt(sq);
           end
       end
   end
end

% Sweep from left to right, outgoing wave
if (min_loc(i_close_min)==n_array) % the close min on the right end of array, no outgoing wave
else
   for i=i_close_min:min_number
       i_min=min_loc(i);
       i_max=0;
       for j=1:max_number
           if (max_loc(j) > i_min )
              i_max=max_loc(j);
              break;
           end
       end
       if (i_max ~= 0)
          for j=i_min+1:i_max
              s_ave=0.25*(s_inner(j)+s_inner(j-1)+s_outer(j)+s_outer(j-1));
              sq=2*(h*s_ave)^2-(t_inner(j)-t_outer(j-1))^2;
              if (sq < 0.0)
                 t_outer(j)=min(t_inner(j),t_outer(j-1))+h*s_ave;
              else
                 t_outer(j)=t_inner(j-1)+sqrt(sq);
              end
          end
       end
   end
end

end

function [min_number,min_loc]=eik2d_local_min(t_inner)

n_array=length(t_inner);
min_number=0;
min_loc=zeros(n_array,1);
% deal with the left end
if (t_inner(1)<=t_inner(2) ) 
   min_number=min_number+1;min_loc(min_number)=1;
end
% deal with the middle
for i=2:n_array-1
   if ( t_inner(i) <= t_inner(i-1) && t_inner(i) <= t_inner(i+1) )
       min_number=min_number+1;min_loc(min_number)=i;
   end
end
% deal with the right end
if ( t_inner(n_array) <= t_inner(n_array-1) )
    min_number=min_number+1;min_loc(min_number)=n_array;
end
min_loc=min_loc(1:min_number);

end

function [max_number,max_loc]=eik2d_local_max(t_inner)

n_array=length(t_inner);
max_number=0;
max_loc=zeros(n_array,1);
% deal with the left end
if (t_inner(1)>=t_inner(2) ) 
   max_number=max_number+1;max_loc(max_number)=1;
end
% deal with the middle
for i=2:n_array-1
   if ( t_inner(i) >= t_inner(i-1) && t_inner(i) >= t_inner(i+1) )
       max_number=max_number+1;max_loc(max_number)=i;
   end
end
% deal with the right end
if ( t_inner(n_array) >= t_inner(n_array-1) )
    max_number=max_number+1;max_loc(max_number)=n_array;
end
max_loc=max_loc(1:max_number);

end


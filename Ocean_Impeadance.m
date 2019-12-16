clear all
close all
clc

% Load month variation - 1D profiles (2002 - 2015)
load('Month_variation.mat')

% Sound velocity
figure;
contourf(X,Z,C,30)
xticks([1:12:totalFiles])
grid on
colorbar
xlabel('Month index')
ylabel('z [m]')
title('c [m/s]')

%Impedance
Imp = C.*RHO;

figure;
contourf(X,Z, Imp,30)
xticks([1:12:totalFiles])
grid on
colorbar
xlabel('Month index')
ylabel('z [m]')
title('Impedance [kg/m^2.s]')


% Reflection and Transmition Coefficients
for iN=1:length(C(:,1))-1
    
    Rcoef(iN,:) = (Imp(iN+1,:)-Imp(iN,:))./(Imp(iN+1,:)+Imp(iN,:));
    Tcoef(iN,:) = (2*Imp(iN,:))./(Imp(iN+1,:)+Imp(iN,:)); 
    
    % t = t0 + 2*dz/c;
    if iN==1
        tt(iN,1:totalFiles)=0;
    else
        tt(iN,:) = tt(iN-1,:) + 2.*abs(Z(iN+1,:)-Z(iN,:))./(C(iN+1,:)+C(iN,:))./2;
    end
end

figure;
subplot(1,2,1)
contourf(X(1:end-1,:),Z(1:end-1,:),Rcoef,30)
xticks([1:12:totalFiles])
grid on
colorbar
xlabel('Month index')
ylabel('z [m]')
title('R coefficient [ ]')

subplot(1,2,2)
contourf(X(1:end-1,:),Z(1:end-1,:),Tcoef,30)
xticks([1:12:totalFiles])
grid on
colorbar
xlabel('Month index')
ylabel('z [m]')
title('T coefficient [ ]')

% figure
for imonth=1:totalFiles
    
    % For 1D trace - ricker wavelet
    [w,tw]=ricker(45,tt(2,imonth),length(tt(:,1)));
    
%     plot(tw,w)
%     hold on
    
    trace(:,imonth) = conv(w,Rcoef(:,imonth));
    
end

% % select month to plot (1 to 168)
% imonth = 80;

gifname = '1D_trace.gif';

h = figure
for imonth=1:totalFiles
    
    subplot(2,3,[1 2 3])
    contourf(X(1:end-1,:),tt,C(1:end-1,:),30);
    xticks([1:12:totalFiles])
    hold on
    plot([imonth imonth],[0 tt(end,imonth)],'r--')
    grid on
    colorbar
    xlabel('Month index')
    ylabel('t [s]')
    title('c [m/s]')
    set(gca,'Ydir','reverse')

    subplot(2,3,4)
    plot(Rcoef(1:end,imonth),tt(1:end,imonth))
    ylim([0 tt(end,imonth)])
    xlabel('R coef.')
    ylabel('t=t_0 + 2*\Deltaz/c_{avg} [s]')
    set(gca,'Ydir','reverse')
    set(gca,'XAxisLocation','top')

    subplot(2,3,5)
    plot(w,tw)
    ylim([0 tt(end,imonth)])
    xlabel('Ricker wavelet')
    ylabel('t [s]')
    set(gca,'Ydir','reverse')
    set(gca,'XAxisLocation','top')
    title(['Reflection Coef. / Wavelet / Convolution (imonth = ', num2str(imonth), ')'])


    subplot(2,3,6)
    plot(trace( 1:length(tt(:,imonth)) , imonth), ...
         tt(1:end,imonth));
    ylim([0 tt(end,imonth)])
    xlabel('Convolution')
    ylabel('t [s]')
    set(gca,'Ydir','reverse')
    set(gca,'XAxisLocation','top')
    
    drawnow
%     % Capture the plot as an image 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if imonth == 1 
%       imwrite(imind,cm,gifname,'gif','DelayTime',0.3, 'Loopcount',inf); 
%     else 
%       imwrite(imind,cm,gifname,'gif','DelayTime',0.3,'WriteMode','append'); 
%     end 
    
end

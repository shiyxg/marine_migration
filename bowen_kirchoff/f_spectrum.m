% To display the frequency spectrum of a time series vector
% @version 1 2014-10-02
% @author Bowen Guo

function f_spectrum(v,dt)

[nt,ng]=size(v);

f=fft(v);

ff=1/dt/nt*(0:floor(nt/2)-1);

if (ng>1)
   figure;imagesc((1:ng),ff,log(abs(f(1:floor(nt/2),:))));
   xlabel('trace No.','fontsize',14);
   ylabel('frequency','fontsize',14);
else
   figure;plot(ff,abs(f(1:floor(nt/2))));
   xlabel('frequency','fontsize',14);
end


end




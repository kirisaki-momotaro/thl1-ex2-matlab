clear all;
close all;

F0=2; %rondomly decided cos frecuency

t=1:0.01:pi; %define time
%plot 5 instances of the given function
figure(24)
for i =1:5
   X= normrnd(0,1) %gaussian distribution (0,1)
   FI= 2*pi*rand(1) %uniform distribution (0,2pi)
   Yt=X*cos(2*pi*F0*t+FI); 
   
   plot(t,Yt)
   hold on;
end
title('5 instances of X*cos(2*pi*F0*t+FI)')

%power spectral density MAYBE WRONG!!!!!!!!!!!!!
K=500; %number of experiments
power_fftshift_signal_sum=zeros(1,2048);
for i = 1:K
    X= normrnd(0,1) %gaussian distribution (0,1)
    FI= 2*pi*rand(1) %uniform distribution (0,2pi)
    Yt=X*cos(2*pi*F0*t+FI); 

    %fft signal
   
    fftshift_signal = fftshift(fft(Yt,2048));
    power_fftshift_signal = abs(fftshift_signal).^2;     % zero-centered power
    power_fftshift_signal_sum=power_fftshift_signal_sum+power_fftshift_signal; %add up all waveforms
    
end
figure(25)
power_fftshift_signal_sum_normal=power_fftshift_signal_sum/K; %divide by K to normalize
semilogy(power_fftshift_signal_sum_normal)
grid on;



% save images
% FolderName = ('C:\Users\chris\Desktop\THL I\ex2\thl1-ex2-matlab\images');   % using my directory
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = num2str(get(FigHandle, 'Number'));
%   set(0, 'CurrentFigure', FigHandle);
% %   saveas(FigHandle, strcat(FigName, '.png'));
%   saveas(FigHandle, fullfile(FolderName,strcat(FigName, '.png'))); % specify the full path
% end



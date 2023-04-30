clear all;
close all;
F0=2;

t=1:0.01:2*pi;
figure()
for i =1:5
   X= normrnd(0,1)
   FI= 2*pi*rand(1)
   Yt=X*cos(2*pi*F0*t+FI);
   
   plot(t,Yt)
   hold on;
end
%% Plot data from E4

clc
clear all

% task 1
% load the data file
distData = importdata('distribution.data');
set(gcf,'renderer','painters','PaperPosition',[0 0 12 6]);
Size = size(distData);

%plot
figure(1);
clf
subplot(2,2,1)
hist(distData(:,1));
title('Uniform','fontsize',12);

%plot
subplot(2,2,2)
hist(distData(:,2));
title('Uniform','fontsize',12);

%plot
subplot(2,2,3)
hist(distData(:,3));
title('Gauss','fontsize',12);

%plot
subplot(2,2,4)
hist(distData(:,4));
title('Gauss','fontsize',12);

print(gcf,'-depsc2','distributions.eps')

%% Task 2
clc 
clear all

% load data file
trajData05 = importdata('trajectory05.data');
trajData5 = importdata('trajectory5.data');

% plot data
figure(2);
clf
plot(linspace(0,length(trajData05)*0.05,length(trajData05)),trajData05,'g', linspace(0,length(trajData5)*0.05,length(trajData5)),trajData5,'b');
hold on
axis([0 100 -0.1 0.1]);
xlabel('Time [s]','fontsize',12);
ylabel('Trajectory [m]','fontsize',12,'Interpreter','latex');

l = legend('Trajectory for $\eta = 0.05 \omega$','Trajectory for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','trajectory.eps')

%% Plot powerspectrum
clc
clear all

% load data file
trajData05 = importdata('trajectory05.data');

fftData05 = abs(fft(trajData05(1:1001)));
fftData05 = fftshift(fftData05);

for i = 101:1000:length(trajData05)-1000
    temp = abs(fft(trajData05(i:i+1000)));
    fftData05 = fftData05 + fftshift(temp);
end

fftData05 = fftData05./100;
freq05 = linspace(-pi/(0.05),pi/(0.05),length(fftData05));

% load data file
trajData5 = importdata('trajectory5.data');

fftData5 = abs(fft(trajData5(1:1001)));
fftData5 = fftshift(fftData5);

for i = 101:1000:length(trajData5)-1000
    temp = abs(fft(trajData5(i:i+1000)));
    fftData5 = fftData5 + fftshift(temp);
end

fftData5 = fftData5./100;
freq5 = linspace(-pi/(0.05),pi/(0.05),length(fftData5));

figure(3);
clf
plot(freq05, fftData05,'g',freq5, fftData5,'b');
xlim([-5 5]);

% labels
xlabel('Frequency [THz]','fontsize',12);
ylabel('Amplitude','fontsize',12);
title('Powerspectrum of the trajectory','fontsize',12);

l = legend('Powerspectrum for $\eta = 0.05 \omega$','Powerspectrum for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','powerspectrum.eps')



%% Corr func
clc

corrData5 = importdata('corrfunc5.data');
corrData05 = importdata('corrfunc05.data');

figure(5);
clf
plot(linspace(0,0.05*length(corrData05),length(corrData05)),corrData05,'g-',linspace(0,0.05*length(corrData5),length(corrData5)),corrData5,'b-');
hold on
%axis([0 25 0.0015 0.005]);
xlabel('Time lag [ps]','fontsize',12);
ylabel('Amplitude','fontsize',12);
title('Correlation function','fontsize',12);

l = legend('Time correlation function for $\eta = 0.05 \omega$','Time correlation function for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','corrFunc.eps')


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
%%
% plot data
figure(2);
clf
subplot(2,1,1)
plot(linspace(0,length(trajData05)*0.05,length(trajData05)),trajData05,'g', linspace(0,length(trajData5)*0.05,length(trajData5)),trajData5,'b');
axis([0 50 -0.01 0.01]);
xlabel('Time [s]','fontsize',12);
ylabel('Trajectory [m]','fontsize',12,'Interpreter','latex');

l = legend('Trajectory for $\eta = 0.05 \omega$','Trajectory for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','trajectory.eps')

subplot(2,1,2)
plot(linspace(0,length(trajData5)*0.05,length(trajData5)),trajData5,'b');
axis([10 50 -0.000055 0.000055]);
xlabel('Time [s]','fontsize',12);
ylabel('Trajectory [m]','fontsize',12,'Interpreter','latex');

%% Plot powerspectrum
clc

fftData05 = abs(fft(trajData05(1:1001)));
fftData05 = fftshift(fftData05);

for i = 1001:1000:length(trajData05)-1000
    temp = abs(fft(trajData05(i:i+1000)));
    fftData05 = fftData05 + fftshift(temp);
end

fftData05 = fftData05./100;
freq05 = linspace(-pi/(0.05),pi/(0.05),length(fftData05));


fftData5 = abs(fft(trajData5(1:1001)));
fftData5 = fftshift(fftData5);

for i = 1001:1000:length(trajData5)-1000
    temp = abs(fft(trajData5(i:i+1000)));
    fftData5 = fftData5 + fftshift(temp);
end

fftData5 = fftData5./100;
freq5 = linspace(-pi/(0.05),pi/(0.05),length(fftData5));

power05 = fftData05.^2./(length(trajData05)*0.05);
power5 = fftData5.^2./(length(trajData5)*0.05);

figure(3);
clf
plot(freq05, power05,'g',freq5, power5,'b');
xlim([-5 5]);

% labels
xlabel('Frequency [Hz]','fontsize',12);
ylabel('Amplitude','fontsize',12);
title('Powerspectrum of the trajectory','fontsize',12);

l = legend('Powerspectrum for $\eta = 0.05 \omega$','Powerspectrum for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','powerspectrum.eps')


%% x-mas func
clc

corrData05 = real(fftshift(ifft(power05)));
corrData5 = real(fftshift(ifft(power5)));

x05 = linspace(0,length(trajData05)*0.05,length(corrData05)/2);
x5 = linspace(0,length(trajData5)*0.05,length(corrData5)/2);

figure(5);
clf
plot( corrData05(ceil(length(corrData05)/2):end-1), x05,'Color',[0.0 1.0 0.0]);
hold on
plot(corrData5(ceil(length(corrData5)/2):end-1), x5,'Color',[0.8 0.5 0.3]);
%axis([0 5000 -0.0000000004 0.0000000004]);
ylabel('Time lag [s]','fontsize',12);
xlabel('Amplitude','fontsize',12);
title('X-mas function','fontsize',12);

l = legend('Time correlation function for $\eta = 0.05 \omega$','Time correlation function for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
%set(gcf,'renderer','painters','PaperPosition',[0 0 8.3 11.7]);
print(gcf,'-depsc2','xFunc.eps')

%%

clf

fftData05 = abs(fft(trajData05));
power05 = fftData05.^2;

data = ifft(power05);
plot(data,'g');
hold on

fftData5 = abs(fft(trajData5));
power5 = fftData5.^2;

data = ifft(power5);
plot(data,'b');
axis([0 10000 -0.0005 0.0005]);

ylabel('Time lag [s]','fontsize',12);
xlabel('Amplitude','fontsize',12);
title('Correlation function','fontsize',12);

l = legend('Time correlation function for $\eta = 0.05 \omega$','Time correlation function for $\eta = 5 \omega$');
set(l,'Interpreter','latex')
print(gcf,'-depsc2','xFunc.eps')

%%
clf
data = trajData05;
datasq = data.^2;


norm = mean(datasq) - mean(data)^2;
cov =  xcov(data, 500000);
cov =  cov/max(cov);
index = find(cov(500001:end) < exp(-2),1);
phi_s = cov(500+index);

sigma = sqrt(norm/(900000)*index)
plot(cov);
axis([500000 510000 -0.2 0.2]);
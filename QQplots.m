% Script to calculate intra-hourly wind data from hourly wind data. Based
% on 'main.m' written by Erik Jonasson and further developed by Samuel
% Forsberg.

%% Add input data here!
%The input data should be in vector format, and describe the intrahourly
%variability. For example, if "power" is 10 minute resolved power, and
% "powerhourly" is 60 minute resolved power (but the vectors of same
% length, such that powerhourly(1:6) is the same value, and is mean(power(1:6))
% then devs = power - powerhourly.

% Please make sure that there are no NaN-points. If it is, either replace
% these with sourrunding values if use devs = devs(~isnan(devs)) to select
% only the non-NaN points.

%% workspace cleanup
format compact
clear
clc


% Load data - NOTE: use one of the stations below (QPTR1 or NTKM3)
%  data2012_6min = load('data2012_G4_bus33_QPTR1_6min.mat');
%  data2012_1h = load('data2012_G4_bus33_QPTR1_1h.mat');

data2012_6min = load('data2012_G4_bus33_QPTR1_6min.mat');
data2012_1h = load('data2012_G4_bus33_QPTR1_1h.mat');

% Repeat elements
U2020_1h_rep = repelem(data2012_1h.U,10);

% Calculate deviation
devs = data2012_6min.U - U2020_1h_rep; % deviation between 6 min and 1 h resolution for each 6 min time step

% Introduce settings for Maximum Likelihood Algorithm
opt = statset('MaxIter',1e5);

%Use Maximum Likelihood estimation to generate parameters for
%t-scale-location distribuitno (student t) and gaussian
phat = mle(devs,'distribution','tLocationScale','options',opt);
phatgauss = mle(devs);

% "Manual" MLE estimation of parameters for La-place distribution
m = median(devs);
b = mean(abs(devs-m));

% Generates distribution objects of t-locationsclae and normal distribution
pd = makedist('tLocationScale','mu',phat(1),'sigma',phat(2),'nu',phat(3));
pd2 = makedist('normal','mu',phatgauss(1),'sigma',phatgauss(2));

% Just a for-loop to be able to plot the PDF of the laplace distribution
count = 1;
for i = -2:0.01:2
    f(count)=(1/(2*b))*exp(-(abs(i-m)/b));
    x(count)=i;
    count = count+1;
end


%% Plot the distribution and the binned data for visual inspection
figure('Name','Histogram')
h = histogram(devs,'Normalization','pdf'); %Binned data
hold on
plot(pd) %t-location scale
plot(pd2) % Gaussian
plot(x,f) % La-place

legend('Binned data','t-Location','Gaussian','Laplace')
hold off

%% Create simulated wind speed vectors
% from student-t distribution
U = U2020_1h_rep;
Ustudent = zeros(length(U),1);
noise_student = zeros(length(U),1);
for i = 1:length(U)
    noise_student(i) = random(pd);
    Ustudent(i) = U(i) + noise_student(i);
end

% from gaussian
Ugauss = zeros(length(U),1);
noise_gauss = zeros(length(U),1);
for i = 1:length(U)
    noise_gauss(i) = random(pd2);
    Ugauss(i) = U(i) + noise_gauss(i);
end

%from laplacian
% this loop calls 'randraw.m'
% more details: https://se.mathworks.com/matlabcentral/fileexchange/7309-randraw
% credit: Alex Bar-Guy, last update 06 Mar 2013

Ulaplace = zeros(length(U),1);
noise_laplace = randraw('laplace', [m,b], 1e5);
for i = 1:length(U)
    Ulaplace(i) = U(i) + noise_laplace(i);
end

% alternative laplace generation using laprnd.m
% I have commented this out because randraw gives better results
% % % Ulaplace = zeros(length(U),1);
% % % noise_laplace = laprnd(length(U),1,m, b);
% % % for i = 1:length(U)
% % %     Ulaplace(i) = U(i) + noise_laplace(i);
% % % end

%% QQ plots
% student distribution
Q1 = quantile(devs,(0.01:0.01:1));
Q2 = quantile(noise_student,(0.01:0.01:1));

x = linspace(-2,2,100);
figure(1)
%tiledlayout(3,1)
%nexttile
hold on
axis([-2 2 -2 2])
xlabel('Measured quantiles [m/s]') 
ylabel('Simulated quantiles [m/s]') 
plot(Q1,Q2,'+')
plot(x,x)
title('Student t')
hold off

% gaussian distribution
Q1 = quantile(devs,(0.01:0.01:1));
Q2 = quantile(noise_gauss,(0.01:0.01:1));

%nexttile
figure(2)
hold on
xlabel('Measured quantiles [m/s]') 
ylabel('Simulated quantiles [m/s]') 
axis([-2 2 -2 2])
plot(Q1,Q2,'+')
plot(x,x)
title('Gaussian')
hold off

% Laplacian distribution
Q1 = quantile(devs,(0.01:0.01:1));
Q2 = quantile(noise_laplace,(0.01:0.01:1));

%nexttile
figure(3)
hold on
xlabel('Measured quantiles [m/s]') 
ylabel('Simulated quantiles [m/s]') 
axis([-2 2 -2 2])
plot(Q1,Q2,'+')
plot(x,x)
title('Laplace')
hold off


% % %% The idea to generate higher frequency variations
% % targetpowh = % A vector with mean hourly values, and the size for the reference
% % % resolution. targetpow is the vector for the resulting 6-minute
% % % resolution, so size(targetpow) should be equal to size(targetpowh). Hope
% % % I made it clear :) 
% % 
% % %See below, uncomment the one that you want to use to generate the noise. 
% % 
% % for i = 1:length(targetpowh)
% %     temp = random(pd); %Generate a random number for the t-location scale distribution
% %     %temp = random(pd2); %Generate a random number from the gaussian distribution 
% %     %temp = laprnd(1,1,m,sqrt(2)*b); %Generate a random number from the la-placedistrbution
% %     while abs(temp + targetpowh(i)) > 1 || temp + targetpowh(i) < 0 % Constrains the power to 0-100%
% %         temp = random(pd);
% %         temp = random(pd);
% %         temp = laprnd(1,1,m,sqrt(2)*b);
% %     end
% %     targetpow(i) = targetpowh(i) + temp;
% % end

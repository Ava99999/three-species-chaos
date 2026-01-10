% system
param.a1 = 5;
param.a2 = 0.1;
param.b1 = 2.46; %2.55; %2.25; % varies between 2 and 6.2
param.b1_r = 2.2:0.001:3.2; %2.2:0.01:6.2; %2.25:0.001:2.6; %3.6:0.001:4; %2.25:0.001:2.6;  %2.9:0.01:4.2; 
param.b2 = 2; %1.6; %2;
param.b2_r = 1.5:0.001:2.15; %1.5:0.01:3.2;
param.d1 = 0.4; 
param.d2 = 0.01;
param.ic1 = [0.76, 0.16, 9.9]; % specific to reproduction figure 2
param.ic3 = [0.77, 0.17, 10];
param.ic_r = 10; % amount of ic's to try for bifurcation plots
%param.ic1 = [1, 1, 1];

% (reproduction paper)
% p = [5.0, 0.1, 3.0, 2.0, 0.4, 0.01] # values for a1, a2, b1, b2, d1 and d2 

% time steps
param.N = 10000; %10000; 
param.T = 10000; %10000;
param.dt = param.T/param.N;
param.t = [0:(param.N -1)]*param.dt;

% bifurcation plot
param.NB = 1000; %N_Bifurcation
param.NT = param.N/10; %N_01 test 

% 01 test
param.ic2 = [1,1,1];
param.sizec = 100;
help_c = pi/5 + 3*pi/5.*rand(1,param.sizec);
param.c = sort(help_c);


clear all
close all

chaos_parameters

%% defining three species system
f1 = @(u) param.a1*u/(1 + param.b1*u); % functional response f1
f2 = @(u) param.a2*u/(1 + param.b2*u); % functional response f2

%% parameters ode
tspan = param.t; %[t0,T];
ic = param.ic2; % initial condition for the 01-test
options = odeset('RelTol',1e-8,'abstol',1e-8*[1,1,1]);

%% solve ode
% todo might need to make the timesteps a bit smaller
% might need to sample just 10% in case we need to discard a lot (so
% integrate longer)
[t,X] = ode45(@(t,y) superpredators(t,y,f1,f2,param.d1,param.d2),tspan, ic, options);

% disregard transient behaviour and only select the 3rd variable
% functions to create p, q and M


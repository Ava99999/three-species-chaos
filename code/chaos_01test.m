clear all
close all

chaos_parameters

%% defining three species system
f1 = @(u) param.a1*u/(1 + param.b1*u); % functional response f1
f2 = @(u) param.a2*u/(1 + param.b2*u); % functional response f2

%% parameters ode
tspan = param.t; %[t0,T];
ic = param.ic1; % initial condition for the 01-test
options = odeset('RelTol',1e-8,'abstol',1e-8*[1,1,1]);

%% solve ode
% todo might need to make the timesteps a bit smaller
% might need to sample just 10% in case we need to discard a lot (so
% integrate longer)
[t,X] = ode45(@(t,y) superpredators(t,y,f1,f2,param.d1,param.d2),tspan, ic, options);

%% sample
%N = param.N - param.NT;
t_start = param.NT+1;
t_zrange = t(t_start:end);
t_end = 2000; % number of timesteps to plot
Z = X(t_start:end,3); % time series
indices = 1:20:length(Z); % to sample
phi = Z(indices);
mean_phi = mean(phi);
N = length(phi);
test_c = param.c;

%% create sequences
[p,q] = create_p_q(param.c,N,phi);
[M,N0] = create_M(param.c,p,q,N);
[D,V_osc] = create_D(param.c,M,N,mean_phi);
[K_seq,K] = correlation_method(N0,param.c,D); 

%% plots
% time series
figure;
plot(t_zrange(1:t_end),Z(1:t_end));
xlabel('t');
ylabel('z');
axis tight;
%ylim([7.2 10.7]);

% (p,q)-dynamics for a specific c
c_index = 10;
fprintf('Value of c used = %.6g\n', test_c(c_index));
figure;
%plot(p(c_index,1:500),q(c_index,1:500));
p_lot = plot(p(c_index,:),q(c_index,:));
p_lot.LineWidth = 0.1;
xlabel('p');
ylabel('q');
%xlim([-25 15]);
%ylim([-5 35]);
axis tight;

% comparing M and D for a specific c
figure;
plot(1:N0, M(c_index,:),'b-');
hold on;
plot(1:N0, D(c_index,:), 'r--');
plot(1:N0, V_osc(c_index,:), 'g');
hold off;
xlabel('n');
legend('M', 'D', 'V_{osc}');
axis tight;

% correlation coefficient
figure;
plot(param.c, K_seq);
hold on;
scatter(param.c, K_seq, "x", 'MarkerEdgeColor', [0 0.4470 0.7410]);
hold off;
xlabel('c');
ylabel('K_c');
title(sprintf('Median value = %.4g', K));

%% functions
function [p,q] = create_p_q(c,N,phi)
    % c is a row vector
    % phi is a time series of length N
    p = zeros(length(c),N);
    q = zeros(length(c),N);
    for i = 1:length(c)
        p(i,1) = phi(1)*cos(c(i));
        q(i,1) = phi(1)*sin(c(i));
        for n = 2:N
            p(i,n) = p(i,n-1) + phi(n)*cos(n*c(i));
            q(i,n) = q(i,n-1) + phi(n)*sin(n*c(i));
        end
    end
end

function [M,N0] = create_M(c,p,q,N)
    N0 = floor(N/10);
    M = zeros(length(c),N0);
    tic
    for i = 1:length(c)
        for n = 1:N0
            for j = 1:(N-N0)
                M(i,n) = M(i,n) + (p(i,(j+n))- p(i,j))^2 + (q(i,(j+n))-q(i,j))^2;
            end
            M(i,n) = (1/(N-N0))*M(i,n);
        end
    end
    toc
end

% D will probably be unused if it is regular enough
function [D, V_osc] = create_D(c,M,N,mean)
    N0 = floor(N/10);
    V_osc = zeros(length(c),N0);
    for i = 1:length(c)
        for n = 1:N0
            V_osc(i,n) = (mean^2)*((1-cos(n*c(i)))/(1 - cos(c(i))));
        end
    end
    D = M - V_osc;
end

function [K_seq, K] = correlation_method(N0, c, M)
    zeta = double(1:N0); 
    K_seq = zeros(1,length(c));
    for i = 1:length(c)
        R = corrcoef(zeta,M(i,:));
        K_seq(i) = R(1,2); 
    end
    K = median(K_seq);
end
    

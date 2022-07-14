close all;
clear

%% Parametes
Tm = 500;        % [h] Time horizon
M_b = 1e3;       % [-] Stories for bayesian integral loop
M_r = 1e3;       % [-] Stories for reliability MC sim

%% Components' transition rates
components(1,1) = 5e-3; %[1/h]
components(1,2) = 0.1;  %[1/h]
components(2,1) = 4e-3; %[1/h]
components(2,2) = 0.15; %[1/h]
components(3,1) = 3e-3; %[1/h]
components(3,2) = 0.05; %[1/h]
components(4,1) = 1e-3; %[1/h]
components(4,2) = 0.25; %[1/h]

var_coeff = 0.1;        %[-]
m_a = 5e-3;             %[1/h]
m_c = 3e-3;             %[1/h]

%% l_a & l_c distribution
v_a = (var_coeff*m_a)^2;
mu_a = log((m_a^2)/sqrt(v_a+m_a^2));
sigma_a = sqrt(log(v_a/(m_a^2)+1));

v_c = (var_coeff*m_c)^2;
mu_c = log((m_c^2)/sqrt(v_c+m_c^2));
sigma_c = sqrt(log(v_c/(m_c^2)+1));

%% Failure states
is_system_failed = @(state) state(4) || ((state(1)  && state(2)) || (state(2) && state(3)) || (state(1) && state(3)));

%% Monte Carlo simulations for reliability and availability

Rel_b = zeros(M_b,Tm+1);
Avail_b = zeros(M_b,Tm+1);

% Bayesian integral loop
for i=1:M_b
    % Sample a realization of components' transition rates
    l_a = lognrnd(mu_a, sigma_a);
    l_c = lognrnd(mu_c, sigma_c);

    % Save them into the transition table
    components(1,1) = l_a;
    components(3,1) = l_c;

    % MC simulation of the system
    [Time_axis, Rel_b(i,:), ~, Avail_b(i,:), ~] = mc_sim(components, is_system_failed, Tm, M_r, false);
end

Rel = mean(Rel_b,1);
Rel_var = (mean(Rel_b.^2, 1) - Rel.^2)/M_b;

Avail = mean(Avail_b,1);
Avail_var = (mean(Avail_b.^2, 1) - Avail.^2)/M_b;

%% Plots
figure(1)
hold on
title('Case study 3: Reliability')
xlabel('Time [h]');
ylabel('R(t)');
plot(Time_axis, Rel, 'k-')

figure(2)
hold on
title('Case study 3: Var[R(t)]')
xlabel('Time [h]');
ylabel('Var[R(t)]');
plot(Time_axis, Rel_var, 'k-')

figure(3)
hold on
title('Case study 3: Availability')
xlabel('Time [h]');
ylabel('A(t)');
plot(Time_axis, Avail, 'k-')

figure(4)
hold on
title('Case study 3: Var[A(t)]')
xlabel('Time [h]');
ylabel('Var[A(t)]');
plot(Time_axis, Avail_var, 'k-')

%% Display
fprintf('Case study 3:\n')
fprintf('  R from MC sim (95%%): (%f; %f)\n',...
    Rel(end)-2*sqrt(Rel_var(end)),...
    Rel(end)+2*sqrt(Rel_var(end)))
fprintf('  A from MC sim (95%%): (%f; %f)\n',...
    Avail(end)-2*sqrt(Avail_var(end)),...
    Avail(end)+2*sqrt(Avail_var(end)))
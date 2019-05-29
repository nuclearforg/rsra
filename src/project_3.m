close all;
clear

%% Parametes
Tm = 500;        % [h] Time horizon
M_e = 1e3;       % [-] Stories for epistemic uncertainty loop
M_a = 1e3;       % [-] Stories for aleatory uncertainty loop

%% Components' transition rates
components(1,1) = 5e-3; %[1/h]
components(1,2) = 0.1;  %[1/h]
components(2,1) = 4e-3; %[1/h]
components(2,2) = 0.15; %[1/h]
components(3,1) = 3e-3; %[1/h]
components(3,2) = 0.05; %[1/h]
components(4,1) = 1e-3; %[1/h]
components(4,2) = 0.25; %[1/h]

%% l_a & l_c distribution
m_a = 5e-3;
v_a = (0.1*m_a)^2;
mu_a = log((m_a^2)/sqrt(v_a+m_a^2));
sigma_a = sqrt(log(v_a/(m_a^2)+1));

m_c = 3e-3;
v_c = (0.1*m_c)^2;
mu_c = log((m_c^2)/sqrt(v_c+m_c^2));
sigma_c = sqrt(log(v_c/(m_c^2)+1));

%% Failure states
is_system_failed = @(state) state(4) || ((state(1)  && state(2)) || (state(2) && state(3)) || (state(1) && state(3)));

%% Monte Carlo simulations for reliability and availability

Rel_e = zeros(M_e,501);
Avail_e = zeros(M_e,501);

% External epistemic loop
for i=1:M_e
    % Sample a realization of the epistemic variables from the
    % distributions
    l_a = lognrnd(mu_a, sigma_a);
    l_c = lognrnd(mu_c, sigma_c);

    % Save them into the transition table
    components(1,1) = l_a;
    components(3,1) = l_c;

    % Aleatory uncertainty simulation
    [Time_axis, Rel_e(i,:), ~, Avail_e(i,:), ~] = mc_sim(components, is_system_failed, Tm, M_a, false);
end

Rel = mean(Rel_e,1);
Rel_var = mean(Rel_e.^2, 1) - Rel.^2;
Avail = mean(Avail_e,1);
Avail_var = mean(Avail_e.^2, 1) - Avail.^2;

%% Plots
figure(1)
hold on
R_plus_std = Rel + 2*sqrt(Rel_var);
R_minus_std = Rel - 2*sqrt(Rel_var);
fill([Time_axis,fliplr(Time_axis)], [R_plus_std,fliplr(R_minus_std)], 0.8*[1 1 1], 'EdgeColor','none')
plot(Time_axis, Rel, 'k-', 'LineWidth', 2)
plot(Time_axis,R_plus_std,'k-')
plot(Time_axis,R_minus_std,'k-')

figure(2)
hold on
A_plus_std = Avail + 2*sqrt(Avail_var);
A_minus_std = Avail - 2*sqrt(Avail_var);
fill([Time_axis,fliplr(Time_axis)], [A_plus_std,fliplr(A_minus_std)], 0.8*[1 1 1], 'EdgeColor','none')
plot(Time_axis, Avail, 'k-', 'LineWidth', 2)
plot(Time_axis,A_plus_std,'k-')
plot(Time_axis,A_minus_std,'k-')

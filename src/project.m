close all;
clear

bar = waitbar(0, "Initializing...");
%% Parametes
Tm = 500;        % [h] Time horizon
M = 1e5;         % [-] Stories for main MC computation
N_val = 1e3;     % [-] Trials for MTTF validation
M_val = 1e3;     % [-] Stories for MTTF validation

%% Components' transition rates
% Ex. 1
components1(1,1) = 5e-3; %[h]
components1(1,2) = 0;    %[h]
components1(2,1) = 4e-3; %[h]
components1(2,2) = 0;    %[h]
components1(3,1) = 3e-3; %[h]
components1(3,2) = 0;    %[h]
components1(4,1) = 1e-3; %[h]
components1(4,2) = 0;    %[h]

% Ex. 2
components2(1,1) = 5e-3; %[h]
components2(1,2) = 0.1;  %[h]
components2(2,1) = 4e-3; %[h]
components2(2,2) = 0.15; %[h]
components2(3,1) = 3e-3; %[h]
components2(3,2) = 0.05; %[h]
components2(4,1) = 1e-3; %[h]
components2(4,2) = 0.25; %[h]

%% Failure states
is_system_failed = @(state) state(4) || ((state(1)  && state(2)) || (state(2) && state(3)) || (state(1) && state(3)));

%% Analytical solutions (when applicable)

% Ex. 1
% Reliability of the 2 out of 3 subsystem
% R_23 = R_A*R_B*R_C + (1-R_A)*R_B*R_C + R_A*(1-R_B)*R_C + R_A*R_B*(1-R_C)
%      = R_A*R_B + R_B*R_C + R_A*R_C - 2*R_A*R_B*R_C
%
% Reliability of the entire system
% R_sys = R_D * (R_A*R_B*R_C + R)
%
% Symbolic calculus (MTTF is the integral of Reliability from 0 to +inf)
syms t
R1_real_sym = exp(t*-1e-3) .* (exp(t*-9e-3) + exp(t*-8e-3) + exp(t*-7e-3) - exp(t*-12e-3)*2);
MTTF1_real_sym = int(R1_real_sym, t, 0, inf);

% From symbolic to numerical
R1_real_t = matlabFunction(R1_real_sym);
MTTF1_real = double(MTTF1_real_sym);

%% Monte Carlo simulations for reliability and availability
% Ex. 1
bar = waitbar(0.01, bar, "Ex.1 MC sim in progress...");
[Time1, Rel1, Rel1_var, ~, ~] = mc_sim(components1, is_system_failed, Tm, M, true);
% Ex. 2
bar = waitbar(0.07, bar, "Ex.2 MC sim in progress...");
[Time2, Rel2, Rel2_var, Avail2, Avail2_var] = mc_sim(components2, is_system_failed, Tm, M, false);

%% MC syms for MTTF computation
bar = waitbar(0.35, bar, "MTTF sim...");
[MTTF1_MC, MTTF1_MC_var] = mttf(components1, is_system_failed, M);
[MTTF2_MC, MTTF2_MC_var] = mttf(components2, is_system_failed, M);

%% MC syms for MTTF validation
% In this matrix we'll store the result of each trial
MTTF_trials = zeros(N_val,2);

% For each trial we compute the MTTF and its std
bar = waitbar(0.36, bar, "MTTF validation...");
parfor i=1:N_val
    [MTTF_trial, MTTF_trial_var] = mttf(components1, is_system_failed, M_val);
    MTTF_trials(i,:) = [MTTF_trial, sqrt(MTTF_trial_var)];
end

MTTF_val_1sigma = sum(abs(MTTF_trials(:,1) - MTTF1_real) < MTTF_trials(:,2))/N_val;
MTTF_val_2sigma = sum(abs(MTTF_trials(:,1) - MTTF1_real) < 2*MTTF_trials(:,2))/N_val;
MTTF_val_3sigma = sum(abs(MTTF_trials(:,1) - MTTF1_real) < 3*MTTF_trials(:,2))/N_val;

close(bar)
%% Plots
figure(1)
hold on
plot(Time1, Rel1, 'k-');
fplot(R1_real_t, [0 Tm],'LineWidth',2);

figure(2)
plot(Time1,Rel1_var,'k-')

figure(3)
plot(Time2,Rel2,'k-')

figure(4)
plot(Time2,Rel2_var,'k-')

figure(5)
plot(Time2,Avail2,'k-')

figure(6)
plot(Time2,Avail2_var,'k-')

figure(7)
hold on
errorbar(MTTF_trials(1:20,1), 2*MTTF_trials(1:20,2), 'b.')
errorbar(MTTF_trials(1:20,1), MTTF_trials(1:20,2), 'r.')
plot([0 20], [MTTF1_real, MTTF1_real], 'k-')

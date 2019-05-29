close all;
clear

%% Parametes
Tm = 500;        % [h] Time horizon
M = 1e5;         % [-] Stories for main MC computation
N_val = 1e3;     % [-] Trials for MTTF validation
M_val = 1e3;     % [-] Stories for MTTF validation

%% Components' transition rates
% Ex. 1
components1(1,1) = 5e-3; %[1/h]
components1(1,2) = 0;    %[1/h]
components1(2,1) = 4e-3; %[1/h]
components1(2,2) = 0;    %[1/h]
components1(3,1) = 3e-3; %[1/h]
components1(3,2) = 0;    %[1/h]
components1(4,1) = 1e-3; %[1/h]
components1(4,2) = 0;    %[1/h]

% Ex. 2
components2(1,1) = 5e-3; %[1/h]
components2(1,2) = 0.1;  %[1/h]
components2(2,1) = 4e-3; %[1/h]
components2(2,2) = 0.15; %[1/h]
components2(3,1) = 3e-3; %[1/h]
components2(3,2) = 0.05; %[1/h]
components2(4,1) = 1e-3; %[1/h]
components2(4,2) = 0.25; %[1/h]

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
[Time1, Rel1, Rel1_var, ~, ~] = mc_sim(components1, is_system_failed, Tm, M, true);
% Ex. 2
[Time2, Rel2, Rel2_var, Avail2, Avail2_var] = mc_sim(components2, is_system_failed, Tm, M, false);

%% Average availability (integral mean over mission time span)
Avail_i = zeros(1, M);
for i=1:M
    Avail_i(i) = Avail2(round(rand()*Tm) + 1);
end
Avail2_avg = mean(Avail_i);
Avail2_avg_var = (mean(Avail_i.^2) - Avail2_avg^2)/M;

%% MC syms for MTTF computation
[MTTF1_MC, MTTF1_MC_var] = mttf(components1, is_system_failed, M);
[MTTF2_MC, MTTF2_MC_var] = mttf(components2, is_system_failed, M);

%% MC syms for MTTF validation
% In this matrix we'll store the result of each trial
MTTF_trials = zeros(N_val,2);

% For each trial we compute the MTTF and its std
for i=1:N_val
    [MTTF_trial, MTTF_trial_var] = mttf(components1, is_system_failed, M_val);
    MTTF_trials(i,:) = [MTTF_trial, sqrt(MTTF_trial_var)];
end
MTTF_val_1sigma = sum(abs(MTTF_trials(:,1) - MTTF1_real) < MTTF_trials(:,2))/N_val;
MTTF_val_2sigma = sum(abs(MTTF_trials(:,1) - MTTF1_real) < 2*MTTF_trials(:,2))/N_val;
MTTF_val_3sigma = sum(abs(MTTF_trials(:,1) - MTTF1_real) < 3*MTTF_trials(:,2))/N_val;

%% Reliability validation
% In this matrix we'll store the result of each trial
Rel_trials = zeros(N_val,2);

% For each trial we compute the MTTF and its std
for i=1:N_val
    [~, Rel_trial, Rel_trial_var]= mc_sim(components1, is_system_failed, Tm, M_val, true);
    Rel_trials(i,:) = [Rel_trial(end), sqrt(Rel_trial_var(end))];
end

Rel_val_1sigma = sum(abs(Rel_trials(:,1) - R1_real_t(Tm)) < Rel_trials(:,2))/N_val;
Rel_val_2sigma = sum(abs(Rel_trials(:,1) - R1_real_t(Tm)) < 2*Rel_trials(:,2))/N_val;
Rel_val_3sigma = sum(abs(Rel_trials(:,1) - R1_real_t(Tm)) < 3*Rel_trials(:,2))/N_val;

%% Plots
figure(1)
hold on
title('Case study 1: Reliability')
xlabel('Time [h]');
ylabel('R(t)');
plot(Time1, Rel1, 'k-');

figure(2)
hold on
title('Case study 1: Var[R(t)]')
xlabel('Time [h]');
ylabel('Var[R(t)]');
plot(Time1,Rel1_var,'k-')

figure(3)
hold on
title('Case study 2: Reliability')
xlabel('Time [h]');
ylabel('R(t)');
plot(Time2,Rel2,'k-')

figure(4)
hold on
title('Case study 2: Var[R(t)]')
xlabel('Time [h]');
ylabel('Var[R(t)]');
plot(Time2,Rel2_var,'k-')

figure(5)
hold on
title('Case study 2: Availability')
xlabel('Time [h]');
ylabel('A(t)');
plot(Time2,Avail2,'k-')

figure(6)
hold on
title('Case study 2: Var[A(t)]')
xlabel('Time [h]');
ylabel('Var[A(t)]');
plot(Time2,Avail2_var,'k-')

figure(7)
hold on
title('MTTF validation')
xlabel('Trial')
ylabel('MTTF')
errorbar(MTTF_trials(1:20,1), 2*MTTF_trials(1:20,2), 'b.')
errorbar(MTTF_trials(1:20,1), MTTF_trials(1:20,2), 'r.')
plot([0 20], [MTTF1_real, MTTF1_real], 'k-')
legend('\alpha = 0.955','\alpha = 0.683','MTTF_{true}')

figure(8)
hold on
title('Rel validation')
xlabel('Trial')
ylabel('R(Tm)')
errorbar(Rel_trials(1:20,1), 2*Rel_trials(1:20,2), 'b.')
errorbar(Rel_trials(1:20,1), Rel_trials(1:20,2), 'r.')
plot([0 20], [R1_real_t(Tm), R1_real_t(Tm)], 'k-')
legend('\alpha = 0.955','\alpha = 0.683','R_{true}(Tm)')
function [time_axis, rel, rel_var, avail, avail_var] = mc_sim(components, is_failing, Tm, M, stop_if_failing)
Dt=1;
n_comp = size(components,1);
time_axis=0:Dt:Tm;
counter_Ft=zeros(1,length(time_axis));
counter_q=zeros(1,length(time_axis));

for i=1:M
    t = 0;
    state = zeros(1,n_comp);
    prev_failed = false;
    has_failed = false;
    lower_b = 1;

    while t<Tm
        % Extract transition rates out of current state
        l_states = zeros(1,n_comp);
        for k=1:n_comp
            l_states(k) = components(k,state(k)+1);
        end
        l_state = sum(l_states);

        % Sample time
        t = t-log(1-rand())/l_state;

        % Find affected component
        p = [0, cumsum(l_states./l_state)];
        R = rand();
        affected = find(p<R, 1, 'last');

        % Toggle status
        state(affected) = ~state(affected);

        % Check if system is failing in new state
        if is_failing(state)
            lower_b=find(time_axis>=t,1,'first');

            prev_failed = true;

            if has_failed == false
                has_failed = true;
                counter_Ft(lower_b:end)=counter_Ft(lower_b:end)+1;
            end

            if stop_if_failing==true
                t=Tm;
            end

        elseif prev_failed==true
            upper_b=find(time_axis<t,1,'last');
            counter_q(lower_b:upper_b)=counter_q(lower_b:upper_b)+1;
            prev_failed = false;
        end
    end
end

unrel = counter_Ft/M;
rel = 1 - unrel;
rel_var = (unrel - unrel.^2)/M;
unav = counter_q/M;
avail = 1 - unav;
avail_var = (unav - unav.^2)/M;
end
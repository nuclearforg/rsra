function [MTTF, MTTF_var] = mttf(components, is_failing, M)
n_comp = size(components,1);
MTTF_i = zeros(1,M);

for i=1:M
    t = 0;
    state = zeros(1,n_comp);
    failed = false;

    while failed == false
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
            MTTF_i(i) = t;
            failed = true;
        end
    end
end

MTTF = mean(MTTF_i);
MTTF_var = (mean(MTTF_i.^2)-MTTF^2)/M;
end
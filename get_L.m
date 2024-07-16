function L = get_L(FRAME, PRIMARIES)

% -----------------------------------------------------------------%
% This function computes the mean distance of two given primaries over a
% timespan of 50 years 
% -----------------------------------------------------------------%

body_1 = PRIMARIES{1};
body_2 = PRIMARIES{2};
ti = '1950 JAN 2 TDB';
tf = '2049 DEC 31 TDB';

et_start = cspice_str2et(ti);
et_end = cspice_str2et(tf);

time_step_coarse = 86400; % in seconds
coarse_times = et_start:time_step_coarse:et_end;
num_coarse_times = length(coarse_times);
for i = 1:num_coarse_times
    [state, ~] = cspice_spkpos(body_2, coarse_times(i), FRAME, 'NONE', body_1);
    coarse_distances(i) = norm(state);
end

L = mean(coarse_distances);

% STEPS = 10000;
% times = linspace(start_time_et, end_time_et, STEPS);
% 
% [state, ~] = cspice_spkezr(body_2, times, FRAME, 'NONE', body_1); 
% pos = state(1:3);
% L = mean(norm(pos))

end
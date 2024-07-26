function get_residual_acceleration(rtbp_times, rtbp_states)
global mu
global PRIMARIES
global FRAME
global OBSERVER
global L
global BODIES
global interpolators

% --------------------------------------------------------------------------------------------------------------------------------%
% THIS Function Computes THE RESIDUAL ACCELERATION (AS IN Solar_system_models_with_a_selected_set_Gomez_Masdemont.pdf) 
% in the case of HTBP, RTBP and ETBP vs FF model. 

% Inputs are: rtbp_times and rtbp_states (as retrieved after a run of full_force model)
% Outputs are the residual accelerations for the three different models
% --------------------------------------------------------------------------------------------------------------------------------%


% --------------- LOAD KERNELS -------------
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %
% ---------------------------------- Compute trajectory length ------------------------------------------------------%
x_array = rtbp_states.'
D_matrix = norm(diff(x_array(1:6,:)));
D = sum(D_matrix);
% --------------------------------------------------------------------------------------------------------------------%

% ------------------Number of points to sample for the sum computation------------------------%

num_samples = 100;
% -----------------------------------------------------------------------------------------------------%

% ------------------------ Preallocate result vectors for the sampled points---------------------------%

res_acc_r3b = zeros(1, num_samples);
res_acc_h3b = zeros(1, num_samples);
res_acc_e3b = zeros(1, num_samples);
% -----------------------------------------------------------------------------------------------------%

% -----------------------------------------------------------------------------------------------------%
% Generate indices to sample points
sample_indices = round(linspace(1, length(rtbp_times), num_samples));
ff = zeros(6, num_samples);
r3b = zeros(6, num_samples);
h3b = zeros(6, num_samples);
e3b = zeros(6, num_samples);

t = rtbp_times(sample_indices);
x = x_array(:,sample_indices);
v = x_array(4:6,sample_indices);
% -----------------------------------------------------------------------------------------------------%


% ------------------- Computing Vector Fields at all points t, x along the  orbit  --------------------- %
parfor idx = 1:length(t)
    ff(:,idx) = full_force(t(idx,1),x(:, idx));
    r3b(:,idx) = rtbp(t(idx,1),x(:, idx));
    h3b(:,idx) = hill3b(t(idx,1),x(:, idx));
    e3b(:,idx) = etbp(t(idx,1),x(:, idx)); 
    res_acc_r3b(:,idx) = norm(ff(:, idx) - r3b(:, idx)) * norm(v(:, idx)) / (norm(ff(:,idx)) * D);
    res_acc_h3b(:,idx) = norm(ff(:, idx) - h3b(:, idx)) * norm(v(:, idx)) / (norm(ff(:,idx)) * D);
    res_acc_e3b(:,idx) = norm(ff(:, idx) - e3b(:, idx)) * norm(v(:, idx)) / (norm(ff(:,idx)) * D);
end
% -------------------------------------------------------------------------------------------------------%


% -------------- Computing the sum of residual acceleration over the curve ----------------------------%
sum(res_acc_r3b)
sum(res_acc_h3b)
sum(res_acc_e3b)
% -----------------------------------------------------------------------------------------------------%


%---------------------------------------Plotting Res_acc(t) -----------------------------------------------------%
figure
hold on

plot(t, res_acc_h3b, 'DisplayName', 'R3BP')
plot(t, res_acc_h3b, 'DisplayName', 'H3BP')
plot(t, res_acc_e3b, 'DisplayName', 'E3BP')

%ratio_r3b_h3b = sum(res_acc_r3b) / sum(res_acc_h3b)
%ratio_r3b_e3b = sum(res_acc_r3b) / sum(res_acc_e3b)

hold off
legend('show')
% -----------------------------------------------------------------------------------------------------%


cspice_kclear()

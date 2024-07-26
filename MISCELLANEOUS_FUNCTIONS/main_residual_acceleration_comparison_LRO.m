global mu
global PRIMARIES
global FRAME
global OBSERVER
global L
global BODIES
global interpolators

% ---- THIS MAIN COMPUTES THE RESIDUAL ACCELERATION (AS IN Solar_system_models_with_a_selected_set_Gomez_Masdemont.pdf) --- %

% In the case of HTBP and RTBP vs FF model, and for the LRO (SPICE-retrieved) POSITIONS AND VELOCITIES over the spacecraft's trajectory. 

% The code also ouputs the real LRO orbit (from SPICE) versus the one computed starting from the SPCIE initial conditons.
% --------------------------------------------------------------------------------------------------------------------------%



% --------------- LOAD KERNELS -------------
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %

% ----- Define the parameters to enter the convert_object conversion function -----
ti = '2009 SEP 17 TDB';
tf = '2009 SEP 18 TDB';
FRAME = 'J2000';
OBSERVER = 'SOLAR SYSTEM BARYCENTER';
BODIES = [{'MARS'}, {'SUN'}, {'JUPITER BARYCENTER'}]
PRIMARIES = [{'EARTH'}, {'MOON'}];
L = get_L(FRAME, PRIMARIES);

hmin=1.e-6; hmax=1.e0; tol=1.e-10; h=0.001;
ntall=1; iorb=-1;
% computes LRO RTBP-LIKE coordinates
[rtbp_times, rtbp_state_primaries, body_1, body_2, object, rtbp_state_object, interpolators_obj] = convert_object(ti,tf, 'J2000', 'SOLAR SYSTEM BARYCENTER', [{'EARTH'}, {'MOON'}], 384400, {'LRO'});

x_array = rtbp_state_object.("LRO");
x_body = x_array(1:6,:);
% -----------------------------------------------------------------------------------------------------%


% ---------------------------------- Compute trajectory length ------------------------------------------------------%

D_matrix = norm(diff(x_array(1:6,:)));
D = sum(D_matrix);
% --------------------------------------------------------------------------------------------------------------------%


% --------------------   DEFINE BODIES THAT NEED TO ENTER THE FF COMPUTATION ---------------------------%
% 

[rtbp_times, rtbp_state_primaries, body_1, body_2, BODIES, rtbp_state_bodies, interpolators] = new_spice_to_rtbp(ti,tf, 'J2000', 'SOLAR SYSTEM BARYCENTER', [{'EARTH'}, {'MOON'}], 384400, BODIES);
% -----------------------------------------------------------------------------------------------------%

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
etb = zeros(6, num_samples);

t = rtbp_times(sample_indices);
x = x_array(:,sample_indices);
v = x_array(4:6,sample_indices);
% -----------------------------------------------------------------------------------------------------%


% ------------------- Computing Vector Fields at all points t, x along the  orbit  --------------------- %
parfor idx = 1:length(t)
    ff(:,idx) = full_force(t(idx,1),x(:, idx));
    r3b(:,idx) = rtbp(t(idx,1),x(:, idx))
    h3b(:,idx) = hill3b(t(idx,1),x(:, idx))
    e3b(:,idx) = etbp(t(idx,1),x(:, idx))    
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

ratio_r3b_h3b = sum(res_acc_r3b) / sum(res_acc_h3b)
ratio_r3b_e3b = sum(res_acc_r3b) / sum(res_acc_e3b)

hold off
legend('show')
% -----------------------------------------------------------------------------------------------------%


% ------------------------------ Plotting Real LRO vs propagated one -------------------------- %

x_lro_init = [x_array(1,1), x_array(2,1), x_array(3,1), x_array(4,1), x_array(5,1), x_array(6,1)];
[tf_lro, xf_lro] = ode45(@full_force, [rtbp_times(1), rtbp_times(end)], x_lro_init);

figure 
hold on

fprintf('Initial LRO conditions in the Synodic System: %f\n', x_lro_init)
plot3(xf_lro(:,1), xf_lro(:,2), xf_lro(:,3))
plot3(x_body(1,:), x_body(2,:), x_body(3,:))
% -----------------------------------------------------------------------------------------------------%



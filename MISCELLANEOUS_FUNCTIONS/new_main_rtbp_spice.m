global e
global nu_dot
global interpolators

% --------------- LOAD KERNELS -------------
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %

FRAME = 'ECLIPJ2000';
OBSERVER = 'EARTH-MOON BARYCENTER';
BODIES = [{'MARS'}, {'VENUS'}, {'MERCURY'}, {'JUPITER BARYCENTER'}, {'SATURN BARYCENTER'}];
PRIMARIES = [{'EARTH'}, {'MOON'}];
L = get_L(FRAME, PRIMARIES);

mu = get_mu(PRIMARIES);
%xi_L1 = [(mu-1) + 0.001, 0,0,0.01,0.01,0];
xi_L1 = [(mu-1)+sqrt(mu/3), 0, 0, 0.02, 0.02, 0.02];
hmin=10^-6; hmax=1; tol=1.e-10; h=0.001;
ntall=1; iorb=-1;
rtbp_ti = 0;
rtbp_tf = 2*pi;

MODEL = '@rtbp';
% ------------------------------------------------------------ %
acc_flag = 0; %if acc_flag = 1 computes also the numerical acceleration, besides the analytical one
% ------------------------------------------------------------ %

% ---- Compute RTBP orbit in RTBP coordinates (synodic)  ------- %
if strcmp(MODEL, '@etbp') 
    BASE_MODEL = '@rtbp';
else
    BASE_MODEL = MODEL;
end

[t, x] = ode78(BASE_MODEL, [rtbp_ti, rtbp_tf], xi_L1);

STEPS = length(t);

ti = t(1);
tf = t(end);
n_rev = tf/(2*pi);

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nAdimensional Start Simulation Epoch: %f\n', ti)
fprintf('Adimensional End Simulation Epoch: %f\n\n', tf)
fprintf('This corresponds to approximately %.1f revolutions of the system of two primaries.\n', n_rev)

rtbp_times = linspace(ti, tf, STEPS);

rtbp_pos = x(:,1:3).';
rtbp_vel = x(:,4:6).';
rtbp_acc = zeros(3, STEPS);

for i = 1:STEPS
    rtbp_t = rtbp_times(i);
    rtbp_p = rtbp_pos(:,i);
    rtbp_v = rtbp_vel(:,i);
    results = rtbp(rtbp_t, [rtbp_p, rtbp_v]);
    rtbp_acc(:,i) = results(4:6,:);
end

%-------- Conversion from adimensional to dimensional time epochs ---------------
n = get_n(PRIMARIES, L); %L is the distance between the PRIMARIES (averaged over 50 years)
times = (rtbp_times(:).')/n;
% -----------------------------------------------------------------%

%--------- Retrieve INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
ti_utc = cspice_et2utc(times(1), 'C', 0);
tf_utc = cspice_et2utc(times(end), 'C', 0);
cspice_kclear()

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nUTC Start Simulation Epoch: %s\n', ti_utc)
fprintf('UTC End Simulation Epoch: %s\n', tf_utc)
fprintf('-------------------------------------------------------------------------------------------------------------------\n')
fprintf('MAIN\nCalling get_ephemeris function...\n\n')


META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
[inertial_state_primaries, ~, interpolators] = get_ephemeris(times, PRIMARIES, BODIES, FRAME, OBSERVER);

% -----------------------------------------------------------------%
%This part takes care of the ETBP model. In that case, the adimensional
%time is rtbp_times/nu_dot. Therefore, we compute an average nu_dot and use
%the formula rtbp_times/nu_dot = times*n/nu_dot = rtbp_times*n/(n*nu_dot)
if strcmp(MODEL, '@etbp')
    avg_e = get_avg_e(inertial_state_primaries, PRIMARIES);
    avg_nu_dot = get_avg_nu_dot(avg_e, inertial_state_primaries, PRIMARIES, L);

    fprintf('Compare avg_nu_dot with n:\n')
    fprintf('End Comparison.\n')
    times = times*n/avg_nu_dot;
end
% -----------------------------------------------------------------%


% ------------------------ RTBP variables ---------------------------%
[mu, body_1, body_2] = get_mu(PRIMARIES);
GM_1 = get_GM_body(body_1);
GM_2 = get_GM_body(body_2);
%---------------------------------------------------------------------%


%-----------------------------------------------------------------------------------------------------%

% --------------------- GO FROM SYNODIC (RTBP) TO INERTIAL REFERENCE SYSTEM (FOR THE SPACECRAFT) -------------------------

PRIMARY_1 = PRIMARIES{1};
PRIMARY_2 = PRIMARIES{2};
PRIMARY_str_1 = regexprep(PRIMARY_1, [{'\s+'}, {'-'}], '_');
PRIMARY_str_2 = regexprep(PRIMARY_2, [{'\s+'}, {'-'}], '_');

pos_1 =inertial_state_primaries.(PRIMARY_str_1).position;
pos_2 =inertial_state_primaries.(PRIMARY_str_2).position;

rel_pos = pos_2 - pos_1;
vel_1 =inertial_state_primaries.(PRIMARY_str_1).velocity;
vel_2 =inertial_state_primaries.(PRIMARY_str_2).velocity;
rel_vel = vel_2 - vel_1;
rel_acc = zeros(3,STEPS);
over_acc = zeros(3,STEPS);
acc_1 = zeros(3,STEPS);

SEb_pos = pos_1 + mu*rel_pos;
SEb_vel = vel_1 + mu*rel_vel;

% Compute the Relative Accelerations (needed for velocity conversion)
rel_acc(1,1:end-1) = diff(rel_vel(1,:))./diff(times);
rel_acc(2,1:end-1) = diff(rel_vel(2,:))./diff(times);
rel_acc(3,1:end-1) = diff(rel_vel(3,:))./diff(times);
rel_acc(:,1) = rel_acc(:,2);
rel_acc(:,end) = rel_acc(:,end-1);


% Compute the Relative Over Accelerations (needed for acceleration conversion)
over_acc(1,1:end-1) = diff(rel_acc(1,:))./diff(times);
over_acc(2,1:end-1) = diff(rel_acc(2,:))./diff(times);
over_acc(3,1:end-1) = diff(rel_acc(3,:))./diff(times);
over_acc(:,1) = over_acc(:,2);
over_acc(:,end) = over_acc(:,end-1);


%Initialize Inertial Positions and Velocities of spacecraft
check_rtbp_pos = zeros(3, length(rel_pos));
check_rtbp_vel = zeros(3, length(rel_pos));

% ----------------- Acceleration of Barycenter is needed. Here we compute it ----------%
acc_1(1,1:end-1) = diff(vel_1(1,:))./diff(times);
acc_1(2,1:end-1) = diff(vel_1(2,:))./diff(times);
acc_1(3,1:end-1) = diff(vel_1(3,:))./diff(times);
acc_1(:,1) = acc_1(:,2);
acc_1(:,end) = acc_1(:,end-1);
SEb_acc = acc_1 + mu*rel_acc;



% ----- Convert the Computed Orbit's (pos, vel) into the inertial reference frame ----- %
[inertial_pos_spacecraft, inertial_vel_spacecraft, inertial_acc_spacecraft] = go_inertial(rel_pos,rel_vel, rel_acc, over_acc, SEb_pos, SEb_vel, SEb_acc, rtbp_pos, rtbp_vel, rtbp_acc, n);
% -------------------------------------------------------------------------------------------%
%Concatenate Inertial Velocities and Accelerations to get the full vector field (xdot) in output
inertial_vf_spacecraft = [inertial_vel_spacecraft; inertial_acc_spacecraft];

%Same for positions and velocities, to get the full state
inertial_state_spacecraft = [inertial_pos_spacecraft; inertial_vel_spacecraft];


% ----- Compute the full_force vector field values along the RTBP trajectory in the inertial system ----- %

% Generate indices to sample points
num_samples = 100;
sample_indices = round(linspace(1, length(t), num_samples));

sampled_times = times(sample_indices); % sampled times in the inertial physical frame
sampled_states = inertial_state_spacecraft(:,sample_indices); %sampled rtbp states in the inertial physical frame

full_force_xdot = zeros(6, length(sampled_times));

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nComputing Full Ephemeris Model...\n\n');
fprintf('Function: full_force\nComputing all (PRIMARIES + BODIES) gravitational contributions on the spacecraft...\n')
parfor s = 1:num_samples
    full_force_xdot(:,s) = full_force(sampled_times(s),sampled_states(:,s)); %full_force vector field physical values 
end
fprintf('Done.\n')
fprintf('-----------------------------------------------------------\n')

% ------------------------------------------------------------ %

% ------ Compute the Residual accelerations along the trajectory as in Solar_system_models_with_a_selected_set_Gomez_Masdemont.pdf ------- %

%1) Get difference in acceleration between the two models

rtbp_xdot = inertial_vf_spacecraft(:,sample_indices); %sampled rtbp vector field converted in inertial physical coordinates

res_state = rtbp_xdot-full_force_xdot;

%2) Get length of trajectory
D_matrix = norm(diff(inertial_state_spacecraft(1:3,sample_indices)));
D = sum(D_matrix);

%3) Implement TRAPEZOIDAL RULE as in wikipedia, with deltaxk = delta_t, the
%differences in inertial times  (this is because, given t, we associate a
%point on the curve gamma. this means we have the gamma as a function of t,
%and therefore our f(x) depends ultimately on t only. 

fx = zeros(1, length(sample_indices)); %function to be integrated 
for i = 1:length(sample_indices)
fx(1,i) = norm(res_state(4:6,i))*(norm(rtbp_xdot(1:3,i)))/norm(rtbp_xdot(4:6,i));
end

delta_times = diff(times);
res_acceleration = 0;
for i = 1:length(sample_indices)-1
    tmp = 0.5*(fx(i) + fx(i+1))*delta_times(i)/D; %trapezoidal sum 
    res_acceleration = res_acceleration + tmp; 
end

% --------------------------------------------------------------------- %

% ----- Plot the spacecraft's trajectory in the inertial reference system ----- %
figure
%%hold on
plot3(inertial_state_spacecraft(1 ,:),inertial_state_spacecraft(2,:),inertial_state_spacecraft(3,:),'DisplayName', 'Spacecraft Orbit')
%plot(inertial_state_spacecraft_3(1 ,:),inertial_state_spacecraft_3(2,:), 'DisplayName', 'Spacecraft Orbit 3')


%plot3(inertial_state_primaries.('MOON' ).position(1, :),inertial_state_primaries.('MOON' ).position(2, :), inertial_state_primaries.('MOON' ).position(3, :), 'DisplayName', 'EM Barycenter Orbit')

%plot3(inertial_state_primaries.('EARTH').position(1, :),inertial_state_primaries.('EARTH').position(2, :),inertial_state_primaries.('EARTH').position(3, :), 'DisplayName', 'SUN Orbit')
legend('show')
% ------------------------------------------------------------ %
% Going synodic again to make sure that the the go_synodic function works well 
% and gives back the original [t,x] of the rtbp simulation

inertial_pos = inertial_state_spacecraft(1:3,:);
inertial_vel = inertial_state_spacecraft(4:6,:);
inertial_acc = inertial_vf_spacecraft(4:6,:);

[rtbp_pos_spacecraft, rtbp_vel_spacecraft, rtbp_acc_spacecraft] = go_synodic(rel_pos, rel_vel, rel_acc, over_acc, SEb_pos, SEb_vel, SEb_acc, inertial_pos, inertial_vel, inertial_acc, n);


%--------------------------------------------------------------------------------------------%
% Compute the Relative Accelerations (needed for velocity conversion)
if acc_flag == 1 %If flag is on, initialize and compute the numerical acceleration of the spacecraft
fprintf('-------------------------------------------------------------------------------------------------------------------\n')
fprintf('MAIN\nRetrieving inertial numerical acceleration for the spacecraft\n')
inertial_acc_spacecraft_numerical = zeros(3, length(inertial_acc_spacecraft));
inertial_acc_spacecraft_numerical (1,1:end-1) = diff(inertial_vel_spacecraft(1,:))./diff(times);
inertial_acc_spacecraft_numerical (2,1:end-1) = diff(inertial_vel_spacecraft(2,:))./diff(times);
inertial_acc_spacecraft_numerical (3,1:end-1) = diff(inertial_vel_spacecraft(3,:))./diff(times);
inertial_acc_spacecraft_numerical (:,1) = inertial_acc_spacecraft_numerical (:,2);
inertial_acc_spacecraft_numerical (:,end) = inertial_acc_spacecraft_numerical (:,end-1);
fprintf('Retrieved.\n')


% Here, first we compute the absolute value of the differences between
% numerical vs computed accelerations (so we get a 3 X STEPS matrix)
fprintf('-------------------------------------------------------------------------------------------------------------------\n')
fprintf('Function: rtbp_to_spice\nComparing the numerical acceleration versus the computed acceleration as a check...\n')
abs_numerical_vs_computed = abs(inertial_acc_spacecraft_numerical(:,:)-inertial_acc_spacecraft(:,:));
%Then we compute the mean value of the difference for every coordinate
%x,y,z so as to obtain a mean 3 X 1 vector
mean_numerical_vs_computed = zeros(3,1);
mean_numerical_vs_computed(1) = mean(abs_numerical_vs_computed(1,:));
mean_numerical_vs_computed(2) = mean(abs_numerical_vs_computed(2,:));
mean_numerical_vs_computed(3) = mean(abs_numerical_vs_computed(3,:));

%In the end, we take the norm of this vector
norm_numerical_vs_computed = norm(mean_numerical_vs_computed)*1000; %1000 is km --> m conversion factor

fprintf('The mean difference between numerical acceleration and computed acceleration is %.e m/s^2 .\n', norm_numerical_vs_computed);
end
%--------------------------------------------------------------------------------------------%


% Perform a quick check that the obtained inertial positions and velocities
% can be convertedback into the original rtbp positions and velocities and
% they are the same if converted back to rtbp coordinates

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nChecking whether the backward conversion INERTIAL --> SYNODIC has been successfull...\n')
rtbp_state_spacecraft = [rtbp_pos_spacecraft; rtbp_vel_spacecraft];
array_to_check = rtbp_state_spacecraft - x.';
check = 0;
for i = 1:3
    for j = 1:length(array_to_check)
       if array_to_check(i,j) > 1e-10
           check = 1 ;
       end
    end
end

if check == 1
    error('Attention. Conversion check has not been successfull. Aborting...\n')
    fprintf('-----------------------------------------------------------\n')
else
    fprintf('All good.\n')
    fprintf('-----------------------------------------------------------\n\n')

    fprintf('The residual acceleration over %d sampled points for the model %s is: %.10e\n\n', num_samples, MODEL, res_acceleration)
end

% --------------- %
cspice_kclear()
% -------------- %
%function [f_value, correction] = parallel_shooting(N, orbit_file)

clear all; close all; clc;

global inertial_state_primaries
global interpolators
global n_rtbp
global n_anomalistic
global PRIMARIES
global mu
global BODIES
global L
global eclipse_date_et
global FRAME
global OBSERVER

% --------------- LOAD KERNELS -------------
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %

%  --------------- SET THE MODEL -------------
FRAME = 'ECLIPJ2000';
OBSERVER = 'EARTH-MOON BARYCENTER';
BODIES = [];
PRIMARIES = [{'EARTH'}, {'MOON'}];
L = get_L(FRAME, PRIMARIES); % Computes the mean distance of two given primaries over a timespan of 50 years 
%L = 384601.25606767; %value form Gomez et al
mu = get_mu(PRIMARIES); %Compute grav. parameter for the system
MODEL = '@etbp'; 
%-------- Compute mean anomaly (global variable, needed for time converison) ---------------%
n_anomalistic = 2.639394888546285e-06; %computed with the function n_peri_peri, corresponding to the anomalistic month
n_rtbp = get_n(PRIMARIES,L); %computed with Kepler's 3rd Law
%------------------------------------------------------------------------------%

%-----Set LUNAR ECLIPSE of 21 JAN 2000 as starting point in time (taken from almanacs and verified with a SPICE plot) -----%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
eclipse_date_UTC= '21 JAN 2000 04:05:02';
eclipse_date_et = cspice_str2et(eclipse_date_UTC);
t_list = linspace(eclipse_date_et,eclipse_date_et + 100*pi/n_rtbp, 27*100); %propagation times for get_ephemeris function.
cspice_kclear()
%-------------------------------------------------------------------------------------------------------------%

%--------- Retrieve INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
[inertial_state_primaries, inertial_state_bodies, interpolators] = get_ephemeris(t_list, PRIMARIES, BODIES, FRAME, OBSERVER);
N = 20; %number of nodes

mu = 0.0122;
%orbit_file = 'Halo.txt';
orbit_file = 'orbit_coordinates.txt';
%mu = 0.1 %Masde initial conditions
%orbit_file = 'orbit_data.txt' %Masde orbits

% Retrieve (read) rtbp orbit
[t, x] = read_orbit(orbit_file);
%rtbp times
ti = t(1); % 
tf = t(end);  

% figure
% plot3(x(:,1), x(:,2), x(:,3))

%convert rtbp times into inertial times (eclipse date = ti)
ti_conversion = (ti/n_anomalistic + eclipse_date_et);
tf_conversion = (tf/n_anomalistic + eclipse_date_et);

%sample the rtbp times
t_sampled = linspace(ti, tf, N);
%sample the inertial times
t_sampled_inertial = linspace(ti_conversion, tf_conversion, N);

%interpolate the rtbp orbit along rtbp times
x_sampled = zeros(6, length(t_sampled));
for dim = 1:6
   interpol.spline{dim} = spline(t, x(:,dim).');
   x_sampled(dim,:) = ppval(interpol.spline{dim}, t_sampled);
end

% Evaluate the splines for the primaries
primary_name = PRIMARIES{1};
secondary_name = PRIMARIES{2};
primary_str = regexprep(primary_name, [{'\s+'}, {'-'}], '_');
secondary_str = regexprep(secondary_name, [{'\s+'}, {'-'}], '_');

% Initialize arrays to store interpolated position and velocity
rp = zeros(3, N);
vp = zeros(3, N);
rs = zeros(3, N);
vs = zeros(3, N);
R = zeros(N,3);
V = zeros(N,3);

for dim = 1:12
    t_sampled;
    interpolated_p = ppval(interpolators.(primary_str).spline{dim}, t_sampled_inertial);
    interpolated_s = ppval(interpolators.(secondary_str).spline{dim}, t_sampled_inertial);
    if dim <= 3
        rp(dim, :) = interpolated_p;
        rs(dim,:) = interpolated_s;
    elseif dim>=4 && dim<=6 
        vp(dim-3, :) = interpolated_p;
        vs(dim-3, :) = interpolated_s;
    elseif dim>=7 && dim<=9
        ap(dim-6, :) = interpolated_p;
        as(dim-6, :) = interpolated_s; 
    else
        oap(dim-9, :) = interpolated_p;
        oas(dim-9, :) = interpolated_s;            
    end
end
%retrieving primaries relative acceleration at each rtbp time
rs_rp = rp-rs;
vs_vp = vp-vs;
as_ap = ap-as;
oas_oap = oap-oas;


%barycenter
SEb_pos = rp - mu*rs_rp;
SEb_vel = vp - mu*vs_vp;
SEb_acc = ap - mu*as_ap;

rtbp_pos = x_sampled(1:3,:);
rtbp_vel = x_sampled(4:6,:);
rtbp_acc = zeros(3, N);

for i = 1:N
    rtbp_t = t_sampled(i);
    rtbp_p = rtbp_pos(:,i);
    rtbp_v = rtbp_vel(:,i);
    results = rtbp(rtbp_t, [rtbp_p, rtbp_v]);
    rtbp_acc(:,i) = results(4:6,:);
end

[inertial_pos_spacecraft, inertial_vel_spacecraft, ~] = go_inertial(rs_rp, vs_vp, as_ap, oas_oap, SEb_pos, SEb_vel, SEb_acc, rtbp_pos, rtbp_vel, rtbp_acc, n_rtbp);

Q0 = [inertial_pos_spacecraft;inertial_vel_spacecraft];
phi_Q_list = []
for iteration = 1:15
fprintf('iteration %f\n', iteration)
t_list_ = [];
F_list = [];
df = cell(N-1, N);

fprintf('Q0 to recover %f\n', inertial_pos_spacecraft(1,1))
fprintf('Q0 old %f\n', Q0(1,1))


for i = 1:N-1
    xiv=eye(6,6);

    [t_, phi_Q_tot] = ode78(@new_full_force_var_vectorized, [t_sampled_inertial(i), t_sampled_inertial(i+1)], [Q0(:,i);xiv(:)]);
    phi_Q = phi_Q_tot(:,1:6);

    if iteration == 15
        phi_Q_list = [phi_Q_list; phi_Q];
    end
    phi_Q_var = phi_Q_tot(:,7:42);
    stm_x_noised = phi_Q_var(end,1:36);
    stm_6x6 = reshape(stm_x_noised, [6,6]);
    %scatter3(Q0(1,:),Q0(2,:), Q0(3,:), 'filled')
    %plot3(phi_Q(:,1), phi_Q(:,2), phi_Q(:,3), 'Color', 'Black', 'LineWidth',1)
    t_list_ = [t_list_; t_];
    F = phi_Q(end,:).' - Q0(:,i+1);
    F_list = [F_list; F];
    df{i,i} = stm_6x6;
    df{i,i+1} = -eye(6);
    for j = 1:N
        if j ~= i && j ~= i+1
            df{i,j} = zeros(6);
        end
    end
end

% Initialize M with zeros (or appropriate initial values)
[row_size, col_size] = size(df{1, 1});  % Assuming all matrices are of the same size
DF = zeros(size(df, 1) * row_size, size(df, 2) * col_size);

% Populate DF with data from df
for i = 1:N-1
    for j = 1:(N)
        % Compute the starting index for each matrix in DF
        start_row = (i - 1) * row_size + 1;
        start_col = (j - 1) * col_size + 1;

        % Assign the matrix from df to the corresponding position in DF
        df{i,j};
        DF(start_row:start_row + row_size - 1, start_col:start_col + col_size - 1) = df{i, j};
    end
end

delta_Q = - DF.' * inv(DF*DF.') * F_list;
delta_Q = reshape(delta_Q, [6,N]);
norm(delta_Q.')
fprintf('Q0 old %f\n', Q0(1,1))
Q0 = Q0 + delta_Q;
fprintf('Q0 new %f\n', Q0(1,1))
%scatter3(Q0(1,:),Q0(2,:), Q0(3,:), 'filled')


end
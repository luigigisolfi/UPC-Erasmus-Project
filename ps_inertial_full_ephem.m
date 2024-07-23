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
FRAME = 'J2000';
OBSERVER = 'EARTH';
BODIES = [{'JUPITER BARYCENTER'}, {'MARS BARYCENTER'},{'MOON'}, {'SATURN BARYCENTER'}, {'MERCURY'}, {'VENUS'}, {'URANUS BARYCENTER'}, {'PLUTO BARYCENTER'}, {'NEPTUNE BARYCENTER'}];
PRIMARIES = [{'SUN'}, {'EARTH'}];
L = get_L(FRAME, PRIMARIES); % Computes the mean distance of two given primaries over a timespan of 50 years 
mu = get_mu(PRIMARIES); %Compute grav. parameter for the system
MODEL = '@etbp'; 
%-------- Compute mean anomaly (global variable, needed for time converison) ---------------%
%n_anomalistic = 2.639394888546285e-06 %computed with the function n_peri_peri, corresponding to the anomalistic month FOR MOON EARTH
n_anomalistic = 1.991059558724508e-07 %FOR EARTH SUN;
n_rtbp = get_n(PRIMARIES,L); %computed with Kepler's 3rd Law
%------------------------------------------------------------------------------%

%-----Set LUNAR ECLIPSE of 21 JAN 2000 as starting point in time (taken from almanacs and verified with a SPICE plot) -----%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
eclipse_date_UTC= '21 JAN 2000 04:05:02';
eclipse_date_et = cspice_str2et(eclipse_date_UTC);
t_list = linspace(eclipse_date_et,eclipse_date_et + 100*pi/n_rtbp, 27*100); %propagation times for get_ephemeris function.
%-------------------------------------------------------------------------------------------------------------%

%--------- Retrieve INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
[inertial_state_primaries, inertial_state_bodies, interpolators] = get_ephemeris(t_list, PRIMARIES, BODIES, FRAME, OBSERVER);
N = 50; %number of nodes

%mu = 0.0122; %EM
mu = 3.040423398444176E-06; %SE
%orbit_file = 'lissa_cut.txt' %Masde orbits
%orbit_file = 'circular_orbit_Earth_Sun.txt'
orbit_file = 'halo_demo_master.txt';

% Retrieve (read) rtbp orbit
[t, x] = read_orbit(orbit_file);
%rtbp times
ti = t(1); % 
tf = t(end);  
% 
% figure
% hold on
% plot3(x(:,1), x(:,2), x(:,3))
% scatter3(mu-1,0,0)
% return

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
rs_rp = rp - rs;
vs_vp = vp -vs;
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
[rtbp_pos_spacecraft_] = go_synodic_pos_only(rs_rp, vs_vp, SEb_pos,inertial_pos_spacecraft.');
% 
% figure
% hold on
% axis equal
% plot3(inertial_pos_spacecraft(1,:),inertial_pos_spacecraft(2,:),inertial_pos_spacecraft(3,:) )
% scatter3(rs(1,:), rs(2,:), rs(3,:))
% return
%plot3(rp(1,:), rp(2,:), rp(3,:))
Q0 = [inertial_pos_spacecraft;inertial_vel_spacecraft];
phi_Q_list = [];
t_inertial_list = [];
%-----------------------------------------------------------------------------------------------------------------------%
for iteration = 1:5
fprintf('iteration %f\n', iteration)
t_list_ = [];
F_list = [];
df = cell(N-1, N);

fprintf('Q0 to recover %f\n', inertial_pos_spacecraft(1,1))
fprintf('Q0 old %f\n', Q0(1,1))
    xiv=eye(6,6);

for i = 1:N-1
    xiv=eye(6,6);
    [t_, phi_Q_tot] = ode113(@new_full_force_var_vectorized, [t_sampled_inertial(i), t_sampled_inertial(i+1)], [Q0(:,i);xiv(:)]);
    phi_Q = phi_Q_tot(:,1:6);
    if iteration == 5
        t_inertial_list = [t_inertial_list; t_];
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

delta_Q = - DF.' * pinv(DF*DF.')*F_list;
delta_Q = reshape(delta_Q, [6,N]);
norm(delta_Q.')
fprintf('Q0 old %f\n', Q0(1,1))
Q0 = Q0 + delta_Q;
fprintf('Q0 new %f\n', Q0(1,1))
end
%------------------------------------------------------------------------------------------------------------------------%

% ------------------------------------------------------------------------------------------------------------------------%
% Prepare data matrix for saving inertial Lissajous data in a text file. 
data = [t_inertial_list(:), phi_Q_list(:,1),phi_Q_list(:,2),phi_Q_list(:,3),phi_Q_list(:,4),phi_Q_list(:,5),phi_Q_list(:,6)]; % Combine time with x, y, z coordinates

% Save data to Lissajous inertial text file
filename_end = 'lissajous_inertial.txt';
fid = fopen(filename_end, 'w');
fprintf(fid, 'Time (s)\tX (km)\tY (km)\tZ (km)\tvx (km)\tvy (km)\tvz (km)\n');
fprintf(fid, '%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n', data.');
fclose(fid);
disp(['Orbit coordinates saved to ' filename_end]);
%------------------------------------------------------------------------------------------------------------------------%

[inertial_state_primaries_new, inertial_state_bodies_new, interpolators_new] = get_ephemeris(t_list, PRIMARIES, BODIES, FRAME, OBSERVER);
t_list(1)
t_list(end)
t_list_(1)
t_list_(end)

% Initialize arrays to store interpolated position and velocity
len_t_list_ = length(t_list_(:));
rp_new = zeros(3, len_t_list_);
vp_new = zeros(3, len_t_list_);
rs_new = zeros(3, len_t_list_);
vs_new = zeros(3, len_t_list_);
ap_new = zeros(3, len_t_list_);
as_new = zeros(3, len_t_list_);

for dim = 1:9
    interpolated_p_new = ppval(interpolators_new.(primary_str).spline{dim}, t_list_(:));
    interpolated_s_new = ppval(interpolators_new.(secondary_str).spline{dim}, t_list_(:));
    if dim <= 3
        rp_new(dim, :) = interpolated_p_new;
        rs_new(dim,:) = interpolated_s_new;
    elseif dim>=4 && dim<=6 
        vp_new(dim-3, :) = interpolated_p_new;
        vs_new(dim-3, :) = interpolated_s_new;      
    elseif dim>=7 && dim<=9
        ap_new(dim-6, :) = interpolated_p_new;
        as_new(dim-6, :) = interpolated_s_new; 
    end
    
end


%retrieving primaries relative acceleration at each rtbp time
rs_rp_new = rp_new-rs_new;
vs_vp_new = vp_new-vs_new;
as_ap_new = ap_new - as_new;

%barycenter
SEb_pos_new = rp_new - mu*rs_rp_new;
SEb_vel_new = vp_new - mu*vs_vp_new;


%figure
%hold on
%axis equal
%plot3(phi_Q_list(:,1), phi_Q_list(:,2), phi_Q_list(:,3))
%plot3(rp_new(1,:), rp_new(2,:), rp_new(3,:))
%scatter3(rs_new(1,:), rs_new(2,:), rs_new(3,:))


[rtbp_pos_spacecraft, rtbp_vel_spacecraft] = go_synodic_pos_vel_only(rs_rp_new,vs_vp_new, as_ap_new, SEb_pos_new, SEb_vel_new, phi_Q_list, n_rtbp);
% Save data to Lissajous inertial text file
filename_end = 'test_circular_for_comparison.txt';
size(rtbp_pos_spacecraft)
size(rtbp_vel_spacecraft)
size(t_inertial_list(:))
data_new = [t_inertial_list(:)*n_rtbp, rtbp_pos_spacecraft(1,:).',rtbp_pos_spacecraft(2,:).',rtbp_pos_spacecraft(3,:).',rtbp_vel_spacecraft(1,:).',rtbp_vel_spacecraft(2,:).',rtbp_vel_spacecraft(3,:).']; % Combine time with x, y, z coordinates
fid = fopen(filename_end, 'w');
fprintf(fid, 'Time (s)\tX (km)\tY (km)\tZ (km)\tvx (km)\tvy (km)\tvz (km)\n');
fprintf(fid, '%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n', data_new.');
fclose(fid);
disp(['Orbit coordinates saved to ' filename_end]);
%return
%figure
% hold on
%plot3(rtbp_pos_spacecraft(1,:), rtbp_pos_spacecraft(2,:), rtbp_pos_spacecraft(3,:))%
%scatter3(mu-1,0,0)

% Create the 3D plot
figure;
hold on
axis equal
h = plot3(rtbp_pos_spacecraft(1,:), rtbp_pos_spacecraft(2,:), rtbp_pos_spacecraft(3,:), 'LineWidth', 1.2, 'DisplayName','F. Ephem Jup + Mar + Moon');
h2 = plot3(x(:,1), x(:,2), x(:,3), 'LineWidth', 1, 'LineStyle','--', 'DisplayName', 'RTBP');
scatter3(mu-1, 0, 0, 'filled', 'MarkerFaceColor', 'r');  % Plot the fixed point
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Orbit');
legend('show')

% Set up the view angles
az = linspace(0, 90, 90);  % Azimuth angles
el = linspace(0, 50, 90);   % Constant elevation angle

% Capture frames while rotating the view
filename = '3d_circular_spacecraft.gif';
for i = 1:length(az)
    view(az(i), el(i));
    drawnow;
    
    % Capture the current frame
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Write to the GIF file
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.01);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.01);
    end
end

disp(['GIF saved as ', filename]);

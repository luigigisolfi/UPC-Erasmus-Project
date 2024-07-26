clear all; close all; clc;

global avg_e
global inertial_state_primaries
global interpolators
global n_rtbp
global n_anomalistic
global PRIMARIES
global mu
global BODIES
global L
global eclipse_date_et
global nu


% --------------- LOAD KERNELS -------------
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %

%  --------------- SET THE MODEL -------------
FRAME = 'ECLIPJ2000';
OBSERVER = 'EARTH-MOON BARYCENTER';
BODIES = [{'SUN'},{'JUPITER BARYCENTER'}];
PRIMARIES = [{'EARTH'}, {'MOON'}];
L = get_L(FRAME, PRIMARIES); % Computes the mean distance of two given primaries over a timespan of 50 years 
%L = 384601.25606767; %value form Gomez et al
mu = get_mu(PRIMARIES) %Compute grav. parameter for the system
MODEL = '@etbp'; 
%-------- Compute mean anomaly (global variable, needed for time converison) ---------------%
n_anomalistic = 2.639394888546285e-06 %computed with the function n_peri_peri, corresponding to the anomalistic month
n_rtbp = get_n(PRIMARIES,L) %computed with Kepler's 3rd Law
%------------------------------------------------------------------------------%

%-----Set LUNAR ECLIPSE of 21 JAN 2000 as starting point in time (taken from almanacs and verified with a SPICE plot) -----%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
eclipse_date_UTC= '21 JAN 2000 04:05:02';
eclipse_date_et = cspice_str2et(eclipse_date_UTC);
t_list = linspace(eclipse_date_et,eclipse_date_et + 100*pi/n_rtbp, 27*100);
cspice_kclear()
%-------------------------------------------------------------------------------------------------------------%

%--------- Retrieve INERTIAL state vectors of the two PRIMARIES in the ECLIPJ2000 frame (from SPICE) ---------%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
[inertial_state_primaries, inertial_state_bodies, interpolators] = get_ephemeris(t_list, PRIMARIES, BODIES, FRAME, OBSERVER);
% VERIFY ECLIPSE
% figure
% hold on
% scatter3(inertial_state_primaries.EARTH.position(1),inertial_state_primaries.EARTH.position(2), inertial_state_primaries.EARTH.position(3), 'DisplayName','EARTH')
% scatter3(inertial_state_primaries.MOON.position(1),inertial_state_primaries.MOON.position(2),inertial_state_primaries.MOON.position(3))
% scatter3(inertial_state_bodies.SUN.position(1),inertial_state_bodies.SUN.position(2),inertial_state_bodies.SUN.position(3))
% %quiver3(inertial_state_primaries.EARTH.position(1), inertial_state_primaries.EARTH.position(2), inertial_state_primaries.EARTH.position(3), inertial_state_bodies.SUN.position(1),inertial_state_bodies.SUN.position(2), inertial_state_bodies.SUN.position(3), 'r', 'LineWidth', 2, 'DisplayName','JUPITER BARYCENTER');
% cspice_kclear()
% legend('show')
% return
%-------------------------------------------------------------------------------------------------------------%

% ------ Critical Part: Defining the ephemeris times (inertial, common for all models) ------ %
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
ti_utc = cspice_et2utc(t_list(1), 'C', 0);
tf_utc = cspice_et2utc(t_list(end), 'C', 0);
cspice_kclear()

fprintf('-----------------------------------------------------------\n')
fprintf('MAIN\nUTC Start Simulation Epoch: %s\n', ti_utc)
fprintf('UTC End Simulation Epoch: %s\n', tf_utc)
fprintf('-------------------------------------------------------------------------------------------------------------------\n')
fprintf('MAIN\nCalling get_ephemeris function...\n\n')
% -----------------------------------------------------------------%


% ------------------------ RTBP variables ---------------------------%
[mu, body_1, body_2] = get_mu(PRIMARIES);
GM_1 = get_GM_body(body_1);
GM_2 = get_GM_body(body_2);
%---------------------------------------------------------------------%

% ----- Define the evaluation grid --------- %
%x1 = linspace(mu-1-0.01, mu-1+0.01, 10); % Range for x(1)
x1 = linspace((mu-1)-sqrt((mu)/3), (mu-1)+sqrt((mu)/3), 4); % Range for x(1)
%x1 = linspace((mu-1)-sqrt((mu)/3), (mu)+sqrt((mu)/3), 20); % Range for x(1)
x2 = linspace(-0.07,0.07, 4); % Range for x(2)
[X1, X2] = meshgrid(x1, x2);
%---------------------------------------------------------------------%

% ---- Initialize and compute acceleration values for contour plot ----- %
Z = zeros(size(X1));
Z_rtbp = zeros(size(X1));
Z_etbp = zeros(size(X1));

avg_e = get_avg_e(inertial_state_primaries, PRIMARIES);
nu = 0.5412;
model_indices = zeros(size(X1));
indexes = [1,2,3,4,5]
for p = 1:length(indexes)
 if strcmp(MODEL, '@etbp')
    p
    vector_field = zeros(3, length(X1),length(X1));
    index = indexes(p);
    t_HFEM = t_list(index);
    META = 'kernels_to_load.tm'; %initialize required kernels
    cspice_furnsh(META); %furnish kernels
    t_peri_peri = (t_HFEM - eclipse_date_et)*n_anomalistic;
    t_rtbp = (t_HFEM - eclipse_date_et)*n_rtbp;
    t_utc = cspice_et2utc(t_HFEM, 'C',3);
    t_etbp_utc = cspice_et2utc(t_peri_peri/n_anomalistic + eclipse_date_et, 'C',3);
    t_rtbp_utc = cspice_et2utc(t_rtbp/n_rtbp + eclipse_date_et, 'C',3);
    cspice_kclear()

% Define a state vector for each grid point (considering only 2D for simplicity)
    for i = 1:length(x1)
        i
        for j = 1:length(x2)
            x = [X1(i,j); X2(i,j); 0; 0.01; 0.01; 0]; % Example state vector, modify as needed
            % Compute the vector fields
            f1_val = HFEM(t_peri_peri, x) %translation of time not needed: we want HFEM to start from the eclipse time, which is t_HFEM = t_list(1)*n
            f2_val_full = HFEM_etbp_b_coefficients(t_peri_peri, x); %translation of time needed: we want etbp to start from zero at t_HFEM. Note, however, the avg_nu_dot conversion
            f2_val = f2_val_full(1:6)
            f3_val = HFEM_rtbp(t_rtbp, x)
            f1_val-f2_val
            f4_val_full = HFEM_residual_etbp((t_HFEM)*n_anomalistic,x);
            f4_val = f4_val_full(1:6)
            nu_dot_res = f4_val_full(7);
            nu_dot = f2_val_full(7);
            Z_rtbp(i,j) =  abs(norm(f1_val(4:6))- norm(f3_val(4:6))) %Difference of magnitudes of the acceleration vectors HFEM - RTBP
            Z_etbp(i,j) =  abs(norm(f1_val(4:6))- norm(f2_val(4:6)))%Difference of magnitudes of the acceleration vectors HFEM - ETBP
            vector_field(:,i,j) = f1_val(4:6);
        end
    end
 else
    p
    t_HFEM = t_list(index);
% Define a state vector for each grid point (considering only 2D for simplicity)
    for i = 1:length(x1)
        for j = 1:length(x2)
            x = [X1(i,j); X2(i,j); 0; 0; 0; 0]; % Example state vector, modify as needed
            % Compute the vector fields
            f1_val = HFEM(t_HFEM*n, x);
            f2_val = rtbp(t_HFEM*n, x);
            Z(i,j) = abs(norm(f1_val(4:6)) - norm(f2_val(4:6))); % Difference in acceleration (or any other component)
        end
    end     
 end

nu = nu_dot*t_peri_peri + nu %update the true anomaly assuming linearity (works quite well, same was done in Beom Park et al)
nu_res = nu_dot_res*t_peri_peri + nu


 % ------- Initialize and compute index matrix, which will tell us if Z_rtbp >= Z_etbp or viceversa ----- %
index_matrix = zeros(size(Z_rtbp));
for i = 1:size(Z_rtbp, 1)
    for j = 1:size(Z_rtbp, 2)
        if Z_rtbp(i, j) < Z_etbp(i, j)
            Z(i, j) = Z_rtbp(i, j);
            index_matrix(i, j) = 1;  % Indicates the value is from Z_rtbp
        elseif Z_etbp(i, j) < Z_rtbp(i, j)
            Z(i, j) = Z_etbp(i, j);
            index_matrix(i, j) = 2;  % Indicates the value is from Z_etbp
        else 
            Z(i,j) = Z_etbp(i,j);
            index_matrix(i,j) = 3;
        end
    if isnan(Z(i,j))
        fprintf('found a NaN at coordinates: %d, %d', X1(i,j),X2(i,j))
        Z(i,j) = 0;
    end

    end
end
% ------------------------------------------------------------------------------------------------%

% ----- Retrieve Additional Bodies Positions ---------------------------------------%
%inertial_pos_j =zeros(3,1);
%inertial_pos_j(1) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{1}, t_HFEM);
%inertial_pos_j(2) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{2}, t_HFEM);
%inertial_pos_j(3) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{3}, t_HFEM);
inertial_pos_m =zeros(3,1);
inertial_pos_m(1) =  ppval(interpolators.('SUN').spline{1}, t_HFEM);
inertial_pos_m(2) =  ppval(interpolators.('SUN').spline{2}, t_HFEM);
inertial_pos_m(3) =  ppval(interpolators.('SUN').spline{3}, t_HFEM);
%----------------------------------------------------------------------------------------%

% ----------------Evaluate the splines for the primaries
% ------------------------%
primary_name = PRIMARIES{1};
secondary_name = PRIMARIES{2};
primary_str = regexprep(primary_name, [{'\s+'}, {'-'}], '_');
secondary_str = regexprep(secondary_name, [{'\s+'}, {'-'}], '_');

% Initialize arrays to store interpolated position and velocity
rp = zeros(3, 1);
vp = zeros(3, 1);

rs = zeros(3, 1);
vs = zeros(3, 1);


for dim = 1:12
    interpolated_p = ppval(interpolators.(primary_str).spline{dim}, t_HFEM);
    interpolated_s = ppval(interpolators.(secondary_str).spline{dim}, t_HFEM);
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
%----------------------------------------------------------------%

%retrieving primaries relative acceleration at time t
rs_rp = rp-rs;
vs_vp = vp-vs;
as_ap = ap-as;
oas_oap = oap-oas;

SEb_pos = rp - mu*rs_rp;
SEb_vel = vp - mu*vs_vp;
SEb_acc = ap - mu*as_ap;

b = SEb_pos;
b_dot = SEb_vel;
b_ddot = SEb_acc;
k = norm(rs_rp);

C = construct_C(rs_rp, vs_vp);
%rtbp_pos_j =  C\(inertial_pos_j-b)/k;
rtbp_pos_m =  C\(inertial_pos_m-b)/k;
% Plot the difference as a contour plot

figure;
hold on
axis equal
contourf(X1, X2, Z);
colorbar;

% Find the indices where index_matrix is 1
[row, col] = find(index_matrix == 1);

% Extract the corresponding X1, X2, and Z values
crossX1 = X1(sub2ind(size(X1), row, col));
crossX2 = X2(sub2ind(size(X2), row, col));
crossZ = Z(sub2ind(size(Z), row, col));
% Overlay the crosses on the contour plot
scatter3(crossX1, crossX2, crossZ, 'xk','DisplayName', 'RTBP');

[row_, col_] = find(index_matrix == 2);

% Extract the corresponding X1, X2, and Z values
crossX1_ = X1(sub2ind(size(X1), row_, col_));
crossX2_ = X2(sub2ind(size(X2), row_, col_));
crossZ_ = Z(sub2ind(size(Z), row_, col_));
% Overlay the crosses on the contour plot
scatter3(crossX1_, crossX2_, crossZ_, 'o', 'DisplayName', 'ETBP');

[row__, col__] = find(index_matrix == 3);
% Extract the corresponding X1, X2, and Z values
crossX1__ = X1(sub2ind(size(X1), row__, col__));
crossX2__ = X2(sub2ind(size(X2), row__, col__));
crossZ__ = Z(sub2ind(size(Z), row__, col__));
% Overlay the crosses on the contour plot
scatter3(crossX1__, crossX2__, crossZ__, 's', 'DisplayName', 'SAME');

%quiver3(mu-1, 0, 0, rtbp_pos_j(1)/(norm(rtbp_pos_j)*100), rtbp_pos_j(2)/(norm(rtbp_pos_j)*100), 0, 'r', 'LineWidth', 2, 'DisplayName','JUPITER BARYCENTER');
quiver3(mu-1, 0, 0, rtbp_pos_m(1)/(norm(rtbp_pos_m)*100), rtbp_pos_m(2)/(norm(rtbp_pos_m)*100), 0, 'y', 'LineWidth', 2, 'DisplayName','SUN DIRECTION');

for g = 1:(length(X1))
    for d = 1:length(X1)
        quiver3(X1(g,d), X2(g,d),0, vector_field(1,g,d)*0.01/(norm(vector_field(:,g,d))),vector_field(2,g,d)*0.01/norm(vector_field(:,g,d)),0, 'HandleVisibility','off', 'LineWidth',1, 'Color','red')
    end
end

scatter(mu-1,0, 'filled', 'DisplayName','MOON')
th = 0:pi/50:2*pi;
x_hill = sqrt((mu)/3) * cos(th) + mu-1;
y_hill = sqrt((mu)/3) * sin(th) + 0;
plot(x_hill,y_hill, 'LineWidth',2)

if strcmp(MODEL, '@etbp')
    title_str = sprintf('a - difference between HFEM and ETBP (t = %s)', t_utc(1:11));
else
    title_str = sprintf('a - difference between HFEM and RTBP (t = %s)', t_utc(1:11));
end

title(title_str)
xlabel('x (adimensional)');
ylabel('y (adimensional)');
legend('show')
% Save the figure as a PNG
fileout = sprintf('%s.png', t_utc(1:11)); % 'plot_1.png', 'plot_2.png', ...
saveas(gcf, fileout); % Save as PNG

end
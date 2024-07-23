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
FRAME = 'J2000';
OBSERVER = 'SUN';
BODIES = [{'JUPITER BARYCENTER'}, {'MARS BARYCENTER'},{'MOON'}, {'SATURN BARYCENTER'}, {'MERCURY'}, {'VENUS'}, {'URANUS BARYCENTER'}, {'PLUTO BARYCENTER'}, {'NEPTUNE BARYCENTER'}];
PRIMARIES = [{'SUN'}, {'EARTH'}];
L = get_L(FRAME, PRIMARIES); % Computes the mean distance of two given primaries over a timespan of 50 years 
%L = 384601.25606767; %value form Gomez et al
mu = get_mu(PRIMARIES) %Compute grav. parameter for the system
MODEL = '@etbp'; 
%-------- Compute mean anomaly (global variable, needed for time converison) ---------------%
%n_anomalistic = 2.639394888546285e-06 %computed with the function n_peri_peri, corresponding to the anomalistic month
n_anomalistic = 1.991059558724508e-07 %FOR EARTH SUN;
n_rtbp = get_n(PRIMARIES,L); %computed with Kepler's 3rd Law
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
% x1 = linspace((mu-1)-sqrt((mu)/3), (mu-1)+sqrt((mu)/3), 4); % Range for x(1)
% %x1 = linspace((mu-1)-sqrt((mu)/3), (mu)+sqrt((mu)/3), 20); % Range for x(1)
% x2 = linspace(-0.07,0.07, 4); % Range for x(2)
filename = 'test_lissajous_for_comparison.txt';
%filename = 'circular_orbit_Earth_Sun_Full_Ephem.txt'
[t_orb,x_orb] = read_orbit(filename);
%[X, Y] = meshgrid(x1, x2);
%---------------------------------------------------------------------%

% ---- Initialize and compute acceleration values for contour plot ----- %
Z = zeros(size(x_orb(:,1)));
Z_rtbp = zeros(size(x_orb(:,1)));
Z_etbp = zeros(size(x_orb(:,1)));
Z_hill3b = zeros(size(x_orb(:,1)));

avg_e = get_avg_e(inertial_state_primaries, PRIMARIES);
%nu = 0.5412; %for moon earth at eclipse_date_et,  computed with spkezr + rv2cel
nu = 0.3218; %for sun earth at eclipse_date_et, computed with spkezr + rv2cel
vector_field = zeros(3, length(x_orb(:,1)));
% Define a state vector for each grid point (considering only 2D for simplicity)
for i = 1:1:length(x_orb(:,1))
    length(x_orb(:,1))
    i
        t_HFEM = t_orb(i)/n_rtbp;
        t_peri_peri = (t_HFEM - eclipse_date_et)*n_anomalistic;
        t_rtbp = (t_HFEM - eclipse_date_et)*n_rtbp;
        x = [x_orb(i,1), x_orb(i,2), x_orb(i,3), x_orb(i,4), x_orb(i,5), x_orb(i,6)].'; % Example state vector, modify as needed
        % Compute the vector fields
        f1_val = HFEM(t_peri_peri, x); %translation of time not needed: we want HFEM to start from the eclipse time, which is t_HFEM = t_list(1)*n
        
        f2_val_full = HFEM_etbp_b_coefficients(t_peri_peri, x); %translation of time needed: we want etbp to start from zero at t_HFEM. Note, however, the avg_nu_dot conversion
        f2_val = f2_val_full(1:6);

        f3_val = HFEM_rtbp(t_rtbp, x);

        f4_val = hill3b(t_rtbp, x);
        
        %f4_val_full = HFEM_residual_etbp((t_HFEM)*n_anomalistic,x);
        %f4_val = f4_val_full(1:6);

        %nu_dot_res = f4_val_full(7); %HFEM residual etbp gives nu as output, besides the vector state (6 + 1 elements)
        nu_dot = f2_val_full(7); %HFEM etbp b coefficientd also gives the nu value. Here, we compute them to compare them

        Z_rtbp(i) =  abs(norm(f1_val(4:6) - f3_val(4:6))); %Difference of magnitudes of the acceleration vectors HFEM - RTBP
        Z_etbp(i) =  abs(norm(f1_val(4:6) - f2_val(4:6)));%Difference of magnitudes of the acceleration vectors HFEM - ETBP
        Z_hill3b(i) = abs(norm(f1_val(4:6) - f4_val(4:6)));%Difference of magnitudes of the acceleration vectors HFEM - HILL3B

        %Z_etbp(i) =  (norm(f4_val(4:6)));%Difference of magnitudes of the acceleration vectors HFEM - ETBP
        
        vector_field(:,i) = f1_val(4:6);
        nu = nu_dot*t_peri_peri + nu; %update the true anomaly assuming linearity (works quite well, same was done in Beom Park et al)

        %nu_res = nu_dot_res*t_peri_peri + nu;
end

 % ------- Initialize and compute index matrix, which will tell us if Z_rtbp >= Z_etbp or viceversa ----- %
index_matrix = zeros(size(Z_rtbp));
for i = 1:size(Z_rtbp, 1)
        if Z_rtbp(i) < Z_etbp(i) && Z_rtbp(i) < Z_hill3b(i)
            Z(i) = Z_rtbp(i);
            index_matrix(i) = 1;  % Indicates the value is from Z_rtbp
        elseif Z_etbp(i) < Z_rtbp(i) && Z_etbp(i) < Z_hill3b(i)
            Z(i) = Z_etbp(i);
            index_matrix(i) = 2;  % Indicates the value is from Z_etbp
        elseif Z_hill3b(i) < Z_rtbp(i) && Z_hill3b(i) < Z_etbp(i)
            Z(i) = Z_hill3b(i);
            index_matrix(i) = 3;
        else 
            Z(i) = Z_etbp(i);
            index_matrix(i) = 4;
        end
    if isnan(Z(i))
        fprintf('found a NaN at coordinates: %d, %d, %d.\n', x_orb(i,1),x_orb(i,2),x_orb(i,3))
        Z(i) = 0;
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
[X,Y] = meshgrid(x_orb(:,1), x_orb(:,2));
% create corresponding Z values, assume z = 0 for locations with no z data
Z_grid = zeros(length(x_orb(:,1)),length(x_orb(:,2))) ;
for i = 1:length(x_orb(:,1))
    for j = 1:length(x_orb(:,2))
        if i==j % z data exist for only for x(n) y(n) location, n = 1,2,3...
        Z_grid(i,j) = Z(i);
        end
    end
end

% contourf(X,Y,Z_grid)
% title('Contour plot of interpolated data');
% xlabel('x');
% ylabel('y');
% colorbar;

% Normalize L for color scaling
Z_min = min(Z);
Z_max = max(Z);
Z_norm = (Z - Z_min)/ (Z_max - Z_min); % Normalized to [0, 1]

% Create colormap from blue to red
cmap = jet; % Jet colormap from blue to red
% Get the number of colors in the colormap
nColors = size(cmap, 1);
% Map the normalized Z values to colormap indices
colorIdx = round(Z_norm * (nColors - 1)) + 1;
colors = cmap(colorIdx, :);
opacity = 0.6;
scatter3(x_orb(:,1), x_orb(:,2), x_orb(:,3), 36, colors, 'filled', 'MarkerFaceAlpha', opacity, 'MarkerEdgeAlpha', opacity);
% Add colorbar for reference
colormap(cmap);
c = colorbar;
c.Ticks = linspace(0, 1, 11); % 11 ticks from 0 to 1
c.TickLabels = num2str(linspace(Z_min, Z_max, 11)', '%.2f'); % Corresponding values of L
c.Label.String = 'Value of Z';

% Find the indices where index_matrix is 1
[row] = find(index_matrix == 1);

% Extract the corresponding X, Y, and Z values
crossX1 = x_orb(row,1);
crossX2 = x_orb(row,2);
crossX3 = x_orb(row,3);
% Overlay the crosses on the contour plot
scatter3(crossX1, crossX2,crossX3,80,  'xk','DisplayName', 'RTBP');

[row_] = find(index_matrix == 2);

% Extract the corresponding X, Y, and Z values
crossX1_ = x_orb(row_,1);
crossX2_ = x_orb(row_,2);
crossX3_ = x_orb(row_,3);
% Overlay the crosses on the contour plot
scatter3(crossX1_, crossX2_,crossX3_,80,'o','DisplayName', 'ETBP');

[row__] = find(index_matrix == 3);
% Extract the corresponding X, Y, and Z values
crossX1__ = x_orb(row__,1);
crossX2__ = x_orb(row__,2);
crossX3__ = x_orb(row__,3);
% Overlay the crosses on the contour plot
scatter3(crossX1__, crossX2__, crossX3__,80, 's', 'DisplayName', 'hill3b');

[row___] = find(index_matrix == 4);
% Extract the corresponding X, Y, and Z values
crossX1___ = x_orb(row___,1);
crossX2___ = x_orb(row___,2);
crossX3___ = x_orb(row___,3);
% Overlay the crosses on the contour plot
%scatter3(crossX1___, crossX2___,crossX3___,'s',45, 'DisplayName', 'SAME');
legend('show')

%compute best model
alpha = 0.5;
avg_rtbp = mean(Z(row))
avg_etbp = mean(Z(row_))
avg_hill3b = mean(Z(row__))

score_rtbp = length(row)/length(Z)*(1-alpha) + min(min(avg_rtbp, avg_etbp), avg_hill3b)*alpha/avg_rtbp

score_etbp = length(row_)/length(Z)*(1-alpha) + min(min(avg_rtbp, avg_etbp), avg_hill3b)*alpha/avg_etbp

score_hill3b = length(row__)/length(Z)*(1-alpha) + min(min(avg_rtbp, avg_etbp), avg_hill3b)*alpha/avg_hill3b

return

%quiver3(mu-1, 0, 0, rtbp_pos_j(1)/(norm(rtbp_pos_j)*100), rtbp_pos_j(2)/(norm(rtbp_pos_j)*100), 0, 'r', 'LineWidth', 2, 'DisplayName','JUPITER BARYCENTER');
quiver3(mu-1, 0, 0, rtbp_pos_m(1)/(norm(rtbp_pos_m)*100), rtbp_pos_m(2)/(norm(rtbp_pos_m)*100), 0, 'y', 'LineWidth', 2, 'DisplayName','SUN DIRECTION');

for g = 1:(length(X))
    for d = 1:length(X)
        quiver3(X(g,d), Y(g,d),0, vector_field(1,g)*0.0001/(norm(vector_field(:,g))),vector_field(2,g)*0.0001/norm(vector_field(:,g)),0, 'HandleVisibility','off', 'LineWidth',1, 'Color','red')
    end
end
return


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

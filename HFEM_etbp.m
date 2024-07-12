function xdot = HFEM_etbp(t,x)

global PRIMARIES
global BODIES
global interpolators
global nu_dot
global mu
global n
global nu
global r
global L
%here, t = rtbp_times. But interpolators are given in inertial physical
%coordinates. therefore, if we want to retrieve rs_rp, we need to write
%interpolators(inertial_t) where inertial_t = t*n = rtbp_times*n


META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
inertial_t = t/n;
utc_t = cspice_et2utc(inertial_t, 'C', 3);
% Evaluate the splines for the primaries
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
    interpolated_p = ppval(interpolators.(primary_str).spline{dim}, inertial_t);
    interpolated_s = ppval(interpolators.(secondary_str).spline{dim}, inertial_t);
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

%retrieving relative acceleration at time t
rs_rp = rs-rp;
vs_vp = vs-vp;
as_ap = as-ap;
oas_oap = oas-oap;

SEb_pos = rp + mu*rs_rp;
SEb_vel = vp + mu*vs_vp;
SEb_acc = ap + mu*as_ap;

b = SEb_pos;
b_dot = SEb_vel;
b_ddot = SEb_acc;
k = norm(rs_rp);
k_dot = dot(rs_rp, vs_vp)/k;
k_ddot = (dot(vs_vp, vs_vp) + dot(rs_rp, as_ap) - k_dot^2)/k;
h = norm(cross(rs_rp, vs_vp));
hp = dot(cross(rs_rp, as_ap), cross(rs_rp, vs_vp))/h;

%
b1 = -b_ddot(1)/(nu_dot^2*k);
b2 = -b_ddot(2)/(nu_dot^2*k);
b3 = -b_ddot(3)/(nu_dot^2*k);
b4 = -k_dot/(2*nu_dot*k);
b5 = 2*h/(nu_dot*k^2);
b6 = 2*k*as_ap(3)/(nu_dot*h);
b7 = -k_ddot/(nu_dot^2*k) + h^2/(nu_dot^2*k^4);
b8 = - as_ap(3)/(nu_dot^2*k);
b9 = hp/(nu_dot^2*k^2);
b10 = -k_ddot/(nu_dot^2*k) + h^2/(nu_dot^2*k^4) + k^2*(as_ap(3))^2/(nu_dot^2*h^2);
b11 = (3*h*k_dot - 2*k*hp)*(as_ap(3))/(nu_dot^2*h^2) + k*oas_oap(3)/(nu_dot^2*h);
b12 = - k_ddot/(nu_dot^2*k) + k^2*(as_ap(3))^2/(nu_dot^2*h^2);
b13 = 1;

first = [b1;b2;b3];
second = [b4,b5,0;-b5,b4,b6;0,-b6,b4];
third = [b7,b9,b8;-b9,b10,b11;b8,-b11,b12];

all_bodies = [BODIES];
n_bodies = length(all_bodies);
acc_list = zeros(3, n_bodies);
xdot = zeros(6,1);

rtbp_pos = x(1:3);
rtbp_vf = rtbp(t,x);
rtbp_vel = x(4:6);

C = construct_C(rs_rp, vs_vp);
C_dot = construct_C_dot(C, rs_rp, vs_vp, as_ap);
C_ddot = construct_C_ddot(C, C_dot, rs_rp, vs_vp, as_ap, oas_oap);


% This was just a check. Not needed
% rtbp_acc = rtbp_vf(4:6);
%inertial_rtbp_acc = b_ddot + (k_ddot*C + 2*k_dot*C_dot+ k*C_ddot)*rtbp_pos + (2*k_dot*C + 2*k*C_dot)*rtbp_vel*n + k*C*n^2*rtbp_acc;
%new_rtbp_acc = (C*k*n^2)\(inertial_rtbp_acc - b_ddot - (k_ddot*C + 2*k_dot*C_dot+ k*C_ddot)*rtbp_pos - (2*k_dot*C + 2*k*C_dot)*rtbp_vel*n);

inertial_pos_spacecraft = k*C*rtbp_pos+ b ;
inertial_vel_spacecraft = b_dot + k_dot*C*rtbp_pos + k*(C_dot*rtbp_pos+ nu_dot*C*rtbp_vel);
x_inertial = [inertial_pos_spacecraft, inertial_vel_spacecraft];

GM_1 = get_GM_body(PRIMARIES{1});
GM_2 = get_GM_body(PRIMARIES{2});
% % 
% % 
% for i = 1:n_bodies % PARALLELIZED COMPUTATION ON ALL BODIES
%     BODY = all_bodies{i};
%     inertial_pos_body =zeros(3,1);
%     BODY_str = regexprep(BODY, [{'\s+'}, {'-'}], '_');
%     mu_body = get_GM_body(BODY)/(GM_1+GM_2);
%     inertial_pos_body(1) =  ppval(interpolators.(BODY_str).spline{1}, inertial_t);
%     inertial_pos_body(2) =  ppval(interpolators.(BODY_str).spline{2}, inertial_t);
%     inertial_pos_body(3) =  ppval(interpolators.(BODY_str).spline{3}, inertial_t);
% 
%     rtbp_pos_body =  C\(inertial_pos_body-b)/k; 
%     x_sb = x(1:3)-rtbp_pos_body;
%     rho_sb = (x(1)-rtbp_pos_body(1))^2+(x(2)-rtbp_pos_body(2))^2+(x(3)-rtbp_pos_body(3))^2;
%     rho_sb3 = rho_sb*sqrt(rho_sb);
%     synodic_acc = -mu_body*x_sb/rho_sb3;
%      %now compute the gravitational_acceleration_synodic(all bodies)
%      %the function gravitational_acceleration gives acc in the inertial frame, so we need to divide the GMj contribution by (GM_1+GM_2) to get mu_j
%     acc_list(:, i) =synodic_acc;
% end

Delta_Omega = sum(acc_list, 2);

inertial_pos_body_1_x =  ppval(interpolators.(PRIMARIES{1}).spline{1}, inertial_t);
inertial_pos_body_1_y =  ppval(interpolators.(PRIMARIES{1}).spline{2}, inertial_t);
inertial_pos_body_1_z =  ppval(interpolators.(PRIMARIES{1}).spline{3}, inertial_t);
inertial_pos_body_1 = [inertial_pos_body_1_x; inertial_pos_body_1_y; inertial_pos_body_1_z];
inertial_pos_body_2_x =  ppval(interpolators.(PRIMARIES{2}).spline{1}, inertial_t);
inertial_pos_body_2_y =  ppval(interpolators.(PRIMARIES{2}).spline{2}, inertial_t);
inertial_pos_body_2_z =  ppval(interpolators.(PRIMARIES{2}).spline{3}, inertial_t);
inertial_pos_body_2 = [inertial_pos_body_2_x; inertial_pos_body_2_y; inertial_pos_body_2_z];

rtbp_pos_body_1 = -(C\(inertial_pos_body_1-b)/k);%should I normalize the lengths so as to  have homogenous comparison in the adimensional coordinates?
rtbp_pos_body_2 = -(C\(inertial_pos_body_2-b)/k); %should I normalize the lengths so as to  have homogenous comparison in the adimensional coordinates?

% fprintf('HFEM_etbp POSITIONS PRIMARIES NORMALIZED TO RTBP LENGTH UNITS')
% rtbp_pos_body_1
% rtbp_pos_body_2

x_s1 = x(1:3)-rtbp_pos_body_1;
x_s2 = x(1:3)-rtbp_pos_body_2;

rho_s1 = (x(1)-rtbp_pos_body_1(1))^2+(x(2))^2+(x(3))^2;
rho_s2 = (x(1)-rtbp_pos_body_2(1))^2+(x(2))^2+(x(3))^2;
rho_s13 = rho_s1*sqrt(rho_s1);
rho_s23 = rho_s2*sqrt(rho_s2);

synodic_acc_primaries = -(1-mu)*x_s1/rho_s13 - (mu)*x_s2/rho_s23; %primaries contribution to synodic acceleration

if isempty(BODIES) == 0
synodic_acc_bodies  = 0;
for i = 1:n_bodies % PARALLELIZED COMPUTATION ON ALL BODIES
    BODY = all_bodies{i};
    inertial_pos_body =zeros(3,1);
    BODY_str = regexprep(BODY, [{'\s+'}, {'-'}], '_');
    mu_body = get_GM_body(BODY)/(GM_1+GM_2);
    inertial_pos_body(1) =  ppval(interpolators.(BODY_str).spline{1}, inertial_t);
    inertial_pos_body(2) =  ppval(interpolators.(BODY_str).spline{2}, inertial_t);
    inertial_pos_body(3) =  ppval(interpolators.(BODY_str).spline{3}, inertial_t);

    rtbp_pos_body =  (-C\(inertial_pos_body-b)/k); 
    x_sb = x(1:3)-rtbp_pos_body;
    rho_sb = (x(1)-rtbp_pos_body(1))^2+(x(2)-rtbp_pos_body(2))^2+(x(3)-rtbp_pos_body(3))^2;
    rho_sb3 = rho_sb*sqrt(rho_sb);
    synodic_acc_bodies = -mu_body*x_sb/rho_sb3 + synodic_acc_bodies; %primaries + bodies contribution to synodic acceleration
     %now compute the gravitational_acceleration_synodic(all bodies)
     %the function gravitational_acceleration gives acc in the inertial frame, so we need to divide the GMj contribution by (GM_1+GM_2) to get mu_j
    acc_bodies(:, i) =synodic_acc_bodies;
    
end

Delta_Omega = sum(acc_bodies,2) + synodic_acc_primaries;

else
    Delta_Omega = synodic_acc_primaries;
    fprintf('NO BODIES\n')
end
     %now compute the gravitational_acceleration_synodic(all bodies)
     %the function gravitational_acceleration gives acc in the inertial frame, so we need to divide the GMj contribution by (GM_1+GM_2) to get mu_j


xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

xdot(4:6) = first + second*x(4:6) + third*x(1:3) + b13*Delta_Omega;


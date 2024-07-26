global mu
%------- THIS MAIN IS TO SHOW THE FRAME OF REFERENCE CHANGES PERFORMED BY <spice_to_rtbp.m> -------


% --------------- LOAD KERNELS -------------
META = '/Users/luigigisolfi/Documents/mcodes/mcodes/kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------- %

% ----- Define the parameters to enter the spice_to_rtbp conversion function -----
ti = '2000 SEP 1 TDB';
tf = '2000 OCT 1 TDB';
FRAME = 'ECLIPJ2000';
OBSERVER = 'SOLAR SYSTEM BARYCENTER';
BODIES = [{'MERCURY'}, {'VENUS'}];
PRIMARIES = [{'SUN'}, {'EARTH-MOON BARYCENTER'}];
L = get_L(FRAME, PRIMARIES, OBSERVER);
% ------------------------------------------------------------------------------

% ------ Retrieve SPICE positions and epochs in the RTBP-like coordinates --------
xi_L1 = [(mu-1)+sqrt(mu/3), 0, 0, 0.02, 0.02, 0.02];
hmin=1.e-6; hmax=1.e0; tol=1.e-10; h=0.01;
ntall=1; iorb=-1;

[rtbp_times, rtbp_state_primaries, body_1, body_2, BODIES, rtbp_state_bodies, interpolators]  = new_spice_to_rtbp(ti,tf, FRAME, OBSERVER, PRIMARIES, L, BODIES);

% % Check the first interpolated position against the actual data
% disp('First interpolated vel:');
% disp(planet_positions.("MERCURY")(:, 1));
% disp('Actual first vel:');
% disp(rtbp_state_bodies.MERCURY(:, 1));
% figure
% hold on
% scatter(planet_positions.("MERCURY")(1,:),planet_positions.("MERCURY")(2,:), 'DisplayName', 'Mercury')
% 
% scatter(planet_positions.("VENUS")(1,:),planet_positions.("VENUS")(2,:), 'DisplayName', 'Venus')
% legend('show')

rtbp_ti = rtbp_times(1);
rtbp_tf = rtbp_times(end);

x_b = rtbp_pos_bodies(:,:,1)

% ---------------------------------------------------------------------------------

% ----- Propagate your own spacecraft. this one starts close to L1 ------

[tfa, xf] = propTITF(rtbp_ti,xi_L1,rtbp_tf,@rtbp,hmin,hmax,tol,1);
xf_pos = xf(:,1:3);
xf_vel = xf(:,4:6);
% ---------------------------------------------------------------------------------

% ----- PLOT IN THE SYNODIC FRAME OF THE 2 PRIMARIES ------- %
figure
axis equal
hold on
%scatter(0,0)

scatter3(mu,0,0, 'filled', 'DisplayName',body_1)
scatter3(mu-1,0,0, 'filled', 'DisplayName',body_2)

%text(mu,0,0, body_1)
%text(mu-1,0,0,body_2)
plot3(rtbp_pos_1(1,:), rtbp_pos_1(2,:), rtbp_pos_1(3,:), 'DisplayName', 'Sun')
plot3(rtbp_pos_2(1,:), rtbp_pos_2(2,:), rtbp_pos_2(3,:), 'DisplayName', 'Earth-Moon Barycenter')

for j = 1:length(BODIES)
    BODY = BODIES{j};
    plot3(rtbp_(1,:,j),rtbp_pos_bodies(2,:,j),rtbp_pos_bodies(3,:,j), 'DisplayName', BODY)
end

plot3(xf_pos(:,1), xf_pos(:,2), xf_pos(:,3), 'DisplayName', 'Orbit')

%plot3(rtbp_rtbp_pos_1(1,:), rtbp_rtbp_pos_1(2,:), rtbp_rtbp_pos_1(2,:))
%plot3(rtbp_rtbp_pos_2(1,:), rtbp_rtbp_pos_2(2,:), rtbp_rtbp_pos_1(3,:))
str_title = [body_1, ' / ', body_2, ' SYNODIC REFERENCE FRAME'];
title(str_title)
xlabel('x (AU)')
ylabel('y (AU)')
zlabel('z (AU)')
legend('show')
%--------------------------------------------------------------------------%

% ----- PLOT IN THE ECLIPTIC FRAME ------- %
%figure
%hold on
%axis equal
%plot3(pos_1(1,:), pos_1(2,:), pos_1(3,:), 'DisplayName', 'Sun')
%plot3(pos_2(1,:), pos_2(2,:), pos_2(3,:), 'DisplayName', 'Earth-Moon Barycenter')

%for j = 1:length(BODIES)
 %   BODY = BODIES(j);
 %   plot3(states_bodies(1,:,j),states_bodies(2,:,j), states_bodies(3,:,j))
%end
%str_title = [body_1, ' / ', body_2, ' ECLIPTIC REFERENCE FRAME'];
%title(str_title)
%xlabel('x (km')
%ylabel('y (km)')
%zlabel('z (km)')
%legend('show')
% -------------------------------------------- %

cspice_kclear()



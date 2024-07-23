%-------------------------------------------------------------------
% Example of a propagation of an intial RTBP state (pos+vel) using
% the RK45 integrator and to stop at a given final DT.
%-------------------------------------------------------------------
clear all; close all; clc;
global mu;
%mu=0.0122; %for Earth-Moon
mu = 3.0035e-06 %for Earth Sun
%--------------------------------------------------------------------------
 ti=0;
 r_geo = 2.818791946308725e-04 %corresponding to 42000 km orbit, GEO
 xi_circular = [mu-1-r_geo 0 0 0 sqrt(mu/r_geo) 0]
 %xi_L1=[-0.866                  0  -0.000291323690820                   0  +0.083326303745315                   0];
 %xi_L3 = [1.02 0.000000000E+00 .0000000000E+00 ...
 %    0.0000000000E+00  -0.111258721E-00 .0000000000E+00];
 %xi_L4 = [0.54+mu-1 sqrt(3)/2 0 0 0.05 0]
 
 hmin=1.e-6; hmax=1.e0; tol=1.e-10; h=0.001;
 ntall=1; iorb=-1;
 %[tf_rtbp, xf_rtbp] = propTallSec(ti,xi,h,hmin,hmax,tol,@rtbp,@gdgYSec,ntall,iorb);
 %[tf_bfbp, xf_bfbp] = propTallSec(ti,xi,h,hmin,hmax,tol,@bfbp,@gdgYSec,ntall,iorb);
%-------------------------------------------------------------------------
 %cji=cjrtbp(xi_L3,mu);
 tf = 1
 %tf_rtbp = 2*tf_rtbp(end);
 %tf_bfbp = 2*tf_bfbp(end);
 %xi_rtbp = x_sol_rtbp(1);
 %xi_bfbp = x_sol_bfbp(1);

 [tfa_bfbp_circ, xf_bfbp_circ] = ode78(@rtbp, [0,0.1], xi_circular)
 axis equal
 hold on
 scatter3(xf_bfbp_circ(1,1), xf_bfbp_circ(1,2), xf_bfbp_circ(1,3))
 plot3(xf_bfbp_circ(:,1), xf_bfbp_circ(:,2),xf_bfbp_circ(:,3))
 scatter3(mu-1, 0,0)

 filename_end = 'circular_orbit_Earth_Sun'
data_new = [tfa_bfbp_circ, xf_bfbp_circ(:,1), xf_bfbp_circ(:,2),xf_bfbp_circ(:,3),xf_bfbp_circ(:,4), xf_bfbp_circ(:,5),xf_bfbp_circ(:,6)]; % Combine time with x, y, z coordinates
fid = fopen(filename_end, 'w');
fprintf(fid, 'Time (s)\tX (km)\tY (km)\tZ (km)\tvx (km)\tvy (km)\tvz (km)\n');
fprintf(fid, '%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\t%.12f\n', data_new.');
fclose(fid);
disp(['Orbit coordinates saved to ' filename_end]);
return
 [tfa_rtbp_L1, xf_rtbp_L1] = propTITF(ti,xi_L1,tf,@rtbp,hmin,hmax,tol,1);
 [tfa_etbp_L1, xf_etbp_L1] = propTITF(ti,xi_L1,tf,@etbp,hmin,hmax,tol,1);
 
 [tfa_bfbp_L3, xf_bfbp_L3] = propTITF(ti,xi_L3,tf,@bfbp,hmin,hmax,tol,1);
 [tfa_rtbp_L3, xf_rtbp_L3] = propTITF(ti,xi_L3,tf,@rtbp,hmin,hmax,tol,1);
 [tfa_etbp_L3, xf_etbp_L3] = propTITF(ti,xi_L3,tf,@etbp,hmin,hmax,tol,1);

 [tfa_bfbp_L4, xf_bfbp_L4] = propTITF(ti,xi_L4,tf,@bfbp,hmin,hmax,tol,1);
 [tfa_rtbp_L4, xf_rtbp_L4] = propTITF(ti,xi_L4,tf,@rtbp,hmin,hmax,tol,1);
 [tfa_etbp_L4, xf_etbp_L4] = propTITF(ti,xi_L4,tf,@etbp,hmin,hmax,tol,1);

m3 = 332998; %Msun*mu/m* = Msun*mu/(ME+MM)
a3 = 149597870.691/384400; %AU/EM distance ratio
n3 = sqrt((1 + m3)/a3^3);
w3 = n3 - 1;
theta3 = w3*tfa_bfbp_L1;
Omega = 0; %values of Omega can be retrieved applying a linear interpolation to the dates corresponding to past eclipses ????!!!!. 
inc = deg2rad(-5.16); %inc of the sun wrt EM system is -5.16°. Often it can be assumed to be 0. In that case: "Sun in Plane" B4BP
x3 = a3*(cos(theta3 - Omega)*cos(Omega) - sin(theta3 - Omega)*sin(Omega)*cos(inc)); %Sun coordinates as in Boudad_Purdue
y3 = a3*(cos(theta3 - Omega)*sin(Omega) + sin(theta3 - Omega)*cos(Omega)*cos(inc)); %Sun coordinates as in Boudad_Purdue
z3 = a3*sin(theta3 - Omega)*sin(inc);
%--------------------------------------
% Computing Body-Sun Position for BFBP
% -------------------------------------
r32_L1 = (xf_bfbp_L1(1) -x3).^2 + (xf_bfbp_L1(2) - y3).^2 + (xf_bfbp_L1(3) - z3).^2; %r32: square of distance to Sun (3rd primary)
r32_L3 = (xf_bfbp_L3(1) -x3).^2 + (xf_bfbp_L3(2) - y3).^2 + (xf_bfbp_L3(3) - z3).^2;
r32_L4 = (xf_bfbp_L4(1) -x3).^2 + (xf_bfbp_L4(2) - y3).^2 + (xf_bfbp_L4(3) - z3).^2;%r32: square of distance to Sun (3rd primary)
r31_L1 = sqrt(r32_L1);
r31_L3 = sqrt(r32_L3);
r31_L4 = sqrt(r32_L4);

fine_L1 = min(length(xf_bfbp_L1), min(length(xf_rtbp_L1), length(xf_etbp_L1)))
fine_L3 = min(length(xf_bfbp_L3), min(length(xf_rtbp_L3), length(xf_etbp_L3)))
fine_L4 = min(length(xf_bfbp_L4), min(length(xf_rtbp_L4), length(xf_etbp_L4)))

abs_vals_L1 = sqrt((xf_bfbp_L1(1:fine_L1,1)- xf_rtbp_L1(1:fine_L1,1)).^2 + (xf_bfbp_L1(1:fine_L1,2)- xf_rtbp_L1(1:fine_L1,2)).^2);
abs_vals_L3 = sqrt((xf_bfbp_L3(1:fine_L3,1)- xf_rtbp_L3(1:fine_L3,1)).^2 + (xf_bfbp_L3(1:fine_L3,2)- xf_rtbp_L3(1:fine_L3,2)).^2);
abs_vals_L4 = sqrt((xf_bfbp_L4(1:fine_L4,1)- xf_rtbp_L4(1:fine_L4,1)).^2 + (xf_bfbp_L4(1:fine_L4,2)- xf_rtbp_L4(1:fine_L4,2)).^2);

abs_vals_L1_E = sqrt((xf_etbp_L1(1:fine_L1,1)- xf_rtbp_L1(1:fine_L1,1)).^2 + (xf_etbp_L1(1:fine_L1,2)- xf_rtbp_L1(1:fine_L1,2)).^2);
abs_vals_L3_E = sqrt((xf_etbp_L3(1:fine_L3,1)- xf_rtbp_L3(1:fine_L3,1)).^2 + (xf_etbp_L3(1:fine_L3,2)- xf_rtbp_L3(1:fine_L3,2)).^2);
abs_vals_L4_E = sqrt((xf_etbp_L4(1:fine_L4,1)- xf_rtbp_L4(1:fine_L4,1)).^2 + (xf_etbp_L4(1:fine_L4,2)- xf_rtbp_L4(1:fine_L4,2)).^2);

x3_vals_L1 = x3(1:fine_L1);
x3_vals_L3 = x3(1:fine_L3)
x3_vals_L4 = x3(1:fine_L4)

%-----------------------------------------
% Jacobi Computation (not implemented yet)
%-----------------------------------------
cjf=cjrtbp(xf_rtbp_L3,mu);
fprintf('Jacobi ct: ini, final and diff: %e %e %5.2e\n',cji,cjf,cjf-cji);

%----------------------
% Plots
%----------------------
collinear_lagrange;
plot(xf_bfbp_L1(:,1), xf_bfbp_L1(:,2), 'red')
plot(xf_rtbp_L1(:,1), xf_rtbp_L1(:,2), 'blue')
plot(xf_etbp_L1(:,1), xf_etbp_L1(:,2), 'green')
plot(xf_bfbp_L3(:,1), xf_bfbp_L3(:,2), 'red')
plot(xf_rtbp_L3(:,1), xf_rtbp_L3(:,2), 'blue')
plot(xf_etbp_L3(:,1), xf_etbp_L3(:,2), 'green')
plot(xf_bfbp_L4(:,1), xf_bfbp_L4(:,2), 'red')
plot(xf_rtbp_L4(:,1), xf_rtbp_L4(:,2), 'blue')
plot(xf_etbp_L4(:,1), xf_etbp_L4(:,2), 'green')


yline(0, 'k')
legend('$Earth$','$Moon$', '$L_1$', '$L_2$', '$L_3$', '$L_4$', '$L_5$','$Bicircular$', '$R3BP$', '$E3BP$','Interpreter', 'latex')

figure
hold on
title('Just an Insight on Sun Relative Position')
quiver(xf_bfbp_L1(1:fine_L1,1), xf_bfbp_L1(1:fine_L1,2),x3(1:fine_L1)- xf_bfbp_L1(1:fine_L1,1), y3(1:fine_L1) -xf_bfbp_L1(1:fine_L1,2), 0.2)
quiver(xf_bfbp_L3(1:fine_L3,1), xf_bfbp_L3(1:fine_L3,2),x3(1:fine_L3)- xf_bfbp_L3(1:fine_L3,1), y3(1:fine_L3) -xf_bfbp_L3(1:fine_L3,2), 0.2)
quiver(xf_bfbp_L4(1:fine_L4,1), xf_bfbp_L4(1:fine_L4,2),x3(1:fine_L4)- xf_bfbp_L4(1:fine_L4,1), y3(1:fine_L4) -xf_bfbp_L4(1:fine_L4,2), 0.2)

plot(xf_bfbp_L1(:,1), xf_bfbp_L1(:,2), 'red')
plot(xf_rtbp_L1(:,1), xf_rtbp_L1(:,2), 'blue')
plot(xf_bfbp_L3(:,1), xf_bfbp_L3(:,2), 'red')
plot(xf_rtbp_L3(:,1), xf_rtbp_L3(:,2), 'blue')
plot(xf_bfbp_L4(:,1), xf_bfbp_L4(:,2), 'red')
plot(xf_rtbp_L4(:,1), xf_rtbp_L4(:,2), 'blue')
xlabel('x')
ylabel('y')
legend('$Body-Sun Pos. Vec. L1$','$Body-Sun Pos. Vec. L1$','$Bicircular, L1$', '$R3BP, L1$', '$Bicircular, L3$', '$R3BP, L3$', 'Interpreter', 'latex')

figure 
hold on
title('Comparison Between Different Models (Earth-Moon)')
plot(tfa_bfbp_L1(1:fine_L1), abs_vals_L1, 'red')
plot(tfa_etbp_L1(1:fine_L1), abs_vals_L1_E, 'green')
plot(tfa_bfbp_L3(1:fine_L3), abs_vals_L3, 'red')
plot(tfa_bfbp_L4(1:fine_L4), abs_vals_L4, 'red')
plot(tfa_etbp_L3(1:fine_L3), abs_vals_L3_E, 'green')
plot(tfa_etbp_L4(1:fine_L4), abs_vals_L4_E, 'green')

text(tfa_bfbp_L1(fine_L1), abs_vals_L1(fine_L1), 'L1')
text(tfa_bfbp_L3(fine_L3), abs_vals_L3(fine_L3), 'L3')
text(tfa_bfbp_L4(fine_L4), abs_vals_L4(fine_L4), 'L4')

text(tfa_etbp_L1(fine_L1), abs_vals_L1_E(fine_L1), 'L1')
text(tfa_etbp_L3(fine_L3), abs_vals_L3_E(fine_L3), 'L3')
text(tfa_etbp_L4(fine_L4), abs_vals_L4_E(fine_L4), 'L4')

xlabel('Time')
ylabel('Absolute Difference')
legend('BFBP-RTBP', 'ETBP-RTBP')

%axis equal

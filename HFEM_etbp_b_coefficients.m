function xdot = HFEM_etbp_b_coefficients(t,x)

global PRIMARIES
global BODIES
global interpolators
global mu
global avg_e
global L
global r
global nu


rp = [(mu);0;0]; %should I normalize the lengths so as to  have homogenous comparison in the adimensional coordinates?
rs = [(mu-1);0;0];

b1=0;
b2=0;
b3=0;
b4=-avg_e*sin(nu)/(2*sqrt(1+avg_e*cos(nu)))
b5 = 2*(sqrt(1+avg_e*cos(nu)))
b6=0
b7 = 1;
b8=0;
b9=0;
b10 = 1;
b11=0;
b12=-avg_e*cos(nu);
b13 = 1;

xdot = zeros(7,1);

first = [b1;b2;b3];
second = [b4,b5,0;-b5,b4,b6;0,-b6,b4];
third = [b7,b8,b9;-b8,b10,b11;b9,-b11,b12];

rtbp_pos_body_1 = rp;
rtbp_pos_body_2 = rs;

% fprintf('ELLIPTICAL POSITIONS PRIMARIES NORMALIZED TO RTBP LENGTH UNITS')
% rtbp_pos_body_1
% rtbp_pos_body_2


x_s1 = x(1:3)-rtbp_pos_body_1;
x_s2 = x(1:3)-rtbp_pos_body_2;

rho_s1 = (x(1)-rtbp_pos_body_1(1))^2+(x(2))^2+(x(3))^2;
rho_s2 = (x(1)-rtbp_pos_body_2(1))^2+(x(2))^2+(x(3))^2;
rho_s13 = rho_s1*sqrt(rho_s1);
rho_s23 = rho_s2*sqrt(rho_s2);

synodic_acc_primaries = -(1-mu)*x_s1/rho_s13 - (mu)*x_s2/rho_s23; %primaries contribution to synodic acceleration

Delta_Omega = synodic_acc_primaries;

xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4:6) = first + second*x(4:6) + third*x(1:3) + b13*Delta_Omega;
xdot(7) = sqrt(1+avg_e*cos(nu));



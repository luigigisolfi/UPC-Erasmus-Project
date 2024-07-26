function f = bfbp(t, x)
%------------------------------------------------------------------------
% Spatial Bicircular Restricted BRFBP vectorfield (as in Boudad_Purdue).
% Large Mass, 1-mu (Earth)    to the right of the origin at (mu, 0, 0)
% Small Mass, mu   (Moon) to the left at (mu-1, 0, 0).
% Big Distant Mass (Sun) in circular orbit around them.
%
%                      L5
%
% L2 -- Moon -- L1 ----------------- Earth --------------- L3
%
%                      L4
%
% Input variables: t (time), x (3D state, pos+vel) and mu (mass param.)
% Output: f (vectorfield)
%-----------------------------------------------------------------------
global mu;
global m3;
global a3;
global n3;

m3 = 332998; %Msun*mu/m* = Msun*mu/(ME+MM)
a3 = 149597870.691/384400; %AU/EM distance ratio
n3 = sqrt((1 + m3)/a3^3);
r12= (x(1)-mu+1)^2 + x(2)^2 + x(3)^2; % r12: square of distance to P 
r22= (x(1)-mu)^2 + x(2)^2 + x(3)^2;   % r22: square of distance to S

w3 = n3 - 1;
theta3 = w3*t;
Omega = 0; %values of Omega can be retrieved applying a linear interpolation to the dates corresponding to past eclipses ????!!!!. 
inc = deg2rad(-5.16); %inc of the sun wrt EM system is -5.16°. Often it can be assumed to be 0. In that case: "Sun in Plane" B4BP
x3 = a3*(cos(theta3 - Omega)*cos(Omega) - sin(theta3 - Omega)*sin(Omega)*cos(inc)); %Sun coordinates as in Boudad_Purdue
y3 = a3*(cos(theta3 - Omega)*sin(Omega) + sin(theta3 - Omega)*cos(Omega)*cos(inc)); %Sun coordinates as in Boudad_Purdue
z3 = a3*sin(theta3 - Omega)*sin(inc);
r32 = (x(1) -x3)^2 + (x(2) - y3)^2 + (x(3) - z3)^2; %r32: square of distance to Sun (3rd primary)

r13=r12*sqrt(r12);
r23=r22*sqrt(r22);
r33= r32*sqrt(r32);

Ox = x(1) - ((1-mu)*(x(1)-mu)/r23 + mu*(x(1)-mu+1)/r13);
Oy = x(2) - ((1-mu)* x(2)    /r23 + mu* x(2)      /r13);
Oz =      - ((1-mu)* x(3)    /r23 + mu* x(3)      /r13);

xdot = zeros(6,1);
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

xdot(4) = 2*x(5) + Ox - (m3*(x(1)-x3))/r33 - (m3*x3)/a3^3 ;
xdot(5) =-2*x(4) + Oy - (m3*(x(2)-y3))/r33 - m3*y3/a3^3;
xdot(6) = Oz - m3*x(3)/r33 - m3*z3/a3^3;

f = xdot;

end

function f = hill3b(t, x)
%------------------------------------------------------------------------
% Hill 3 Body Problem vectorfield.
% Small mass (Moon) at the origin of the reference system
% Big Mass (Earth) at a distance of 1 from it
%
%                      L5
%
% L2 -- Origin (Moon) -- L1 ----------------- Earth --------------- L3
%
%                      L4
%
% Input variables: t (time), x (3D state, pos+vel) and mu (mass param.)
% Output: f (vectorfield)
%-----------------------------------------------------------------------
global mu;
% The Hill Problem is centered on the Moon
% so if we want to compare the results with the other models 
% (that, instead are centered on the EM barycenter, 
% we need a translation of the initial x condition on the x-axis:
x(1) = x(1) - mu + 1;

r= (x(1))^2 + x(2)^2 + x(3)^2;
r3=r*sqrt(r);

xdot = zeros(6,1);
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

xdot(4) = 3*x(1) - mu*x(1)/r3 + 2*x(5) ;
xdot(5) =-2*x(4) - mu*x(2)/r3 ;
xdot(6) = -mu*x(3)/r3 - x(3);

f = xdot;

end

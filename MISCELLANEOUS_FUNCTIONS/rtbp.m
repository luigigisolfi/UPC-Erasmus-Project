function f = rtbp(t, x)
%------------------------------------------------------------------------
% Spatial Circular RTBP vectorfield.
% Large Mass, 1-mu (Earth)    to the right of the origin at (mu, 0, 0)
% Small Mass, mu   (Moon) to the left at (mu-1, 0, 0).
%
%                      L5
%
% L2 -- Moon-- L1 ----------------- Earth --------------- L3
%
%                      L4
%
% Input variables: t (time), x (3D state, pos+vel) and mu (mass param.)
% Output: f (vectorfield)
%-----------------------------------------------------------------------
global mu;

r12= (x(1)-mu)^2 + x(2)^2 + x(3)^2;   % r22: square of distance to S
r22= (x(1)-mu+1)^2 + x(2)^2 + x(3)^2; % r12: square of distance to P 
r13=r12*sqrt(r12);
r23=r22*sqrt(r22);

Ox = x(1) - ((1-mu)*(x(1)-mu)/r13 + mu*(x(1)-mu+1)/r23);
Oy = x(2) - ((1-mu)* x(2)    /r13 + mu* x(2)      /r23);
Oz =      - ((1-mu)* x(3)    /r13 + mu* x(3)      /r23);


xdot = zeros(6,1);
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

xdot(4) = 2*x(5) + Ox; 
xdot(5) =-2*x(4) + Oy;
xdot(6) = Oz;

f = xdot;

end
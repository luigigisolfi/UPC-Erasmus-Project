function [xf] = new_propTITF_vfield(ti, xi, tf, vfield, hmin, hmax, tol)
% Propagate the state from ti to tf using the vector field 'vfield'
% with adaptive step size control using the ode78 integrator

%options = odeset('MaxStep', h, 'InitialStep', h);
[~, x] = ode45(vfield, [ti,tf], xi);
% Final state after propagation
xf = x(end);
end

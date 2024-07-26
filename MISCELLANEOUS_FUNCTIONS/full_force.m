function f = full_force(t, x)

%t and x are given in the INERTIAL reference frame. Therefore, this vector
%field is technically "different" than the other ones (RTBP, HILL3B,
%ecc...) because it does not accept rtbp-adimensional coordinates, but only
%physical ones.
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
global BODIES
global FRAME
global OBSERVER
global PRIMARIES
%--------------------------------------------------------------------------
% Suppose we have (x0,y0,z0,vx0,vy0,vz0) of the particle
% (that is, state vector wrt the Earth (could be another body too)
% 1) Retrieve all relative |rbody - rj| with SPICE
% 2) Compute the total acceleration as sum of N acting bodies (planets,
% minor planets) on the particle
% 3) Perform a Runge Kutta integration step and retrieve (x1,y1,z1, ...)
% 4) Repeat
%
% Please note that we might need a function that changes coordinates
% from whatever RF we have them defined in SPICE, to the EM ROTATING FRAME
%--------------------------------------------------------------------------


% ----- Load required Kernels and ephemeris file ------ %
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ------------------------------------------------------ %

all_bodies = [PRIMARIES, BODIES];
n_bodies = length(all_bodies);
acc_list = zeros(3, n_bodies);
xdot = zeros(6,1);

t_utc = cspice_et2utc(t, 'C', 0);

%fprintf('Epoch: %s\n', t_utc);
parfor i = 1:n_bodies % PARALLELIZED COMPUTATION ON ALL BODIES
    BODY = all_bodies{i};
    GM_body = get_GM_body(BODY);
    acc_list(:, i) = gravitational_acceleration(t, x, BODY, GM_body);
end

    % Sum all accelerations to get the total acceleration

atot = sum(acc_list, 2);

xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);
xdot(4) = atot(1);
xdot(5) = atot(2);
xdot(6) = atot(3);

f = xdot;

cspice_kclear()

end

function f = distance_EM(et)

META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels

% Define the observation time
% Use a standard date-time format (UTC) and convert to ephemeris time

% Define the observer and target
observer = 'EARTH';
target = 'MOON';

% Compute the state vector of the Moon relative to the Earth at the specified time
[state, lt] = cspice_spkgeo(cspice_bodn2c(target), et, 'J2000', cspice_bodn2c(observer));

% Extract the position vector (first 3 elements of state)
position = state(1:3);

% Compute the distance (magnitude of the position vector)
f = norm(position);

% Unload kernels after use
cspice_kclear;

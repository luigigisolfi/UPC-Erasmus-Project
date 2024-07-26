META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% Define the time window
start_time = '2000 JAN 01 00:00:00';
end_time = '2000 JAN 02 00:00:00';

% Convert time to ephemeris time (ET)
et_start = cspice_str2et(start_time);
et_end = cspice_str2et(end_time);

% Define the step size (1 minute)
step = 60;

% Set up search parameters
target = 'MOON';
abcorr = 'NONE';
obsrvr = 'EARTH';
relate = 'LOCMIN'; % Minimum illumination condition (eclipse)

% Search for occultation
result = cspice_occult('SUN', 'POINT', ' ', 'MOON', 'ELLIPSOID', 'ITRF93', 'NONE', 'EARTH', [et_start,et_end, step]);

% Check results
if ~isempty(result)
    for i = 1:length(result)
        eclipse_time = cspice_et2utc(result(i), 'C', 3);
        fprintf('Lunar eclipse occurred at: %s\n', eclipse_time);
    end
else
    fprintf('No lunar eclipse found in the specified time range.\n');
end

% Unload kernels
cspice_kclear;

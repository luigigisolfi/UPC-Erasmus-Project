function n = get_n(primaries, L)

% ----------------------------------------------------------------------%
% This function uses the kernel Gravity.tpc and computes the mean motion
% of a given system made of two primaries
% ----------------------------------------------------------------------%

% ----------------LOADING KERNELS --------------------------------------%
META = 'kernels_to_load.tm'; %initialize required kernels
cspice_furnsh(META); %furnish kernels
% ----------------------------------------------------------------------%

% -----------------------------------------------%
% The user is allowed to input the primaries in every order 
% This is taken care of in the following, computing the max_GM
body_1 = primaries(1);
body_2 = primaries(end);
GM_1 = cspice_bodvrd(body_1, 'GM', 1);
GM_2 = cspice_bodvrd(body_2, 'GM', 1);
n = sqrt((GM_1 + GM_2)/L^3);
% ----------------------------------------------------------------------%
cspice_kclear()
% ----------------------------------------------------------------------%

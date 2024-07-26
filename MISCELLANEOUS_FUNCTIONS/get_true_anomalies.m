function nu_interpolator = get_true_anomalies(e,n,tspan, nu0)
% Define the differential equation
dydt = @(t, f) sqrt(1 + e * cos(f));

% Set the initial conditions and time span

% Use ode45 to solve the differential equation
[t, f] = ode78(dydt, tspan, nu0);

true_anomalies = struct();

    true_anomalies.x = t;
    true_anomalies.y = f;
    true_anomalies.spline = cell(1,1); 

    true_anomalies.spline{1} = spline(t, f);

nu_interpolator = true_anomalies;

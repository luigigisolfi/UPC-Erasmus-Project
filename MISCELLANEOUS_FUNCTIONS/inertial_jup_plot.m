global interpolators
global n
list_t = linspace(0,2*pi/n, 10000)
inertial_pos_j =zeros(3,10000);
inertial_pos_e =zeros(3,10000);
for i = 1:10000
    t_ = list_t(i);
    inertial_pos_j(1,i) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{1}, t_);
    inertial_pos_j(2,i) =  ppval(interpolators.('JUPITER_BARYCENTER').spline{2}, t_);
    inertial_pos_e(1,i) =  ppval(interpolators.('EARTH_MOON_BARYCENTER').spline{1}, t_);
    inertial_pos_e(2,i) =  ppval(interpolators.('EARTH_MOON_BARYCENTER').spline{2}, t_);
end


figure
hold on
plot(inertial_pos_j(1,:), inertial_pos_j(2,:))
plot(inertial_pos_e(1,:), inertial_pos_e(2,:))


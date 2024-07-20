
G = 6.67384e-11;
mass_earth = 5.97219e24;
radius_earth = 6.371e6;
mass_moon = 7.3477e22;
radius_moon = 1.7371e6;
mass_vehicle = 10e3;
T_earth_moon = 27.321662*24*3600;
a_earth_moon = ((T_earth_moon/2/%pi)^2*G*(mass_earth+mass_moon))^(1/3.0);

g = G*mass_earth/radius_earth^2;

r_perigee = radius_earth+500e3;
r_apogee = a_earth_moon;

e = (r_apogee-r_perigee)/(r_apogee+r_perigee);

v_circular = sqrt(G*(mass_earth+mass_vehicle)/(radius_earth+500e3));

v_approach = sqrt(G*(mass_earth+mass_vehicle)*(1+e)/(radius_earth+500e3));

a_approach = (r_perigee+r_apogee)/2;

T_approach = 2*%pi*sqrt(a_approach^3/G/(mass_earth+mass_vehicle));

thrust = 2*g;

delta_v = v_approach - v_circular;

T_thrust = delta_v/thrust;


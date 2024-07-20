G = 6.67384e-11;
M = 5.97219e24;
m = 1e3;
radius_earth = 6.371e6;
h = 500e3;
a = radius_earth + h;
theta_dot_circular = sqrt(G*(M+m)/a^3);
v = theta_dot_circular*a;
acc = 9.81*2;

r_maneuver = v^2/acc;


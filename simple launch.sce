
G = 6.67384e-11;
M = 5.97219e24;
r = 6.371e6;

g = G*M/r^2;

gt = 2*g;
h = 500e3;

t_thrust = sqrt(h/(0.5*(gt-g)^2/g + 0.5*(gt-g)));
vt = (gt-g)*t_thrust;

t_coast = vt/g;

t_total = t_thrust + t_coast;


G = 6.67384e-11;
M = 5.97219e24;
ecc = 0.5;
function a=r_ddot(theta_dot,r,GM)
    a = theta_dot^2*r - GM/(r^2);
endfunction

function a=theta_ddot(theta_dot,r,r_dot)
    a = -2.0*theta_dot*r_dot/r;
endfunction

function ydot = ode_func(t,y)
//
// y = [r,theta,r_dot,theta_dot]
//
    ydot = [y(3);y(4);y(4)^2*y(1) - G*M/y(1)^2;-2.0*y(4)*y(3)/y(1)];
endfunction

function [theta1,theta_dot1,r1,r_dot1]=integrate_euler(theta0,theta_dot0,r0,r_dot0,GM,dt)
    kr_ddot = r_ddot(theta_dot0,r0,GM);
    ktheta_ddot = theta_ddot(theta_dot0,r0,r_dot0);
    
    r_dot1 = r_dot0 + kr_ddot*dt;
    r1 = r0 + r_dot0*dt;
    
    theta_dot1 = theta_dot0 + ktheta_ddot*dt;
    theta1 = theta0 + theta_dot0*dt;
endfunction

function ydot = rk4_f(y,GM)
    ydot = [y(2),-2*y(4)*y(2)/y(3),y(4),y(2)^2*y(3)-GM/y(3)^2];
endfunction

function [theta1,theta_dot1,r1,r_dot1]=integrate_rk4(theta0,theta_dot0,r0,r_dot0,GM,dt)
    y0 = [theta0,theta_dot0,r0,r_dot0];
    k1 = rk4_f(y0,GM);
    k2 = rk4_f(y0 + dt/2*k1,GM);
    k3 = rk4_f(y0 + dt/2*k2,GM);
    k4 = rk4_f(y0 + dt*k3,GM);
    y1 = y0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    theta1 = y1(1);
    theta_dot1 = y1(2);
    r1 = y1(3);
    r_dot1 = y1(4);
endfunction

function Ek = Ekinetic(theta_dot,r,r_dot,m)
    vtheta = r*theta_dot;
    v2 = vtheta^2 + r_dot^2;
    Ek = 0.5*m*v2;
endfunction

function Ep = Epotential(r,m,GM)
    Ep = -GM*m/r;
endfunction

m = 1e3;
radius_earth = 6.371e6;
h = 500e3;
a = radius_earth + h;
theta_dot_circular = sqrt(G*M/a^3);
T = 2*%pi/theta_dot_circular;
Npoints = 1000;
dt = T/Npoints;
radius_periapsis = a*(1.0-ecc);
v_periapsis = sqrt(G*M*(2/radius_periapsis-1/a));
theta_dot = v_periapsis/radius_periapsis;
theta = 0.0;
r_dot = 0.0;
r = radius_periapsis;
y0 = [r;theta;r_dot;theta_dot];
Ek0 = Ekinetic(theta_dot,r,r_dot,m);
Ep0 = Epotential(r,m,G*M);
Et0 = Ek0 + Ep0;
clear x;
for i=1:Npoints+1
    x(i,1:2) = [cos(theta)*r,sin(theta)*r];
    [theta,theta_dot,r,r_dot] = integrate_rk4(theta,theta_dot,r,r_dot,G*M,dt);
end

Ek1 = Ekinetic(theta_dot,r,r_dot,m);
Ep1 = Epotential(r,m,G*M);
Et1 = Ek1 + Ep1;

disp((Et1-Et0)/Et0);

clf
plot2d(x(:,1),x(:,2),frameflag=4);

t_eval = linspace(0.0,T,Npoints);
t0 = 0.0;
y_result = ode("rk",y0,t0,t_eval,ode_func);

clear x;

for i=1:Npoints
    x(i,1:2) = [cos(y_result(2,i))*y_result(1,i),sin(y_result(2,i))*y_result(1,i)]
end

plot2d(x(:,1),x(:,2),frameflag=4,style=[2]);




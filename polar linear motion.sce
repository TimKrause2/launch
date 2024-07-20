function a=r_ddot(theta_dot,r)
    a = theta_dot^2*r;
endfunction

function a=theta_ddot(theta_dot,r,r_dot)
    a = -2.0*theta_dot*r_dot/r;
endfunction

function [theta1,theta_dot1,r1,r_dot1]=integrate_euler(theta0,theta_dot0,r0,r_dot0,dt)
    kr_ddot = r_ddot(theta_dot0,r0);
    ktheta_ddot = theta_ddot(theta_dot0,r0,r_dot0);
    
    r_dot1 = r_dot0 + kr_ddot*dt;
    r1 = r0 + r_dot1*dt;
    
    theta_dot1 = theta_dot0 + ktheta_ddot*dt;
    theta1 = theta0 + theta_dot1*dt;
endfunction

function ydot = rk4_df(y)
    ydot = [y(2),-2*y(4)*y(2)/y(3),y(4),y(2)^2*y(3)];
endfunction

function [theta1,theta_dot1,r1,r_dot1]=integrate_rk4(theta0,theta_dot0,r0,r_dot0,dt)
    y0 = [theta0,theta_dot0,r0,r_dot0];
    k1 = rk4_df(y0);
    k2 = rk4_df(y0 + dt/2*k1);
    k3 = rk4_df(y0 + dt/2*k2);
    k4 = rk4_df(y0 + dt*k3);
    y1 = y0 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    theta1 = y1(1);
    theta_dot1 = y1(2);
    r1 = y1(3);
    r_dot1 = y1(4);
endfunction

v_vec = [9.0;10.0];
r_vec = [-10.0;-10.0];

r = sqrt(sum(r_vec.^2));
theta = atan(r_vec(2),r_vec(1));
r_dot = v_vec(1)*cos(theta) + v_vec(2)*sin(theta);
theta_dot = (v_vec(2)*cos(theta)-v_vec(1)*sin(theta))/r;

clear x;

N=25;

dt = 1.0/N;

for i=1:2*N
    x(i,1:2) = [cos(theta)*r,sin(theta)*r];
    [theta,theta_dot,r,r_dot] = integrate_rk4(theta,theta_dot,r,r_dot,dt);
end

clf reset

plot2d(x(:,1),x(:,2),frameflag=3,rect=[-12,-12,12,12],style=-1);


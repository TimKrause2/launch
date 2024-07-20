function ydot = ode_func(t,y)
    ydot = [y(3); y(4); y(1)*y(4)^2-1/y(1)^2; -2*y(3)*y(4)/y(1)];
endfunction

Npoints = 100;

function plot_y(yr,N)
    clear x;
    for i=1:N
        x(i,1:2) = [cos(yr(2,i))*yr(1,i),sin(yr(2,i))*yr(1,i)];
    end
    plot2d(x(:,1),x(:,2),frameflag=4);
endfunction

function T = Torbit(y0)
    vr = y0(3);
    vt = y0(1)*y0(4);
    v2 = vr^2 + vt^2;
    E = v2/2.0 - 1.0/y0(1);
    a = -1/(2*E);
    T = 2*%pi*sqrt(a^3);
endfunction

y0 = [1;0;0;1];

tspan = linspace(0,2*%pi,Npoints);

yr = ode(y0,0.0,tspan,ode_func);

plot_y(yr,Npoints);


y0(4) = 1.1;

T = Torbit(y0);

tspan = linspace(0,T,Npoints);

yr = ode(y0,0.0,tspan,ode_func);

plot_y(yr,Npoints);


y0(4) = 0.8;

T = Torbit(y0);

tspan = linspace(0,T,Npoints);

yr = ode(y0,0.0,tspan,ode_func);

plot_y(yr,Npoints);




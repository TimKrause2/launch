function [r] = r_ellip(theta,a,e)
    r = (a*(1-e^2))/(1+e*cos(theta));
endfunction

function [rc] = r_ellip_c(theta,a,e)
    rc = [cos(theta)*r_ellip(theta,a,e),sin(theta)*r_ellip(theta,a,e)];
endfunction

function [drc] = dr_ellip_c(theta,a,e)
    den = 1 + e*cos(theta);
    A = a*(1-e^2);
    drc = [-A*sin(theta)/den^2,(A*(e+cos(theta)))/den^2];
//    drc = [A*e*cos(theta)*sin(theta)/den^2-A*sin(theta)/den,A*e*sin(theta)^2/den^2+A*cos(theta)/den]
endfunction

function [n] = normal(A)
    mag = sqrt(sum(A.*A));
    n = A/mag;
endfunction

a = 1;
e = 0.5;

Npoints = 33;

theta = linspace(0.0,1.0,Npoints)*2*%pi;

i=1;
clear rc_ellipse;
clear rc_segs;
for t = theta do
    rc_ellipse(i,1:2) = r_ellip_c(t,a,e);
    drc = dr_ellip_c(t,a,e);
    drc = normal(drc)*0.5;
    rc_segs((i-1)*2+1,1:2) = rc_ellipse(i,1:2);
    rc_segs((i-1)*2+2,1:2) = rc_ellipse(i,1:2)+drc;
    i = i + 1;
end

plot2d(rc_ellipse(:,1),rc_ellipse(:,2),frameflag=4);
xsegs(rc_segs(:,1),rc_segs(:,2));



function [theta1,theta2]=angular_limits(t1,t2)
    if t1>=0 then
        if t2>=0 then
            if (t2-t1)>0 then
                theta1 = t1;
                theta2 = t2;
            else
                theta1 = t1;
                theta2 = t2+2*%pi;
            end
        else
            theta1 = t1;
            theta2 = t2+2*%pi;
        end
    else
        if t2>=0 then
            theta1 = t1;
            theta2 = t2;
        else
            if (t2-t1)>0 then
                theta1 = t1;
                theta2 = t2;
            else
                theta1 = t1;
                theta2 = t2+2*%pi;
            end
        end
    end
endfunction











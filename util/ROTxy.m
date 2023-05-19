function [Z] = ROTxy(Z, th0, xy0)
    %help function to rotate FRS in xy-plane
    RR = [cos(th0), -sin(th0); sin(th0) cos(th0)];
    tempc = center(Z);
    tempg = generators(Z);
    tempg(1:2,:) = RR*tempg(1:2,:);
    tempc(3) = tempc(3)+th0;
    
    if nargin == 3
        tempc(1:2) = RR*(tempc(1:2) - xy0) + xy0; 
    end
    
%     tempg(4:5,:) = RR*tempg(4:5,:);
    Z = zonotope([tempc,tempg]);
end






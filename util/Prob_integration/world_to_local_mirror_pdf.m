function pdf_out = world_to_local_mirror_pdf(robot_pose, pdf_world, mirror_flag)
% pdf_world: [mu_x_1, mu_y_1, sigma_xx_1, sigma_xy_1, sigma_yy_1, additional_info; 
%             mu_x_2, mu_y_2, sigma_xx_2, sigma_xy_2, sigma_yy_2, additional_info; ...]

    % extract position and heading from input
    x = robot_pose(1,1) ;
    y = robot_pose(2,1) ;
    h = robot_pose(3,1) ;
    
    % prep
    pdf_out = pdf_world ;
    pdf_out(:,1:2) = pdf_world(:,1:2) - [x y];
    ROT = [ cos(h), sin(h) ;
         -sin(h), cos(h) ] ;
    pdf_out(:,1:2) = pdf_out(:,1:2) * ROT';

    % rotate to local frame
    pdf_out = pdf_out';
    pdf_out(3:5,:) = [      cos(h)^2,     2*cos(h)*sin(h),      sin(h)^2;
                      -cos(h)*sin(h), cos(h)^2 - sin(h)^2, cos(h)*sin(h);
                            sin(h)^2,    -2*cos(h)*sin(h),      cos(h)^2] * pdf_out(3:5,:);
    pdf_out = pdf_out';

    if mirror_flag
        pdf_out(:,2) = -pdf_out(:,2);
        pdf_out(:,4) = -pdf_out(:,4);
    end


 end

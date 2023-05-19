function [vehicle] = plot_vehicle(state_x, state_y, state_h, body_color, window_color, opacity)
%plot_vehicle: Plot vehicles for simulation and plots
%   plot_vehicle(state_x, state_y, state_h, color, opacity)

% Vehicle dimensions
car_width = 2.2;
car_length = 4.8;

% Window dimensions
window_width = 1.7;
roof_length = 1.6;
window_offset = [-0.2, 0];
roof_offset = [-0.4;0];

% body_color_blue = [129,158,228]/255; % blue
% body_color_red = [229,116,108]/255; % red
% window_color = [213,226,242]/255; % light blue


radius = [0.3, 0.3, 0.3, 0.3];
[outer_body_x, outer_body_y] = generate_rounded_edges(radius, state_x, state_y, car_length, car_width);
radius = [0.5, 0.3, 0.3, 0.5];
[inner_body_x, inner_body_y] = generate_rounded_edges(radius, state_x, state_y, 3, window_width);
inner_body_x = inner_body_x + window_offset(1);
inner_body_y = inner_body_y + window_offset(2);
window_left_x = inner_body_x(inner_body_x<0);
window_left_y = inner_body_y(inner_body_x<0);
window_right_x = inner_body_x(inner_body_x>0);
window_right_y = inner_body_y(inner_body_x>0);
% Rotate the vehicle
outer_body_coords = [outer_body_x; outer_body_y];
inner_body_coords = [inner_body_x; inner_body_y];
window_left_coords = [[window_left_x; window_left_y], [roof_offset(1)-roof_length/2, roof_offset(1)-roof_length/2; -window_width/2, window_width/2]];
window_right_top = [window_right_x(window_right_y>0); window_right_y(window_right_y>0)];
window_right_bottom = [window_right_x(window_right_y<0); window_right_y(window_right_y<0)];
window_right_coords = [[roof_offset(1)+roof_length/2; -window_width/2], window_right_bottom, window_right_top, [roof_offset(1)+roof_length/2; window_width/2]];
roof_coords = [roof_offset(1)+roof_length/2, roof_offset(1)-roof_length/2, roof_offset(1)-roof_length/2, roof_offset(1)+roof_length/2; window_width/2, window_width/2, -window_width/2, -window_width/2];
rot = [cos(state_h), -sin(state_h);
       sin(state_h),  cos(state_h)];
outer_body_rot = rot*outer_body_coords;
inner_body_rot = rot*inner_body_coords;
window_left_rot = rot*window_left_coords;
window_right_rot = rot*window_right_coords;
roof_rot = rot*roof_coords;

outer_body = outer_body_rot + [state_x; state_y];
inner_body = inner_body_rot + [state_x; state_y];
window_left = window_left_rot + [state_x; state_y];
window_right = window_right_rot + [state_x; state_y];
roof = roof_rot + [state_x; state_y];

% Plot the vehicle
body_polyshape = polyshape({outer_body(1,:), inner_body(1,:)}, {outer_body(2,:), inner_body(2,:)});
roof_polyshape = polyshape(roof(1,:), roof(2,:));
window_left_polyshape = polyshape(window_left(1,:), window_left(2,:));
window_right_polyshape = polyshape(window_right(1,:), window_right(2,:));
outline_polyshape =  polyshape(outer_body(1,:), outer_body(2,:));

body = plot(body_polyshape, 'FaceColor', body_color, 'EdgeColor','none');
hold on;
roof = plot(roof_polyshape, 'FaceColor', body_color, 'EdgeColor','none');
window_left = plot(window_left_polyshape, 'FaceColor', window_color, 'EdgeColor','none');
window_right = plot(window_right_polyshape, 'FaceColor', window_color, 'EdgeColor','none');
outline = plot(outline_polyshape, 'FaceColor','none', 'EdgeColor','k', 'LineWidth', 1); % Change the thickness of vehicle
vehicle = [body, roof, window_left, window_right, outline];

for i = 1:length(vehicle)
    set(vehicle(i),'FaceAlpha',opacity);
    set(vehicle(i),'EdgeAlpha',opacity);
end
end

function [arc_coord_x, arc_coord_y] = generate_arc(radius, corner_num, corner_coordinate)
    % Generates coordinates for the corners of rounded rectangles where the
    % corner_num is 2   1
    %               3   4
    % and the corner_coordinate is coordinate of the vertex
    % Choose the quarter of the circle for the arc
    theta = 0:0.1:pi/2;
    theta = theta + pi*(corner_num-1)/2;
    % Generate the arc
    arc_coord_x = radius*cos(theta);
    arc_coord_y = radius*sin(theta);
    % Translate the arc to the correct vertex
    if corner_num==1 || corner_num==4
        arc_coord_x = arc_coord_x + corner_coordinate(1) - radius;
    else
        arc_coord_x = arc_coord_x + corner_coordinate(1) + radius;
    end
    if corner_num==1 || corner_num==2
        arc_coord_y = arc_coord_y + corner_coordinate(2) - radius;
    else
        arc_coord_y = arc_coord_y + corner_coordinate(2) + radius;
    end
end

function [arc_coord_x, arc_coord_y] = generate_rounded_edges(radius, state_x, state_y, length, width)
    % Generates the four rounded edges for a rectangle
    vertex1 = [length/2; width/2];
    vertex2 = [-length/2; width/2];
    vertex3 = [-length/2; -width/2];
    vertex4 = [length/2; -width/2];

    vertices = [vertex1, vertex2, vertex3, vertex4];
    arc_coord_x = [];
    arc_coord_y = [];

    for i = 1:1:4
        [temp_x, temp_y] = generate_arc(radius(i), i, vertices(:,i));
        arc_coord_x = [arc_coord_x, temp_x];
        arc_coord_y = [arc_coord_y, temp_y];
    end
end

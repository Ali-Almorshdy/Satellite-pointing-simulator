function [pl]=plot_3d_line2(p1, p2)
% Plots a 3D line given two points in 3D space.

% Extract the x, y, and z coordinates of the two points
x = [p1(1), p2(1)];
y = [p1(2), p2(2)];
z = [p1(3), p2(3)];

% Plot the line
pl=plot3(x, y, z, '-r','LineWidth',2);

end

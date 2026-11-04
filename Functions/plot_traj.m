function [p,IC] = plot_traj(r_M,color,dispName,ICname)
% plots 3d trajectories without me typing it out each time
% r_M: Position vector matrix [3xn] where row1 = x, row2 = y, row3 = z
% dispName, name of the plot, good for using legends in long codes that
% include surfaces

switch nargin
    case 5

    case 4
        % Make sure a figure is opened before use, use hold on if Earth Moon plot
        % is being used,
        
        % Initial Condition
        IC = plot3(r_M(1,1),r_M(2,1),r_M(3,1),'^','Color',color,'DisplayName',ICname);
        p = plot3(r_M(1,:),r_M(2,:),r_M(3,:),'Color',color,'DisplayName',dispName);
    case 3
        p = plot3(r_M(1,:),r_M(2,:),r_M(3,:),'Color',color,'DisplayName',dispName);
    case 2
        p = plot3(r_M(1,:),r_M(2,:),r_M(3,:),'Color',color);
    case 1
        p = plot3(r_M(1,:),r_M(2,:),r_M(3,:));
end
end

function [Moon] = plot_M(mu,LU,target)
%% Plots L1 and Moon for CR3BP trajectory in Rotating Frame
% with the Earth, m1 to the left of the origin at (-mu,0)
% and the Moon, or m2, is at (1 - mu, 0)
%
%                L4
%
% -L3----m1--+-------L1--m2--L2-
%
%                L5
%
% From Shane Ross
%
% Plotting
if nargin < 3
figure
end

switch nargin
    case 0
        % Values taken from https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
        mu = 1.215058560962404E-2;
        LU = 389703; % km
        target= gca;
    case 1
        LU = 389703; % km
        target = gca;
    case 2
        target = gca;
end

earthRadius = 6378.137; % km
moonRadius = 1737.1; % km


% Earth
% [xE,yE,zE] = sphere(50); 
% xE = xE*earthRadius/EMdist - mu;
% yE = yE*earthRadius/EMdist;
% zE = zE*earthRadius/EMdist;

% Moon
[xM,yM,zM] = sphere(50);
xM = xM*moonRadius/LU + (1 - mu);
yM = yM*moonRadius/LU;
zM = zM*moonRadius/LU;


% surface(xE,yE,zE)
hold on 
Moon = surface(target,xM,yM,zM,FaceColor=128/255*[1 1 1],EdgeColor=128/255*[1 1 1],DisplayName="Moon");
hold off


% Axis Labels
grid on
xlabel("${{\hat b}}_x \, (DU)$",'Interpreter','latex')
ylabel("${{\hat b}}_y \, (DU)$",'Interpreter','latex')
zlabel("${{\hat b}}_z \, (DU)$",'Interpreter','latex')

axis equal
% Lagrange Points
L1 = [0.83691513 0 0]';
L2 = [1.15568217	0.00000000	0.00000000]';
L3 = [-1.00506265	0.00000000	0.00000000]';
L4 = [0.48784941	0.86602540	0.00000000]';
L5 = [0.48784941	-0.86602540	0.00000000]';

hold on
% [L1_p] = plot3(L1(1),L1(2),L1(3),'diamondg','MarkerFaceColor','g','DisplayName',"L1");
%plot3(L2(1),L2(2),L2(3),'diamondg','LineWidth',1);
%plot3(L3(1),L3(2),L3(3),'diamondg','LineWidth',1);
%plot3(L4(1),L4(2),L4(3),'diamondg','LineWidth',1);
%plot3(L5(1),L5(2),L5(3),'diamondg','LineWidth',1);
hold off

end
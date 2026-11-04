%% Optimal Control Project Testing - Transfer from L2 to L1 Halo Orbits
% Getting minEnergy trajectory Shooting Method
clc;clear;close all

%% Add paths
addpath("Functions\")


%% Plot Settings
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot,'defaultAxesZGrid','on')
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex')
set(groot, 'defaultAxesFontSize',20)
set(groot,"defaultLineLineWidth",1)

cc = orderedcolors("gem12");

%% General Earth-Moon Info
% Taken from https://ssd.jpl.nasa.gov/tools/periodic_orbits.html
% Non-dimensional position and velocity (LU) and (LU/TU)
% LU = 389703; % Km
% TU = 382981; % s

TU = 3.75162997e5; % s Values used in Hecht's Paper
LU = 3.84400000e5; % km
VU = LU/TU; % km/s
VUmm = VU*1e6; % mm/s
VUcm = VU*1e5; % cm/s
day = 24*3600; % s
% mu = 1.215058560962404E-2; % Non-Dimensional Mass parameter (Earth-Moon)
mu = 1.21506038e-2;
Omega = 2*pi/TU;

massM = 73476730924573500000000;
massE = 5972190000000000000000000;
MU = massM + massE;

%% Initial and Final Orbits

ri = [1.1599795702248494 0.009720428035815552 -0.12401864915284157]';
vi = [0.008477705130550553 -0.20786307954141953 -0.010841912833115475]';

xi = [ri; vi]; x0 = xi;

T1 = 14.2*day/TU;

rf = [0.8484736688482315 0.00506488863463682 0.17343680487577373]';
vf = [0.005241131023638693 0.26343491250951045 -0.008541420325316247]';

xf = [rf;vf];

T2 = 11.2*day/TU;

L1_orbit = prop_cr3bp(xi,mu,[0 T1]);
L2_orbit = prop_cr3bp(xf,mu,[0 T2]);

plot_M(mu,LU);
hold on
plot_traj(L1_orbit,"k","$L_1$ Orbit");
plot_traj(L2_orbit,"k","$L_2$ Orbit");

plot3(ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
plot3(rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")
legend()


%% Simulation Values

Isp = 2000/TU;
g0 = 9.81e-3 *TU^2/LU;
mi = 2e3; % kg
Tmax = 1.5e-3 * TU^2/LU/mi; % Normalized to initial Mass of Satellite


c = Isp*g0;

tf = 12.7*day/TU;

m = 3e3;
t = linspace(0,tf,m);

% Optimal Initial Costates from Hecht paper as
% initial guess 

trajType = 'alpha';

switch trajType
    case 'alpha'

        % Alpha Trajectory
        pr0 = [0.12603 -0.07665 -0.05635]';
        pv0 = [0.03999 -0.00518 -0.06410]';
        pm0 = 0.02236;

        gam = 1e-9;

    case 'beta'
        % Beta Trajectory
        pr0 = [-0.01486, 0.01215, -0.07936]';
        pv0 =  [0.01015, 0.04457, 0.01256]';
        pm0 = 0.07632;

        gam = 1e-8;

    case 'gamma'
        % Gamma trajectory
        pr0 = [-0.02195, 0.00659, 0.07490]';
        pv0 = [-0.04314, 0.03615, 0.03842]';
        pm0 = 0.03489;

        gam = 1e-8;

end

m0 = 1;
p0 = [pr0;pv0;pm0];

N = 25;
j = 1:25;

% epj = (j.^2 - 1)./(N^2 - 1);
epj = linspace(0,1,N);
P0s = zeros(7,N);
Us = zeros(N,m);
Ss = zeros(N,m);
Xs = zeros(6,m,N);


for jj = 1:N

eppy = epj(jj);




maxIter = 500;

for k = 1:maxIter

    y0 = [x0;m0; p0];
    [y,t,te,ye] = propState(y0,t,mu,Tmax,c,eppy);


    x = y(1:6,:);
    p = y(8:13,:);
    m = y(7,:);
    pm = y(14,:);

    % Differences
    dxf = xf - x(:,end);
    dpmf = 0 - pm(:,end);

    

    % Termination criteria
    if k> 1
        diffx = abs(x -xm1);
        diffp = abs(p - pm1);

        normy = sum(max(diffx,[],2) + max(diffp,[],2));

        normy2 = norm(dxf);

        if normy < gam && normy2 < gam
            break
        end
    end

    % Store previous trajectory
    xm1 = x;
    pm1 = p;

    
    %% Now Develop correction using STM
    ddt = 1e-10; % Buffer time

    if ~isempty(te)
    

    % Iteration 1
    ym = propState(y0,[0 te(1)-ddt],mu,Tmax,c,eppy);
    yp = propState(y0,[0 te(1)+ddt],mu,Tmax,c,eppy);
    ym = ym(:,end);
    yp = yp(:,end);

    PHItemp = propPHI(y0,[0 te(1)-ddt],mu,Tmax,c,eppy);
    PSItemp = PSI(ym,yp,mu,Tmax,c,eppy);

    mat = PSItemp*PHItemp;

    yp_last = yp;

    for i = 2:length(te)
        ym = propState(y0,[0 te(i)-ddt],mu,Tmax,c,eppy);
        yp = propState(y0,[0 te(i)+ddt],mu,Tmax,c,eppy);
        ym = ym(:,end);
        yp = yp(:,end);


        PHItemp = propPHI(yp_last,[te(i-1)+ddt te(i)-ddt],mu,Tmax,c,eppy);
        Psitemp = PSI(ym,yp,mu,Tmax,c,eppy);

        mat = Psitemp*PHItemp*mat;

        yp_last = yp;


    end

    Phitemp = propPHI(yp_last,[te(end)+ddt tf],mu,Tmax,c,eppy);

    PHI = Phitemp*mat;

    else
        PHI = propPHI(y0,[0 tf],mu,Tmax,c,eppy);
    end



    phi13 = PHI(1:6,8:13);
    phi14 = PHI(1:6,14);

    phi43 = PHI(14,8:13);
    phi44 = PHI(14,14);

    phi = [phi13 phi14; phi43 phi44];

    % Correction
    dp0 = (phi\[dxf;dpmf]);

    p0 = p0 + dp0;


end
% Get u and Switching curve
[u,S] = getu(y,c,eppy);

P0s(:,jj) = p0;
Us(jj,:) = u;
Ss(jj,:) = S;
Xs(:,:,jj) = x;

% Final Difference in End Conditions
DXF = dxf .* [LU*1e6*ones(3,1); VU*1e6*ones(3,1)];

fprintf("For \x03F5 = %.4f:\n",eppy)
fprintf("Iterations: %d\n\n",k)
fprintf("Final Difference in End Conditions:\n")
fprintf("   \x0394x  = %.4e mm\n",DXF(1))
fprintf("   \x0394y  = %.4e mm\n",DXF(2))
fprintf("   \x0394z  = %.4e mm\n",DXF(3))
fprintf("   \x0394Vx = %.4e mm/s\n",DXF(4))
fprintf("   \x0394Vy = %.4e mm/s\n",DXF(5))
fprintf("   \x0394Vz = %.4e mm/s\n",DXF(6))

fprintf("Distance Norm: %.4e mm \n",norm(DXF(1:3)))

delta_m = (1 - m(end))*mi;

fprintf('\x0394m  = %.4f kg\n\n\n',delta_m)




end



fig3 = figure();
ax3 = axes(fig3);

plot_M(mu,LU,ax3);
hold(ax3,"on")
plot3(ax3,L1_orbit(1,:),L1_orbit(2,:),L1_orbit(3,:),'k')
plot3(ax3,L2_orbit(1,:),L2_orbit(2,:),L2_orbit(3,:),'k')

plot3(ax3,ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
plot3(ax3,rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")

cc = jet(25);

for i = 1:N
plot3(ax3,Xs(1,:,i),Xs(2,:,i),Xs(3,:,i),'Color',cc(i,:),'DisplayName',"Transfer Trajectory")
end
axis equal
colormap(cc)
cb = colorbar;
cb.Label.String = "\epsilon";
clim([0 1])

ind = [1 7:6:N];

figure
colororder(cc);
subplot(2,1,1)
plot(t,Us)
ylabel("Throttle u(t)")
% legend(append("$\epsilon =$ ",num2str(epj(ind(1)))),append("$\epsilon =$ ",num2str(epj(ind(2)))),append("$\epsilon =$ ",num2str(epj(ind(3)))),...
    % append("$\epsilon =$ ",num2str(epj(ind(4)))),append("$\epsilon =$ ",num2str(epj(ind(5)))),'Location','best')
subplot(2,1,2)
plot(t,Ss)
ylabel("Switching Curve")






%% Functions
function [sdot] = cr3bp_lowThrust_eqns(t,s,mu,u,alph,Tmax,c)
% s = [r' v' m]'

x = s(1:6);
m = s(7);

xdot = CR3BP_EOM(t,x,mu);

xdot(4:6) = xdot(4:6) + u*Tmax*alph/m;

mdot = -u*Tmax/c;

sdot = [xdot;mdot];



end

function [ydot,u,S] = stateEq(t,y,mu,Tmax,c,eppy)
% y: States and Costates - see Hecht's Paper

x = y(1:6);
m = y(7);
pr = y(8:10);
pv = y(11:13);
pm = y(14);

alph = -pv./norm(pv);

S = -c/m *norm(pv) - pm + 1; % Switching curve


if S > eppy
    u = 0;
elseif S < -eppy
    u = 1;
else
    u = (eppy - S)/(2*eppy);
end

s = [x; m];

sdot = cr3bp_lowThrust_eqns(t,s,mu,u,alph,Tmax,c);

G = getG(x,mu);
H = 2*[0 1 0;-1 0 0; 0 0 0];

ydot = zeros(14,1);

ydot(1:7) = sdot;

ydot(8:10) = -G'*pv;

ydot(11:13) = -pr - H'*pv;

ydot(14) = -norm(pv)*u*Tmax/m^2;


end

function [G] = getG(X,mu)

mu1 = 1 - mu;
mu2 = mu;
x = X(1);
y = X(2);
z = X(3);

r1 = sqrt((x + mu2)^2 + y^2 + z^2);
r2 = sqrt((x - mu1)^2 + y^2 + z^2);

% From differentiating the potential energy term U with respect to position
Uxx = 1 - mu1/r1^3 - mu2/r2^3 + 3*mu1*(x + mu2)^2/r1^5 + 3*mu2*(x - mu1)^2/r2^5;

Uxy = 3*mu1*(x + mu2)*y/r1^5 + 3*mu2*(x-mu1)*y/r2^5;

Uxz = 3*mu1*(x + mu2)*z/r1^5 + 3*mu2*(x-mu1)*z/r2^5;

Uyy = 1 - mu1/r1^3 - mu2/r2^3 + 3*mu1*y^2/r1^5 + 3*mu2*y^2/r2^5;

Uyz = 3*mu1*y*z/r1^5 + 3*mu2*y*z/r2^5;

Uzz = -mu1/r1^3 - mu2/r2^3 + 3*mu1*z^2/r1^5 + 3*mu2*z^2/r2^5;

% Symmetrical nature of the problem
Uyx = Uxy;
Uzx = Uxz;
Uzy = Uyz;

G = [Uxx Uxy Uxz;
    Uyx Uyy Uyz;
    Uzx Uzy Uzz];

end


function [dGpvdr] = dGpvdr(X,pv,mu)

x = X(1);
y = X(2);
z = X(3);

r1 = sqrt((x+mu)^2 + y^2 + z^2);
r2 = sqrt( (x + mu-1)^2 + y^2 + z^2);

% G11
dG11dx = 3*(1 - mu)*(x+ mu)/r1^5 + ( 6*(1-mu)*(x+mu)/r1^5 - 15*(1-mu)*(x+mu)^3/r1^7 ) + 3*mu*(x+mu-1)/r2^5 + ( 6*mu*(x+mu-1)/r2^5 - 15*mu*(x+mu-1)^3/r2^7 );

dG11dy = 3*(1-mu)*y/r1^5 -15*(1-mu)*(x+mu)^2*y/r1^7 + 3*mu*y/r2^5 - 15*mu*(x+mu-1)^2*y/r2^7;

dG11dz = 3*(1-mu)*z/r1^5 -15*(1-mu)*(x+mu)^2*z/r1^7 + 3*mu*z/r2^5 - 15*mu*(x+mu-1)^2*z/r2^7;

% G22
dG22dx = 3*(1-mu)*(x+mu)/r1^5 - 15*(1-mu)*y^2*(x+mu)/r1^7 + 3*mu*(x+mu-1)/r2^5 - 15*mu*y^2*(x+mu-1)/r2^7;

dG22dy = 3*(1-mu)*y/r1^5 + ( 6*(1-mu)*y/r1^5 - 15*(1-mu)*y^3/r1^7 ) + 3*mu*y/r2^5 + ( 6*mu*y/r2^5 - 15*mu*y^3/r2^7 );

dG22dz = 3*(1-mu)*z/r1^5 - 15*(1-mu)*y^2*z/r1^7 + 3*mu*z/r2^5 - 15*mu*y^2*z/r2^7;

% G33

dG33dx = 3*(1-mu)*(x+mu)/r1^5 - 15*(1-mu)*z^2*(x+mu)/r1^7 + 3*mu*(x+mu-1)/r2^5 - 15*mu*z^2*(x+mu-1)/r2^7;

dG33dy = 3*(1-mu)*y/r1^5 - 15*(1-mu)*z^2*y/r1^7 + 3*mu*y/r2^5 - 15*mu*z^2*y/r2^7;

dG33dz = 3*(1-mu)*z/r1^5 + ( 6*(1-mu)*z/r1^5 - 15*(1-mu)*z^3/r1^7 ) + 3*mu*z/r2^5 + ( 6*mu*z/r2^5 - 15*mu*z^3/r2^7 );

% G12

dG12dx = ( 3*(1-mu)*y/r1^5 - 15*(1-mu)*y*(x+mu)^2/r1^7 ) + ( 3*(mu)*y/r2^5 - 15*(mu)*y*(x+mu-1)^2/r2^7 ); 

dG12dy = ( 3*(1-mu)*(x+mu)/r1^5 - 15*(1-mu)*(x+mu)*y^2/r1^7 ) +  ( 3*(mu)*(x+mu-1)/r2^5 - 15*(mu)*(x+mu-1)*y^2/r2^7 );

dG12dz = -15*(1-mu)*(x+mu)*y*z/r1^7 - 15*mu*(x+mu-1)*y*z/r2^7;

% G13

dG13dx = ( 3*(1-mu)*z/r1^5 - 15*(1-mu)*z*(x+mu)^2/r1^7 ) + ( 3*(mu)*z/r2^5 - 15*(mu)*z*(x+mu-1)^2/r2^7 );

dG13dy = dG12dz;

dG13dz = ( 3*(1-mu)*(x+mu)/r1^5 - 15*(1-mu)*(x+mu)*z^2/r1^7 ) +  ( 3*(mu)*(x+mu-1)/r2^5 - 15*(mu)*(x+mu-1)*z^2/r2^7 );

% G23

dG23dx = -15*(1-mu)*y*z*(x+mu)/r1^7 - 15*(mu)*y*z*(x+mu-1)/r2^7;

dG23dy = ( 3*(1-mu)*z/r1^5 - 15*(1-mu)*z*y^2/r1^7) + ( 3*mu*z/r2^5 - 15*mu*z*y^2/r2^7);

dG23dz = ( 3*(1-mu)*y/r1^5 - 15*(1-mu)*y*z^2/r1^7) + ( 3*mu*y/r2^5 - 15*mu*y*z^2/r2^7);

J11 = dG11dx*pv(1) + dG12dx*pv(2) + dG13dx*pv(3);
J12 = dG11dy*pv(1) + dG12dy*pv(2) + dG13dy*pv(3);
J13 = dG11dz*pv(1) + dG12dz*pv(2) + dG13dz*pv(3);

J21 = dG12dx*pv(1) + dG22dx*pv(2) + dG23dx*pv(3);
J22 = dG12dy*pv(1) + dG22dy*pv(2) + dG23dy*pv(3);
J23 = dG12dz*pv(1) + dG22dz*pv(2) + dG23dz*pv(3);

J31 = dG13dx*pv(1) + dG23dx*pv(2) + dG33dx*pv(3);
J32 = dG13dy*pv(1) + dG23dy*pv(2) + dG33dy*pv(3);
J33 = dG13dz*pv(1) + dG23dz*pv(2) + dG33dz*pv(3);


dGpvdr = [J11 J12 J13; J21 J22 J23; J31 J32 J33];


end



function [y,t,te,ye] = propState(y0,timeVec,mu,Tmax,c,eppy)

opt = odeset("AbsTol",1e-12,"RelTol",1e-12,"Events",@(t,y) eventS0(t,y,c,eppy));

[t,y,te,ye] = ode89(@(t,y) stateEq(t,y,mu,Tmax,c,eppy),timeVec,y0,opt);

y = y';
ye = ye';


end


function [value,isterminal,direction] = eventS0(t,y,c,eppy)

% x = y(1:6);
m = y(7);
% pr = y(8:10);
pv = y(11:13);
pm = y(14);

% alph = -pv./norm(pv);

S = -c/m *norm(pv) - pm + 1; % Switching curve

if eppy > 0
value = [S - eppy ;S + eppy];

isterminal = [0 0]';
direction = [0 0]';

else
    value = S;
    isterminal = 0;
    direction = 0;
end

end



function [PHI] = propPHI(y0,timeVec,mu,Tmax,c,eppy)


opt = odeset("AbsTol",1e-12,"RelTol",1e-12);

I14 = eye(14);

v0 = [y0; reshape(I14,[],1)];

[~,v] = ode89(@(t,y) PHIeqns(t,y,mu,Tmax,c,eppy),timeVec,v0,opt);
v = v';
PHI = reshape(v(15:end,end),[14 14]);




end

function [PSI] = PSI(ym,yp,mu,Tmax,c,eppy)

% ys is state at switch
% sw is order of switching, [0 1] or [1 0]

% y: States and Costates - see Hecht's Paper

% x = ys(1:6);
% m = ys(7);
% pr = ys(8:10);
% pv = ys(11:13);
% pm = ys(14);
% 
% alph = -pv./norm(pv);


mm = ym(7);
prm = ym(8:10);
pvm = ym(11:13);


H = 2*[0 1 0;-1 0 0; 0 0 0];




[ydotm] = stateEq(0,ym,mu,Tmax,c,eppy);
[ydotp] = stateEq(0,yp,mu,Tmax,c,eppy);



dSdy = [zeros(1,6) c/mm^2*norm(pvm) zeros(1,3) -c/mm*pvm'/norm(pvm) -1];

Sdot = c/mm*( prm + H'*pvm)'*pvm/norm(pvm);



PSI = eye(14) + (ydotp - ydotm)*(dSdy/Sdot);


end





function [u,S] = getu(y,c,eppy)

n = size(y,2);

u = nan(1,n);
S = nan(1,n);

m = y(7,:);
pv = y(11:13,:);
pm = y(14,:);

for k = 1:n

    S(k) = -c/m(k) *norm(pv(:,k)) - pm(k) + 1; % Switching curve

    if S(k) > eppy
        u(k) = 0;
    elseif S(k) < -eppy
        u(k) = 1;
    else
        u(k) = (eppy - S(k))/(2*eppy);
    end

end

end


function [cost] = qCost(s,dp0,x0,m0,p0,timeVec,mu,Tmax,c,eppy,xf)

p0 = p0 + s*dp0;

y0 = [x0;m0;p0];

y = propState(y0,timeVec,mu,Tmax,c,eppy);
x = y(1:6,:);
pm = y(14,:);

dxf = xf - x(:,end);
dpmf = 0-  pm(end);

cost = norm(dxf) + abs(dpmf); 

end

% function [vdot] = PHIeqns(t,v,mu,Tmax,c,eppy)
% 
% % v: state vector [14 + 14^2 x 1] column
% % v(1:14) is state/costate vector, remaining fill STM matrix
% 
% vdot = zeros(14+14^2,1);
% 
% y = v(1:14);
% 
% x = y(1:6);
% m = y(7);
% pr = y(8:10);
% pv = y(11:13);
% pm = y(14);
% 
% alph = -pv./norm(pv);
% 
% [ydot,u] = stateEq(t,y,mu,Tmax,c,eppy);
% 
% vdot(1:14) = ydot;
% 
% G = getG(x,mu);
% H = 2*[0 1 0;-1 0 0; 0 0 0];
% 
% 
% % Note these are valid only for u = 0 or u = 1 (which is what i've been
% % using since eppy = 0)
% 
% DyF11 = [zeros(3)       eye(3)            zeros(3,1);
%     G             H        pv/norm(pv)*u*Tmax/m^2;
%     zeros(1,3)     zeros(1,3)            0];
% 
% DyF12 = zeros(7); DyF12(4:6,4:6) = -u*Tmax/m*( eye(3)/norm(pv) - pv*pv'/norm(pv)^3);
% 
% DyF21 = zeros(7); DyF21(1:3,1:3) = -dGpvdr(x,pv,mu); DyF21(7,7) = 2*norm(pv)*u*Tmax/m^3;
% 
% DyF22 = -DyF11';
% 
% DyF = [DyF11 DyF12; DyF21 DyF22];
% 
% PHI = reshape(v(15:end),[14 14]);
% 
% PHIdot = DyF*PHI;
% 
% vdot(15:end) = reshape(PHIdot,[],1);
% 
% end

function [vdot] = PHIeqns(t,v,mu,Tmax,c,eppy)

% v: state vector [14 + 14^2 x 1] column
% v(1:14) is state/costate vector, remaining fill STM matrix

vdot = zeros(14+14^2,1);

y = v(1:14);

x = y(1:6);
m = y(7);
pr = y(8:10);
pv = y(11:13);
pm = y(14);

alph = -pv./norm(pv);

[ydot,u,S] = stateEq(t,y,mu,Tmax,c,eppy);

vdot(1:14) = ydot;

G = getG(x,mu);
H = 2*[0 1 0;-1 0 0; 0 0 0];

if (S > eppy) || (S < -eppy)
    % Note these are valid only for u = 0 or u = 1 (which is what i've been
    % using since eppy = 0)

    DyF11 = [zeros(3)       eye(3)            zeros(3,1);
        G             H        pv/norm(pv)*u*Tmax/m^2;
        zeros(1,3)     zeros(1,3)            0];

    DyF12 = zeros(7); DyF12(4:6,4:6) = -u*Tmax/m*( eye(3)/norm(pv) - pv*pv'/norm(pv)^3);

    DyF21 = zeros(7); DyF21(1:3,1:3) = -dGpvdr(x,pv,mu); DyF21(7,7) = 2*norm(pv)*u*Tmax/m^3;

    DyF22 = [zeros(3) -G' zeros(3,1); -eye(3) -H' zeros(3,1); zeros(1,3) -pv'/norm(pv)*u*Tmax/m^2 0];



else
    O1 = pv/norm(pv)*u*Tmax/m^2 + pv*c*Tmax/(2*eppy*m^3);
    O2 = -pv*pv'/norm(pv)^2 * c*Tmax/(2*eppy*m^2) - (u*Tmax/m)*( eye(3)/norm(pv) - pv*pv'/norm(pv)^3);
    O3 = -pv/norm(pv)*Tmax/(2*eppy*m);
    O4 = norm(pv)*Tmax/(2*eppy*m^2);
    O5 = -pv'/norm(pv)*Tmax/(2*eppy*m);
    O6 = -Tmax/(2*eppy*c);
    O7 = 2*norm(pv)*u*Tmax/m^3 + norm(pv)^2*c*Tmax/(2*eppy*m^4);
    O8 = -pv'/norm(pv)*u*Tmax/m^2 - pv'*c*Tmax/(2*eppy*m^3);
    O9 = -norm(pv)*Tmax/(2*eppy*m^2);

    DyF11 = [zeros(3) eye(3) zeros(3,1); G H O1; zeros(1,3) zeros(1,3) O4];

    DyF12 = [zeros(3) zeros(3) zeros(3,1); zeros(3) O2 O3; zeros(1,3) O5 O6];

    DyF21 = zeros(7); DyF21(1:3,1:3) = -dGpvdr(x,pv,mu); DyF21(7,7) = O7;

    DyF22 = [zeros(3) -G' zeros(3,1); -eye(3) -H' zeros(3,1); zeros(1,3) O8 O9];



end

DyF = [DyF11 DyF12; DyF21 DyF22];

PHI = reshape(v(15:end),[14 14]);

PHIdot = DyF*PHI;

vdot(15:end) = reshape(PHIdot,[],1);



end


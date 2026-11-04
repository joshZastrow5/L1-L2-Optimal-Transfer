%% Transfer 2, this time no given solution
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
set(groot,"defaultLineLineWidth",1.5)

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

x0 = [ri; vi];

rf = [0.8484736688482315 0.00506488863463682 0.17343680487577373]';
vf = [0.005241131023638693 0.26343491250951045 -0.008541420325316247]';

xf = [rf;vf];

L1_orbit = prop_cr3bp(xf,mu,[0 11.2*day/TU]);
L2_orbit = prop_cr3bp(x0,mu,[0 14.2*day/TU]);

plot_M(mu,LU);
hold on
plot_traj(L1_orbit,"k","$L_1$ Orbit");
plot_traj(L2_orbit,"k","$L_2$ Orbit");

plot3(ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
plot3(rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")
quiver3(ri(1),ri(2),ri(3),vi(1),vi(2),vi(3),0.5,"filled",'-r',"DisplayName","$v(t_0)$")
quiver3(rf(1),rf(2),rf(3),vf(1),vf(2),vf(3),0.5,"filled",'-b',"DisplayName","$v(t_f)$")
legend()
title("Original Problem")


T1 = 11.2*day/TU;
T2 = 14.2*day/TU;

% % New Problem
% x1 = prop_cr3bp(xf,mu,[0 T1/2]); % L1
% x2 = prop_cr3bp(x0,mu,[0 T2/2]); % L2
% 
% x0 = x1(:,end); ri = x0(1:3); vi = x0(4:6);
% xf = x2(:,end); rf = xf(1:3); vf = xf(4:6);
% 
% plot_M(mu,LU);
% hold on
% plot_traj(L1_orbit,"k","$L_1$ Orbit");
% plot_traj(L2_orbit,"k","$L_2$ Orbit");
% 
% plot3(ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
% plot3(rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")
% quiver3(ri(1),ri(2),ri(3),vi(1),vi(2),vi(3),0.5,"filled",'-r',"DisplayName","$v(t_0)$")
% quiver3(rf(1),rf(2),rf(3),vf(1),vf(2),vf(3),0.5,"filled",'-b',"DisplayName","$v(t_f)$")
% legend()
% title("New Problem")


%% Simulation Values

Isp = 2000/TU;
g0 = 9.81e-3 *TU^2/LU;
mi = 2e3; % kg
Tmax = 1.5e-3 * TU^2/LU/mi; % Normalized to initial Mass of Satellite


c = Isp*g0;

tf = 12.7*day/TU;

M = 3e3;
t = linspace(0,tf,M);
m0 = 1; 

s0 = [x0;m0];

eppy = 1; % Min Energy solution

%% Initial Guess

% u = .75*ones(1,M);
% alph = kron([1 0 0]',ones(1,M));

load("minEnergyControl.mat");

u = .3*cos(pi*t) + .5;

alph = [.1*sin(2*pi*t); cos(2*pi*t); sin(2*pi*t + pi/16)];

alph = alph./vecnorm(alph,2,1);



s = propState_rk4(s0,mu,u,alph,Tmax,c,t);
x = s(1:6,:);
m = s(7,:);


% Plot Trajectory
plot_M(mu,LU);
hold on
plot_traj(L1_orbit,"k","$L_1$ Orbit");
plot_traj(L2_orbit,"k","$L_2$ Orbit");

plot3(ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
plot3(rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")

plot_traj(s,'r','Transfer Trajectory');
legend()


figure
subplot(4,1,1)
plot(t,u)
ylabel("Throttle u(t)")

subplot(4,1,2)
plot(t,alph(1,:))
ylabel("$\alpha_1$")
subplot(4,1,3)
plot(t,alph(2,:))
ylabel("$\alpha_2$")
subplot(4,1,4)
plot(t,alph(3,:))
ylabel("$\alpha_3$")

dxf = xf - x(:,end);

DXF = dxf .* [LU*ones(3,1); VU*1e3*ones(3,1)];
fprintf("Initial Guess Results\n")
fprintf("For \x03F5 = %.4f:\n",eppy)
fprintf("Final Difference in End Conditions:\n")
fprintf("   \x0394x  = %.4e km\n",DXF(1))
fprintf("   \x0394y  = %.4e km\n",DXF(2))
fprintf("   \x0394z  = %.4e km\n",DXF(3))
fprintf("   \x0394Vx = %.4e m/s\n",DXF(4))
fprintf("   \x0394Vy = %.4e m/s\n",DXF(5))
fprintf("   \x0394Vz = %.4e m/s\n",DXF(6))

fprintf("Distance Norm: %.4e km \n",norm(DXF(1:3)))

delta_m = (1 - m(end))*mi;

fprintf('\x0394m  = %.4f kg\n\n\n',delta_m)

%% Steepest Descent 
gam1 = 1e-4;
gam2 = 1e-5;

maxIter = 500;
J = nan(1,maxIter);

fig = figure(100);
ax3 = axes(fig);

n = 7; % States
q = 6; % terminal constraints
w = 4; % control variables

for i = 1:maxIter

    
    pJf = zeros(n,1);
    pif = [eye(q); zeros(1,q)];
    pif = reshape(pif,[],1);
    pf = [pJf; pif];

    % Propagate co-states backwards
    p = propCoState(pf,flip(s,2),mu,flip(u),flip(alph,2),Tmax,flip(t));
    p = flip(p,2);

    lamJ = p(1:n,:);
    lami = reshape(p(n+1:end,:),n,q,[]);

    Lamphi = zeros(1,w,M);
    Lampsi = zeros(q,w,M);
    
    dg = zeros(q,M);
    dQ = zeros(q,q,M);
    for j = 1:M
        dfdu = zeros(n,w);
        dfdu(4:6,:) = [Tmax/m(j)*alph(:,j) u(j)*Tmax/m(j)*eye(3)];
        dfdu(7,1) = -Tmax/c;

        Lamphi(:,:,j) = [Tmax/c*(1 - eppy + 2*eppy*u(j)) zeros(1,3)] + lamJ(:,j)'*dfdu;
        Lampsi(:,:,j) = lami(:,:,j)'*dfdu;

        dg(:,j) = Lampsi(:,:,j)*Lamphi(:,:,j)';
        dQ(:,:,j) = Lampsi(:,:,j)*Lampsi(:,:,j)';


    end

    g = trapz(t,dg,2);
    Q = trapz(t,dQ,3);

    dPhie = -(xf -x(:,end));

    % Bound K value
    fk = @(K) phi(K,u,alph,Q,g,Lamphi,Lampsi,dPhie,s0,xf,mu,c,Tmax,t);
    % [kL,kU] = fibBound(.1,0.001,fk)
    % 
    % if kL <0
    %     kL = 0;
    % end

    % K = fminbnd(fk,0,10)

    K = .1;

    nu = -Q\( g - 1/K * dPhie);

    du = zeros(4,M);
    for k = 1:M
        du(:,k) = -K*( Lamphi(:,:,k) + nu'*Lampsi(:,:,k) )';
    end

    % du = [1 0; 0 .05]* du;
    % du = du./vecnorm(du,2,1);


    U = [u;alph];
    U = U + du;
    u = U(1,:);
    alph = U(2:4,:);

    % Corrections
    u(u>1) = 1;
    u(u<0) = 0;
    alph = alph./vecnorm(alph,2,1);

    s = propState_rk4(s0,mu,u,alph,Tmax,c,t);
    x = s(1:6,:);
    m = s(7,:);

    cla(ax3)
    plot_M(mu,LU,ax3);
    hold(ax3,"on")
    plot3(ax3,L1_orbit(1,:),L1_orbit(2,:),L1_orbit(3,:),'k')
    plot3(ax3,L2_orbit(1,:),L2_orbit(2,:),L2_orbit(3,:),'k')

    plot3(ax3,ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
    plot3(ax3,rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")
    plot3(ax3,ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
    plot3(ax3,rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")

    plot3(ax3,x(1,:),x(2,:),x(3,:),'-r','DisplayName',"Transfer Trajectory")
    title(append("Iteration: ",num2str(i)))

    drawnow

  
    % Termination criteria

    J(i) = nu'*dPhie + Tmax/c*trapz(t,u.^2,2);

    
    DXF = (xf - x(:,end)) .* [LU*ones(1,3) VU*1e3*ones(1,3)]';
    DRF = norm(DXF(1:3));
    DVF = norm(DXF(4:6));


    if i>1
        hnorm = norm(xf - x(:,end));
        jnorm = abs(J(i) - J(i-1));

        if ( hnorm < gam1) || (jnorm < gam2)
            break
        end

    end


end
J = J(~isnan(J));


plot_M(mu,LU);
hold on
plot_traj(L1_orbit,"k","$L_1$ Orbit");
plot_traj(L2_orbit,"k","$L_2$ Orbit");

plot3(ri(1),ri(2),ri(3),'sk',"DisplayName","Starting Position")
plot3(rf(1),rf(2),rf(3),'^k',"DisplayName","Final Position")

plot_traj(s,'r','Transfer Trajectory');



% Plot Control Profile 
figure
subplot(2,1,1)
plot(t,u)
ylabel("Throttle u(t)")

% Plot Mass
subplot(2,1,2)
plot(t,m*mi)
xlabel("Time (days)")
ylabel("Mass (kg)")


% Cost
figure
plot(1:i,J,'-o')
xlabel("Iteration")
ylabel("Cost")


fprintf("Iterations: %d\n",i)
fprintf("Final Difference in End Conditions:\n")
fprintf("   \x0394x  = %.4e km\n",DXF(1))
fprintf("   \x0394y  = %.4e km\n",DXF(2))
fprintf("   \x0394z  = %.4e km\n",DXF(3))
fprintf("   \x0394Vx = %.4e m/s\n",DXF(4))
fprintf("   \x0394Vy = %.4e m/s\n",DXF(5))
fprintf("   \x0394Vz = %.4e m/s\n",DXF(6))

fprintf("Distance Norm: %.4e km \n",norm(DXF(1:3)))

delta_m = (1 - m(end))*mi;

fprintf('\x0394m  = %.4f kg\n\n\n',delta_m)


% Get initial Co-state from Steepest Descent
p = propCoState(pf,flip(s,2),mu,flip(u),flip(alph,2),Tmax,flip(t));
p = flip(p,2);

lamJ = p(1:n,:);
lami = reshape(p(n+1:end,:),n,q,[]);

p0 = lamJ(:,1)' + nu'*lami(:,:,1)';

pr0 = p0(1:3);
pv0 = p0(4:6);
pm0 = p0(7);

fprintf("pr0 = [%.15f %.15f %.15f]';\n",pr0(1),pr0(2),pr0(3))
fprintf("pv0 = [%.15f %.15f %.15f]';\n",pv0(1),pv0(2),pv0(3))
fprintf("pm0 = %.15f;\n\n",pm0)


%% FUNCTIONS

function [J] = phi(K,u,alph,Q,g,Lamphi,Lampsi,dPhie,s0,xf,mu,c,Tmax,timeVec)

nu = -Q\( g - 1/K * dPhie );

M = length(timeVec);

du = zeros(4,M);
for i = 1:M
    du(:,i) = -K*( Lamphi(:,:,i) + nu'*Lampsi(:,:,i) )';
end

% du = du./vecnorm(du,2,1);
% du = [1 0; 0 .05]* du;
U =[u;alph];

U = U + du;
u = U(1,:);
alph = U(2:4,:);

u(u>1) = 1;
u(u<0) = 0;
alph = alph./vecnorm(alph,2,1);


s =propState_rk4(s0,mu,u,alph,Tmax,c,timeVec);


c1 = u.^2;

J = nu'*(xf -s(1:6,end)) + Tmax/c*trapz(timeVec,c1,2);

end

function [s] = propState_rk4(s0,mu,u,alph,Tmax,c,timeVec)

n = length(timeVec);
dt = timeVec(2) - timeVec(1);
s = zeros(length(s0),n);
s(:,1) = s0;
t = 0;
% I guess use RK4 
for k = 1:n-1
    s_k = s(:,k);
    % Runge-Kutta Steps 
    k1 = cr3bp_lowThrust_eqns(t,s(:,k),mu,u(k),alph(:,k),Tmax,c);
    k2 = cr3bp_lowThrust_eqns(t+dt/2,s(:,k)+k1*dt/2,mu,u(k),alph(:,k),Tmax,c);
    k3 = cr3bp_lowThrust_eqns(t+dt/2,s(:,k)+k2*dt/2,mu,u(k),alph(:,k),Tmax,c);
    k4 = cr3bp_lowThrust_eqns(t+dt,s(:,k)+k3*dt,mu,u(k),alph(:,k),Tmax,c);
    % Calculate next position and velocity
    s_kp1 = s_k + (dt/6)*(k1 +2*k2 +2*k3 + k4);
    % Store it in the next index for next loop
    s(:,k+1) = s_kp1;
end


end

function [p] = propCoState(p0,x,mu,u,alph,Tmax,timeVec)
n = length(timeVec);
dt = timeVec(2) - timeVec(1);

p = zeros(length(p0),n);
p(:,1) = p0;
t = 0;
% I guess use RK4 
for k = 1:n-1
    s_k = p(:,k);
    % Runge-Kutta Steps 
    k1 = coStateEq(t,p(:,k),x(:,k),mu,u(k),alph(:,k),Tmax);
    k2 = coStateEq(t+dt/2,p(:,k)+k1*dt/2,x(:,k),mu,u(k),alph(:,k),Tmax);
    k3 = coStateEq(t+dt/2,p(:,k)+k2*dt/2,x(:,k),mu,u(k),alph(:,k),Tmax);
    k4 = coStateEq(t+dt,p(:,k)+k3*dt,x(:,k),mu,u(k),alph(:,k),Tmax);
    % Calculate next position and velocity
    p_kp1 = s_k + (dt/6)*(k1 +2*k2 +2*k3 + k4);
    % Store it in the next index for next loop
    p(:,k+1) = p_kp1;
end


end

function [pdot] = coStateEq(t,p,s,mu,u,alph,Tmax)


pdot = zeros(7+6*7,1); %co state vector

lamJ = p(1:7);
lami = reshape(p(8:end),[7 6]);

m = s(7)';



G = getG(s,mu);
H = 2*[0 1 0;-1 0 0; 0 0 0];

dfds = [zeros(3)       eye(3)            zeros(3,1);
    G             H        -alph*u*Tmax/m^2;
    zeros(1,3)     zeros(1,3)            0];


lamJdot = -(lamJ'*dfds)';
lamidot = -(lami'*dfds)';


pdot(1:7) = lamJdot;
pdot(8:end) = reshape(lamidot,[],1);


end



function [sdot] = cr3bp_lowThrust_eqns(t,s,mu,u,alph,Tmax,c)
% s = [r' v' m]'

x = s(1:6);
m = s(7);

xdot = CR3BP_EOM(t,x,mu);

xdot(4:6) = xdot(4:6) + u*Tmax*alph/m;

mdot = -u*Tmax/c;

sdot = [xdot;mdot];



end

%% OLD FUNCTIONS (SHOOTING METHOD STUFF)
function [ydot,u] = stateEq(t,y,mu,Tmax,c,eppy)
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

opt = odeset("AbsTol",1e-12,"RelTol",1e-12,"Events",@(t,y) eventS0(t,y,c));

[t,y,te,ye] = ode89(@(t,y) stateEq(t,y,mu,Tmax,c,eppy),timeVec,y0,opt);

y = y';
ye = ye';


end


function [position,isterminal,direction] = eventS0(t,y,c)

% x = y(1:6);
m = y(7);
% pr = y(8:10);
pv = y(11:13);
pm = y(14);

% alph = -pv./norm(pv);

S = -c/m *norm(pv) - pm + 1; % Switching curve


position = S;
isterminal = 0;
direction =0;

end

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

[ydot,u] = stateEq(t,y,mu,Tmax,c,eppy);

vdot(1:14) = ydot;

G = getG(x,mu);
H = 2*[0 1 0;-1 0 0; 0 0 0];


% Note these are valid only for u = 0 or u = 1 (which is what i've been
% using since eppy = 0)

DyF11 = [zeros(3)       eye(3)            zeros(3,1);
            G             H        pv/norm(pv)*u*Tmax/m^2;
        zeros(1,3)     zeros(1,3)            0];

DyF12 = zeros(7); DyF12(4:6,4:6) = -u*Tmax/m*( eye(3)/norm(pv) - pv*pv'/norm(pv)^3);

DyF21 = zeros(7); DyF21(1:3,1:3) = -dGpvdr(x,pv,mu); DyF21(7,7) = 2*norm(pv)*u*Tmax/m^3;

DyF22 = -DyF11';

DyF = [DyF11 DyF12; DyF21 DyF22];

PHI = reshape(v(15:end),[14 14]);

PHIdot = DyF*PHI;

vdot(15:end) = reshape(PHIdot,[],1);





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



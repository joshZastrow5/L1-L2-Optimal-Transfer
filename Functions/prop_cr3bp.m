function [xR,t] = prop_cr3bp(x0,mu,timeVec,options)
%% Propagates Trajectory for Circular Restricted 3 Body Problem (Rotating Frame)
% Given initial conditions, vectors in the Rotating frame
% Non-Dimensional Units for all values
arguments
    x0 double
    mu (1,1) double
    timeVec
    options.noiseCov = []
    options.RelTol double = 1e-12;
    options.AbsTol double = 1e-12;

end

% x0 = [x0*1 x0*2 ... x0*N] : Matrix of initial conditions where each
% column is initial condition of [r0_R' v0_R']'
% r0_R = [x0 y0 z0]': initial position vector in rotating frame
% v0_R = [xdot0 ydot0 zdot0]': initial velocity vector in rotating frame
% timeSpan: either in the form [0 finalTime] or exact time vector
% mu: non-dimensional mass term
% noiseCov: Covariance matrix for zero-mean Gaussian additive noise
% RelTol & AbsTol: Tolerance to be used for ode options
% x_R matrix


[dim,N] = size(x0);
X0 = reshape(x0,[],1);

if isempty(options.noiseCov)
   
    opt = odeset("RelTol",options.RelTol,"AbsTol",options.AbsTol);
    [t,xR] = ode89( @(t,x) CR3BP_EOM(t,x,mu,dim),timeVec,X0,opt);

    xR = reshape(xR,length(t),dim,N);

    xR = pagetranspose(xR);

else
    %% with the inclusion of noise, ode89() is messy due to the adaptive time step stuff
    % Bc it calls the function numerous times and the noise would be different
    % for different calls corresponding to same time we have an issue.
    % To resolve, we'll generate the noise for a time step here and propagate
    % forward using RK4


    % In case timeVec doesn't have small enough timestep we will make sure it
    % does

    %% TOD0 - I'd like to implement adaptive time steps for accuracy like ode89() can do
    % I'd also like to interpolate to 

    t = linspace(timeVec(1),timeVec(end),1e4);

    dt = t(2) - t(1);

    xR = zeros(length(X0),length(t));

    xR(:,1) = X0;

    % RK4
    for i = 1:length(t)-1

        noise = mvnrnd(zeros(dim/2,1),options.noiseCov,N)'; % Generate noise for this timeStep

        k1 = CR3BP_EOM(t(i),xR(:,i),mu,dim,noise);

        k2 = CR3BP_EOM(t(i)+dt/2, xR(:,i)+k1*dt/2,mu,dim,noise);

        k3 = CR3BP_EOM(t(i)+dt/2, xR(:,i)+k2*dt/2,mu,dim,noise);

        k4 = CR3BP_EOM(t(i)+dt, xR(:,i)+k3*dt,mu,dim,noise);

        xR(:,i+1) = xR(:,i) + dt/6*( k1 + 2*k2 + 2*k3 + k4);
    end

    xR = reshape(xR',length(t),dim,N);
    xR = pagetranspose(xR);


end

end

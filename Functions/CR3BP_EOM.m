function xdot = CR3BP_EOM(t,x,mu,dim,noise)

arguments
    t double
    x double
    mu (1,1) double
    dim (1,1) int16 = 6
    noise double = zeros(dim/2,1)
end
%% Circular Restricted 3 Body Problem Equations of Motion
% Non-dimensionalized Units
%
% with the LARGER MASS, m1 to the left of the origin at (-mu,0)
% and m2, or the planet (ie. Earth), is at (1 - mu, 0)
%
%                L4
% -L3----m1--+-----L1--m2--L2-
%                L5
% From Shane Ross
%
% x = [x1; x2; x3; ...; xN]
%
% if dim = 6: xi = [X, Y, Z, Xdot, Ydot, Zdot]' [DEFAULT]
% if dim = 4: xi = [X, Y, Xdot, Ydot]'
% mu: Non-dimensional mass term, used to non-dimensionalize everything else

% If doing particle propagation, the x vector is very long column vector

N = length(x)/dim; % First we find the number of particles

x = reshape(x,dim,N); % Then reshape to propagate in tandem

mu1 = 1 - mu; % mass of larger Primary, m1
mu2 = mu; % mass of smaller primary, m2

xdot = zeros(dim,N);

switch dim
    case 6
        X = x(1,:); Y = x(2,:); Z = x(3,:);
        Xdot = x(4,:); Ydot = x(5,:); Zdot = x(6,:);

        r1 = sqrt((X + mu2).^2 + Y.^2 + Z.^2); % Distance from m1, larger primary
        r2 = sqrt((X - mu1).^2 + Y.^2 + Z.^2); % Distance from m2, smaller primary

        xdot(1:3,:) = [Xdot; Ydot; Zdot];

        xdot(4,:) = 2*Ydot + X - mu1*(X + mu2)./(r1.^3) - mu2*(X-mu1)./(r2.^3) + noise(1);
        xdot(5,:) = -2*Xdot + Y - mu1*(Y./r1.^3) - mu2*Y./r2.^3 + noise(2);
        xdot(6,:) = -mu1*Z./r1.^3 - mu2*Z./r2.^3 + noise(3);

    case 4
        X = x(1,:); Y = x(2,:); 
        Xdot = x(3,:); Ydot = x(4,:); 

        r1 = sqrt((X + mu2).^2 + Y.^2); % Distance from m1, larger primary
        r2 = sqrt((X - mu1).^2 + Y.^2); % Distance from m2, smaller primary

        xdot(1:2,:) = [Xdot; Ydot];

        xdot(3,:) = 2*Ydot + X - mu1*(X + mu2)./(r1.^3) - mu2*(X-mu1)./(r2.^3) + noise(1);
        xdot(4,:) = -2*Xdot + Y - mu1*(Y./r1.^3) - mu2*Y./r2.^3 + noise(2);


end

xdot = reshape(xdot,[],1);

end


function [xL,xU] = fibBound(x0,delta,f)
%% Fibonacci Bounding Problem
% xL,xU, are lower and upper bounds
% x0: initial guess value,
% delta: step size
% f: anonymous function in form f = @(x) ax

%% Phase 0 - Direction of delta
if f(x0) > f(x0+delta)
    delta = 1*delta;
elseif f(x0) > f(x0-delta)
    delta = -delta;
else
    xL = x0-delta;
    xU = x0+delta;
    return
end

%% Phase 1 - Bounding
x = [x0 x0+1.618*delta]; % initialize vector
for i = 2:100
    x(i+1) = x(i) + (1.618^i)*delta;
    if (f(x(i-1)) > f(x(i))) && (f(x(i+1)) > f(x(i)))
        if x(i-1) < x(i+1)
            xL = x(i-1);
            xU = x(i+1);
        else
            xL = x(i+1);
            xU = x(i-1);
        end
        break
    end
end
% fprintf('Fibonacci Bounding Table\n')
% fibBoundTab = table(x',f(x)','VariableNames',["x","f(x)"]);
% disp(fibBoundTab)
end
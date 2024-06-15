% The Newton-Raphson method is a zero-finding algorithm that uses a function and
% its derivative to rapidly converge upon the root.  In general, it is necessary to
% try multiple initial locations (x0) to find all roots.  For our case, it is
% sufficient to try in the vicinity of the expected frequency.
%
% Inputs:
%   f : function with a zero that we wish to locate.
%   df : function that is the derivative of f.
%   x0 : initial estimate of the zero location.
%   tolerance : stop early when |f(x')| < tolerance
%   max_iterations : stop after this many steps.  (Should converge quickly!)
% Outputs:
%   root : root location such that f(root) ~= 0.

%function root = newton_raphson(f, df, x0, tolerance, max_iterations)
function root = newton(f, df, x0, tolerance, max_iterations)

   for i = 1:max_iterations
       fx = f(x0);
       dfx = df(x0);
       
       if abs(fx) < tolerance
           root = x0;
           return;
       end
       
       x0 = x0 - fx / dfx;
   end
   
   root = x0;
end

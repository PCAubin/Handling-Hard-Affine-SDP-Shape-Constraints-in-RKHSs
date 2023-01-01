function [reg_min,G,reg_param] = gcv(U,s,b,method)
%GCV Plot the GCV function and find its minimum.
%
% [reg_min,G,reg_param] = gcv(U,s,b,method)
% [reg_min,G,reg_param] = gcv(U,sm,b,method)  ,  sm = [sigma,mu]
%
% Plots the GCV-function
%          || A*x - b ||^2
%    G = -------------------
%        (trace(I - A*A_I)^2
% as a function of the regularization parameter reg_param.
% Here, A_I is a matrix which produces the regularized solution.
%
% The following methods are allowed:
%    method = 'Tikh' : Tikhonov regularization   (solid line )
%    method = 'tsvd' : truncated SVD or GSVD     (o markers  )
%    method = 'dsvd' : damped SVD or GSVD        (dotted line)
% If method is not specified, 'Tikh' is default.
%
% If any output arguments are specified, then the minimum of G is
% identified and the corresponding reg. parameter reg_min is returned.

% Per Christian Hansen, UNI-C, 03/16/93.

% Reference: G. Wahba, "Spline Models for Observational Data",
% SIAM, 1990.

% Set defaults.
if (nargin==3), method='Tikh'; end  % Default method.
npoints = 100;                      % Number of points on the curve.
smin_ratio = 16*eps;                % Smallest regularization parameter.

% Initialization.
[m,n] = size(U); [p,ps] = size(s);
beta = U'*b; beta2 = b'*b - beta'*beta;
% if (nargout > 0), find_min = 1; else find_min = 0; end

  reg_param = zeros(npoints,1); G = reg_param; s2 = s.^2;
  reg_param(npoints) = max([s(p),s(1)*smin_ratio]);
%   ratio = (s(1)/reg_param(npoints))^(1/(npoints-1));
  ratio = 1.2*(s(1)/reg_param(npoints))^(1/(npoints-1));
  for i=npoints-1:-1:1, reg_param(i) = reg_param(i+1)/ratio; end
  delta0 = 0;
  if (m > n) && (beta2 > 0)
      delta0 = beta2; end
  for i=1:npoints
    f1 = (m*reg_param(i))./(s2 + m*reg_param(i));
    fb = f1.*beta(1:p); rho2 = fb'*fb + delta0;
    G(i) = rho2/m/(m - n + sum(f1)/m)^2;
  end 
  [~,minGi] = min(G);
  reg_min = reg_param(minGi);
  figure
  loglog(reg_param,G,'-'), xlabel('lambda'), ylabel('G(lambda)')
  title('GCV function')
%   if (find_min)
    [minG,minGi] = min(G); reg_min = reg_param(minGi);
    HoldState = ishold; hold on;
    loglog(reg_min,minG,'*',[reg_min,reg_min],[minG/1000,minG],':')
    title(['GCV function, minimum at ',num2str(reg_min)])
    if (~HoldState), hold off; end
%   end
end
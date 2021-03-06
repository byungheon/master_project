%% B-Spline Generator with m number of Parameters
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  params       spline coefficients                m*n     (parameter num * joints num)
%  order        order of spline basis functions    1*1
%  horizon      total time horizon                 1*1
%  t            (optional) user specified times    1*k     (1*desired number of frames)

%% Outputs
% [Name]       [Description]                      [Size]
%  sp           spline pos functions               struct n
%  spdot        (optional) spline vel functions    struct n
%  spddot       (optional) spline acc functions    struct n
%  If t is specified, return real-valued vectors   n*k

%% Implementation
function [sp, varargout] = makeSpline(params, order, horizon, varargin)
    n = size(params,2);
    m = size(params,1);
    
    num_knots = m + order;
    knots = linspace(0, horizon, num_knots);
    
    for i = 1:n
        sp(i)     = spmak(knots, params(:,i)');
        spdot(i)  = fnder(sp(i));
        spddot(i) = fnder(spdot(i));
    end
    
    if     nargin == 3
    elseif nargin == 4   % if time specified
        t = varargin{1};
        num_t = size(t,2);
        q = zeros(n, num_t);
        qdot = zeros(n, num_t);
        qddot = zeros(n, num_t);
        
        for j = 1:n
            q(j,:)     = fnval(sp(j), t);
            qdot(j,:)  = fnval(spdot(j), t);
            qddot(j,:) = fnval(spddot(j), t);
        end

        sp = q;
        spdot = qdot;
        spddot = qddot;
    end
    
    if nargout > 1
        varargout{1} = spdot;
        varargout{2} = spddot;
    end
end
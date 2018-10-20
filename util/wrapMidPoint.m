%% Wrap angle around mid point
% 2018 Byungheon Kim

%% Inputs
% [Name]     [Description]                       [Size]
%  q          joint angles                        n*dof
%  q_midpoint mid points of each joint angle      dof * 1

%% Outputs
% [Name]     [Description]                       [Size]
% q_wrapped  wrapped joint angles                n * dof

%% Examples

%% Implementation
function q_wrapped = wrapMidPoint(q,q_midpoint)
    if size(q,2) ~= size(q_midpoint,2)
        error('wrapping error: wrong input size');
    end
    
    n = size(q,1); % number of data
    
    q_mid_mat = repmat(q_midpoint, n,1);
    q_wrapped = wrapToPi(q - q_mid_mat) + q_mid_mat;
end
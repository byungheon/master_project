%% Visualizer for KUKA LWR iiwa R820 
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                    [Size]
%  trajectory   robot trajectory to visualize    7*num_time

%% Implementation
function appVisualizeKUKA(varargin)
    %% Robot Init
    robot = makeKukaR820();

    n = robot.dof;
    q = zeros(n,1);
    T = zeros(4,4,n);
    for i = 1:n
        T(:,:,i) = solveForwardKinematics(q(1:i), robot.A(:,1:i), robot.M(:,:,1:i));
    end
    
    trajectory = zeros(7,1);
    if nargin > 0
        if size(varargin{1},1) == n
            trajectory = varargin{1};
        else
            error('number of the joints and the trajectory does not match');
        end
    end
    if(nargin > 2)
        bPilco = true;
        cost = varargin{2};
        text1 = varargin{3};
        text2 = varargin{4};
    else
        bPilco = false;
    end
    
    num_time = size(trajectory, 2);
    
    %% STL Load
    fv_base = stlread(['link_0.STL']);  % base link
    fv_base.vertices = (T(1:3,1:3,1)*fv_base.vertices' + T(1:3,4,1)*ones(1,size(fv_base.vertices,1)))';
    
    for i = 1:n+1
        fv_zero{i} = stlread(['link_' num2str(i) '.STL']);
    end
    T_temp = [-1 0 0 0.04428;
              0 -1 0 -0.122;
              0 0 -1 0.3678 + 0.125;
              0 0 0 1];
    TempVectices = ones(4,size(fv_zero{n+1}.vertices,1));
    TempVectices(1:3,:) = fv_zero{n+1}.vertices';
    TempVectices2 = inverse_SE3(T_temp) * TempVectices;
%     TempVectices2(3,:) =  TempVectices2(3,:) + 0.9450 * ones(1,size(TempVectices2,2));
    fv_zero{n+1}.vertices = TempVectices2(1:3,:)';
    
    fv = fv_zero;
    j = 1;
    for i = 1:n
        fv{j}.vertices = (T(1:3,1:3,i)*fv{j}.vertices' + T(1:3,4,i)*ones(1,size(fv{j}.vertices,1)))';
        if i == n-1
            j = j + 1;
            T_temp = T(:,:,i) * [1  0  0  0
                      0  0  1  0
                      0  -1 0  0
                      0  0  0  1];
            fv{j}.vertices = (T_temp(1:3,1:3)*fv{j}.vertices' + T_temp(1:3,4)*ones(1,size(fv{j}.vertices,1)))';
        end
        j = j+1;
    end

    %% Render
    % The model is rendered with a PATCH graphics object. We also add some dynamic
    % lighting, and adjust the material properties to change the specular
    % highlighting.
    figure('Name','KUKA LWR iiwa R820','NumberTitle','off','units','pixels','pos',[100 0 1000 1000]);
    hold on;
    axis equal;
    axis([-1 1 -1 1 -0.5 1.1]);
    xlabel('x'); ylabel('y'); zlabel('z');
    if bPilco
        reward = zeros(num_time,1);
        for i = 1:num_time
            reward(i) = 1-cost.fcn(cost, [trajectory(1:6,i); zeros(7,1); trajectory(7,i)],zeros(14));
        end
        text(0,-0.7, text1,'fontsize', 12);
        text(0,-0.9, text2,'fontsize', 12);
    end
        
    % draw base link
    patch(fv_base,'FaceColor',       [0.7 0.7 0.7], ...
                 'EdgeColor',       'none',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);

    % draw 7 links
    color = ones(n+1,3)*0.8;
    color(2,:) = [246 120 40]/255;
    color(6,:) = [246 120 40]/255;
    color(8,:) = ones(1,3)*0.2;

    for i = 1:8
        render_part{i} = patch(fv{i},'FaceColor',  color(i,:), ...
                 'EdgeColor',       'none',        ...
                 'FaceLighting',    'gouraud',     ...
                 'AmbientStrength', 0.15);
    end
    
    % draw end-effector
    end_effector_M = eye(4);
    end_effector_M(3,4) = 0.125;
    end_effector_T = T(:,:,7) * end_effector_M;
    end_effector = plot_SE3(end_effector_T);

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    % Fix the axes scaling, and set a nice view angle
    view([-135 35]);
    getframe;

    %% Animation Loop
%     while(1)
        for time = 1:num_time
            j = 1;
            T0 = eye(4,4);
            for i = 1:n
                T0 = T0 * inverse_SE3(robot.M(:,:,i)) * exp_se3(robot.A(:,i)*trajectory(i,time));
                T(:,:,i) = T0;
                fv{j}.vertices = (T(1:3,1:3,i)*fv_zero{j}.vertices' + T(1:3,4,i)*ones(1,size(fv_zero{j}.vertices,1)))';
                set(render_part{j}, 'Vertices', fv{j}.vertices);
                if i == n-1
                    j = j+1;
                    T_temp = T(:,:,i) * [1  0  0  0
                      0  0  1  0
                      0  -1 0  0
                      0  0  0  1];
                    fv{j}.vertices = (T_temp(1:3,1:3)*fv_zero{j}.vertices' + T_temp(1:3,4)*ones(1,size(fv_zero{j}.vertices,1)))';
                    set(render_part{j}, 'Vertices', fv{j}.vertices);
                end
                j = j+1;
            end

            end_effector_T = T(:,:,7) * end_effector_M;
            plot_SE3(end_effector_T, end_effector);
            
            getframe;
        end
%     end
end
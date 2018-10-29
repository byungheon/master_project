%% Visualizer for KUKA LWR iiwa R820 
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                    [Size]
%  trajectory   robot trajectory to visualize    7*num_time

%% Implementation
function appVisualizeKUKA_planar(varargin)
    %% Robot Init
    robot = makeKukaR820_planar();
    dt_standard = 0.1;
    dt = dt_standard;
    n = robot.dof;
    T = zeros(4,4,8);
    
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
        dt    = varargin{5};
    else
        bPilco = false;
    end
    if dt_standard/dt < 1
        ratio = 1;
    else
        ratio = dt_standard/dt;
    end
    list = 1:ratio:size(trajectory,2);
    trajectory = trajectory(:,list);
    num_time = size(trajectory, 2);
    
    %% STL Load
    fv_base = stlread(['link_0.STL']);  % base link
    T_vase = robot.M_zero(:,:,1) * exp_se3(robot.A_zero(:,1)*robot.q_zero(1));
    fv_base.vertices = (T_vase(1:3,1:3)*fv_base.vertices' + T_vase(1:3,4)*ones(1,size(fv_base.vertices,1)))';
    
    for i = 1:8
        fv_zero{i}  = stlread(['link_' num2str(i) '.STL']);
    end
    T_temp = [-1 0 0 0.04428;
              0 -1 0 -0.122;
              0 0 -1 0.3678 + 0.125;
              0 0 0 1];
    TempVectices = ones(4,size(fv_zero{8}.vertices,1));
    TempVectices(1:3,:) = fv_zero{8}.vertices';
    TempVectices2 = inverse_SE3(T_temp) * TempVectices;
    fv_zero{8}.vertices = TempVectices2(1:3,:)';
    
    T0 = eye(4,4);
    fv = fv_zero;
    j = 1;
    for i = 1:7
        T0 = T0 * robot.M_zero(:,:,i) * exp_se3(robot.A_zero(:,i)*robot.q_zero(i));
        fv{j}.vertices = (T0(1:3,1:3)*fv{j}.vertices' + T0(1:3,4)*ones(1,size(fv{j}.vertices,1)))';
        if i == 6
            j = j + 1;
            T_temp = T0 * [1  0  0  0
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
    if ~ishghandle(1)
        figure('Name','KUKA LWR iiwa R820','NumberTitle','off','units','pixels','pos',[100 0 1000 1000]);
    else
        set(0,'CurrentFigure',1);
    end
    clf(1);
    hold on;
%     view([0 0]);
    axis equal;
    axis([-1 1 -1 1 -1.3 1.1]);
    xlabel('x'); ylabel('y'); zlabel('z');
    if bPilco
        reward = zeros(num_time,1);
        for i = 1:num_time
            reward(i) = 1-cost.fcn(cost, [trajectory(1:2,i); zeros(3,1); trajectory(3,i) + trajectory(1,i) - trajectory(2,i)],zeros(6));
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
    color = ones(8,3)*0.8;
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
    end_effector_T = T0 * end_effector_M;
    end_effector = plot_SE3(end_effector_T);

    % Add a camera light, and tone down the specular highlighting
    camlight('headlight');
    material('dull');

    % Fix the axes scaling, and set a nice view angle
    view([-15 15]);
    getframe;
    
    temp_traj               = zeros(7,num_time);
    temp_traj([2,4,7],:)    = trajectory;
    modified_traj = temp_traj + repmat(robot.q_zero,1,num_time);
    %% Animation Loop
%     while(1)
        for time = 1:num_time
            j = 1;
            T0 = eye(4,4);
            for i = 1:7
                T0 = T0 * robot.M_zero(:,:,i) * exp_se3(robot.A_zero(:,i)*modified_traj(i,time));
                T(:,:,i) = T0;
                fv{j}.vertices = (T(1:3,1:3,i)*fv_zero{j}.vertices' + T(1:3,4,i)*ones(1,size(fv_zero{j}.vertices,1)))';
                set(render_part{j}, 'Vertices', fv{j}.vertices);
                if i == 6
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
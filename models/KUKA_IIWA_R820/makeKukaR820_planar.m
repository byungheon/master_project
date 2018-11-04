function robot = makeKukaR820_planar()
    
    % degrees of freedom
    robot.dof = 3;
    
    robot.q_zero = [0;pi/8;0;-3*pi/8;pi/2;pi/2;0];
    q = robot.q_zero;
    % screws A_i, i-th screw described in i-th frame
    robot.A_zero =       [0 0 0 0 0 0 0 
               0 0 0 0 0 0 0 
               1 1 1 1 1 1 1 
               0 0 0 0 0 0 0 
               0 0 0 0 0 0 0 
               0 0 0 0 0 0 0 ];
    robot.A = robot.A_zero(1:6,1:robot.dof);
    % link frames M_{i,i-1}
    M(:,:,1) = [1  0  0  0
                      0  1  0  0
                      0  0  1  0
                      0  0  0  1];
     
    M(:,:,2) = [1  0  0  0
                      0  0  -1 0
                      0  1  0  0
                      0  0  0  1];
     
    M(:,:,3) = [1  0  0  0
                      0  0  1  0.42
                      0  -1 0  0
                      0  0  0  1];
     
    M(:,:,4) = [1  0  0  0
                      0  0  1  0
                      0  -1 0  0
                      0  0  0  1]; 

    M(:,:,5) = [1  0  0  0
                      0  0  -1 -0.4
                      0  1  0  0
                      0  0  0  1];

    M(:,:,6) = [1  0  0  0
                      0  0  -1 0
                      0  1  0  0
                      0  0  0  1];
     
    M(:,:,7) = [1  0  0  0
                      0  0  1  0
                      0  -1 0  0
                      0  0  0  1];
    robot.M_zero   = M;
    robot.M(:,:,1) = M(:,:,1) * exp_se3(robot.A_zero(:,1)*q(1)) * M(:,:,2) * exp_se3(robot.A_zero(:,2)*q(2));            
    robot.M(:,:,2) = M(:,:,3) * exp_se3(robot.A_zero(:,3)*q(3)) * M(:,:,4) * exp_se3(robot.A_zero(:,4)*q(4));
    robot.M(:,:,3) = M(:,:,5) * exp_se3(robot.A_zero(:,5)*q(5)) * M(:,:,6) * exp_se3(robot.A_zero(:,6)*q(6)) * M(:,:,7) * exp_se3(robot.A_zero(:,7)*q(7));

    for i = 1:robot.dof
        robot.M(:,:,i) = inverse_SE3(robot.M(:,:,i));
    end
    
    % inertia matrix Phi for linear dynamics tau = Y * Phi
%     robot.Phi = [4.10357488000000,0.000381222106352000,0.138656348192301,-0.330854499988890,0.0642805199818296,0.0598897338632779,0.0169802769496190,-1.25611747470647e-05,0.00626773563684459,3.02263830489678e-05,3.94661137000000,-0.000575613268314500,0.232279706554808,0.165660511968071,0.0511613065643967,0.0179414068940495,0.0440363075767375,6.47379952010187e-05,-0.00487830894892654,4.63415856705431e-05,3.17014730000000,-0.00289577105118500,-0.0936995999254780,-0.284062031562163,0.0532424223544108,0.0492097915394241,0.0103580310991201,-0.000129659899551928,-0.00353625169819843,-0.000417566462730458,2.73577100000000,8.67239407000000e-06,-0.184536973498270,0.0937214136906700,0.0328276173863488,0.00933210666178131,0.0288109507795505,5.84982205989516e-07,0.00348231519459983,-2.67096881399424e-07,1.69035100000000,0.000210972708310000,0.0361473278630100,-0.236740773927710,0.0439802849613939,0.0418353983987378,0.00522597922566354,-6.05154799058228e-06,0.00198113544746485,3.03476159939175e-05,1.80624790000000,-4.69624454000000e-06,0.00109610347563600,0.000571280085812000,0.00488548584389870,0.00355788069667578,0.00467648517164339,8.28498690366536e-08,0.000242193324392726,-1.38514671776889e-07,0.307396530000000,1.87819279830000e-06,0,0.0299208500849349,0.00313112573774775,0.00312784574922351,0.000333610011475758,0,0,-2.32816394018952e-07]';
    Phi = [4.10357488000000,0.000381222106352000,0.138656348192301,-0.330854499988890,0.0642805199818296,0.0598897338632779,0.0169802769496190,-1.25611747470647e-05,0.00626773563684459,3.02263830489678e-05,3.94661137000000,-0.000575613268314500,0.232279706554808,0.165660511968071,0.0511613065643967,0.0179414068940495,0.0440363075767375,6.47379952010187e-05,-0.00487830894892654,4.63415856705431e-05,3.17014730000000,-0.00289577105118500,-0.0936995999254780,-0.284062031562163,0.0532424223544108,0.0492097915394241,0.0103580310991201,-0.000129659899551928,-0.00353625169819843,-0.000417566462730458,2.73577100000000,8.67239407000000e-06,-0.184536973498270,0.0937214136906700,0.0328276173863488,0.00933210666178131,0.0288109507795505,5.84982205989516e-07,0.00348231519459983,-2.67096881399424e-07,1.69035100000000,0.000210972708310000,0.0361473278630100,-0.236740773927710,0.0439802849613939,0.0418353983987378,0.00522597922566354,-6.05154799058228e-06,0.00198113544746485,3.03476159939175e-05,1.80624790000000,-4.69624454000000e-06,0.00109610347563600,0.000571280085812000,0.00488548584389870,0.00355788069667578,0.00467648517164339,8.28498690366536e-08,0.000242193324392726,-1.38514671776889e-07,0.3000 ,        0 ,   0.0750  ,  0.0375 ,   0.0234   , 0.0047 ,   0.0187     ,    0  , -0.0094    ,     0]';
    
    T               = inverse_SE3(M(:,:,3) * exp_se3(robot.A_zero(:,3) *q(3)));
    robot.Phi(:,1)  = Phi(11:20) + convertInertiaGToPhi(large_Ad(T)' * convertInertiaPhiToG(Phi(21:30)) * large_Ad(T));
    T5              = inverse_SE3(M(:,:,5) * exp_se3(robot.A_zero(:,5) *q(5)));
    T6              = inverse_SE3(M(:,:,6) * exp_se3(robot.A_zero(:,6) *q(6))) * T5;
    robot.Phi(:,2)  = Phi(31:40) + convertInertiaGToPhi(large_Ad(T5)' * convertInertiaPhiToG(Phi(41:50)) * large_Ad(T5)) + convertInertiaGToPhi(large_Ad(T6)' * convertInertiaPhiToG(Phi(51:60)) * large_Ad(T6));
    robot.Phi(:,3)  = Phi(61:70);
    robot.G = zeros(6,6,robot.dof);
    
    for i = 1:robot.dof
        robot.G(:,:,i) = convertInertiaPhiToG(robot.Phi(:,i));
    end
    
    F = [0.335546 0.181378;0.161421 0.355506;0.083204 0.144099; 0.206558 0.070046; 0.085193 0.040282; 0.101937 0.136860;0.02 0.02];
    F(:,1) = zeros(7,1);
    robot.F = F([2,4,7],:);
    robot.F(1:2,:) = robot.F(1:2,:);
    robot.Sigmoid = 100 * ones(robot.dof,1);

end
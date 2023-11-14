%% 2.74 Final Project Simulation 
function simulate_arm()
    global hit
    %% Definte fixed paramters

    m1 = 0.154;     % CAD forearm is 5.443 oz = 0.154 kg          
    m2 =  0.145;     % CAD hand is 145 g = .145 kg  
    m3 = 0.26;      % volleyball in kg 

    l1 = 0.257;      % CAD forearm length in 10.12 in  =  0.257 m         
    l2 = 0.173;      % CAD hand length in 6.8 in = 0.173 m 

    l1c = 0.125;      % CAD forearm COM length in 5.293 in = 0.125 m       
    l2c = 0.091;      % CAD hand COM length 3.60 in in = 0.091

    I = @(m, l) (1/3)*m*l^2; 

   
    % rotational inertias 
    Ir1 = I(m1, l1);       
    Ir2 = I(m2, l2);

    % ball position parameters
    ball_x_i = -.2; 
    ball_y_i = 0.3; 

    m_B = 0.26;  % volleyball in kg
    r = .1;

    g = -9.81;    
    
    %% Parameter vector
    p   = [m1 m2 Ir1 Ir2 l1 l2 l1c l2c g m3 ball_x_i ball_y_i]'; %need to add m3 to parameters
    hit = false;
    p_B = [m_B g r]';

    %% Perform Dynamic simulation
    dt = 0.001;
    tf = 2;
    num_step = floor(tf/dt);
    tspan = linspace(0, tf, num_step); 

    % arm state vector 
    % z = [theta1, theta2, theta1_dot, theta2_dot] 
    z0 = [pi/12+1; 0; 0; 0];
    z_out = zeros(4,num_step);
    z_out(:,1) = z0;

    % ball state vector 
    % b = [x, y, dx, dy] 
    b0 = [ball_x_i, ball_y_i, 0, 0];
    b_out = zeros(4,num_step);
    b_out(:,1) = b0;

    for i=1:num_step-1
        % dynamics of the system -- how will change in each timestep 
        dz = dynamics(tspan(i), z_out(:,i), p, p_B);
        
        % updates the velocities & positions 
        z_out(3:4,i+1) = z_out(3:4,i) + dz(3:4)*dt;             %velocity
        z_out(1:2,i+1) = z_out(1:2,i) + z_out(3:4,i+1)*dt;      %position

        % ball 
        db = dynamicsB(tspan(i), z_out(:, i), b_out(:, i), p, p_B);

        b_out(3:4,i+1) = b_out(3:4,i) + db(3:4)*dt;             %velocity
        b_out(1:2,i+1) = b_out(1:2,i) + b_out(3:4,i+1)*dt;      %position
        
        % % hand max of 180 degree range of motion 
        % if (z_out(2, i+1) < -pi/2)
        %     z_out(2, i+1) = -pi/2; 
        % elseif (z_out(2, i+1) > pi/2)
        %     z_out(2, i+1) = pi/2; 
        % end         
    end

    %% Below are different graphs to see the energy, position, trajectory, and velocity of the hand and the ball 
    % Some of the graphs are commented so that only a few figures are shown
    % and it doesn't get too crowded

    %% Compute Energy of Hand/Arm System
    E = energy_geometry(z_out,p);
    figure(1); clf
    plot(tspan,E);xlabel('Time (s)'); ylabel('Energy (J)');
    title('Energy of Arm (J)');
    
%     %% Compute Hand Position
%     rE = zeros(2,length(tspan));
%     for i = 1:length(tspan)
%         rE(:,i) = position_hand(z_out(:,i),p);
%     end
%     figure(2); clf;
%     plot(tspan,rE(1:2,:))
%     xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','y'});
%     title('Position vs Time of Arm');
% 
%     %% Trajectory of Hand
%     figure(3); clf;
%     plot(rE(1,:),rE(2,:))
%     xlabel('Time (s)'); ylabel('Position (m)');
%     title('Trajectory of Hand');

    %% Compute Energy of Ball
    E_ball = energy_ball(b_out, p_B);
    figure(4); clf;
    plot(tspan,E_ball);xlabel('Time (s)'); ylabel('Energy (J)');
    title('Energy of Ball (J)');

    %% Compute Ball Position (x and y) 
    rB = zeros(2,length(tspan));
    for i = 1:length(tspan)
        rB(:,i) = position_ball(b_out(:,i),p_B);
    end
%     figure(5); clf;
%     plot(tspan,rB(1:2,:))
%     xlabel('Time (s)'); ylabel('Position (m)'); legend({'x','y'});
%     title('Position vs Time of Ball');
% 
%     %% Compute Ball Velocity (dx and dy) 
%     rV = zeros(2,length(tspan));
%     for i = 1:length(tspan)
%         rV(:,i) = velocity_ball(b_out(:,i),p_B);
%     end
%     figure(6); clf;
%     plot(tspan,rB(1:2,:))
%     xlabel('Time (s)'); ylabel('Velocity (m/s)'); legend({'x','y'});
%     title('Velocity of Ball');

    %% Trajectory of Ball
    figure(7); clf;
    plot(rB(1,:), rB(2,:));
    xlabel('x'); ylabel('y'); ;
    title('Trajectory of Ball');

     %% Animate Solution
     figure(8); clf;
     hold on;
     animateSol(tspan,z_out, b_out, p, p_B);
end

function tau = control_law(t,z,p)
    % Controller gains
    K_q1 = 30; % Spring stiffness q1
    K_q2 = 0.1; % Spring stiffness q2
    D_q1 = 2;  % Damping q1
    D_q2 = 0.05;  % Damping q2
    
    % angles 
    th1 = z(1, :);
    th2 = z(2, :);

    % angular velocities 
    th1dot = z(3, :); 
    th2dot = z(4, :); 

    % desired angles 
    theta1_des = 3*pi/4;
    theta2_des = 0; 

    % forearm control law 
    tau1 = K_q1 * (theta1_des - th1) + D_q1 * (-th1dot);

    % wrist control law 
    tauJ = K_q2*(theta2_des-th2) + D_q2*(-th2dot); 

    % final torques 
    tau = [tau1, tauJ]';

    % if (limit):
        % modify tau:

end  

function Fc = contact_force(z,p,b, p_B)
    global hit;
    %% Fixed parameters for contact -- make these volleyball parameters? 
    % ball characteristics 
    K_c = 1000;
    D_c = 2;

    % ball starting position 
    y_B = b(2); 
    x_B = b(1);  
    r = p_B(end); 
    g = p_B(2);
    mb = p_B(1);
    
    r_E = position_wrist_com(z, p); % position of the COM wrist 
    y = dot(r_E, [0; 1; 0]); % vertical component of the COM wrist 
    x = dot(r_E, [1; 0; 0]); % horizontal component of the COM wrist
    
    C = y - y_B; 
    C_dot = velocity_wrist_vert(z, p); 
   
    H = x - x_B; 
    H_dot = velocity_wrist_hor(z, p); 

    % distance between vall and hand 
    d = sqrt(C^2 + H^2); 

    Fcy = 0; 
    Fcx = 0; 
    
    % Contact with hand
    if d < r && C > 0 && H > 0
        hit = true;
        Fcy = -K_c*C - D_c*C_dot; 
        Fcx = -K_c*H - D_c*H_dot; 
    end

    % % Contact with floor
    K_floor = 1000; % Floor characteristics (stiffness and damping)
    D_floor = 2;

    if y_B -r >= 0 % if the ball is above the floor
        Fc = [[Fcx, Fcy]', [1, 1]']; 
    else %ball hits the floor 
        F_floor = (K_floor*b(2)+D_floor*b(4));
        Fc = [[0, F_floor]', [0, 0]'];
    end
    
end

function dz = dynamics(t,z,p, p_B)   
    % Get mass matrix
    A = A_geometry(z,p);
    
    % Compute Controls
    tau = control_law(t,z, p);
    
    % Get b = Q - V(q,qd) - G(q)
    b = b_geometry(z,tau,p);
     
    % Solve for qdd.
    qdd = A\(b);
    dz = 0*z;
    
    % Form dz
    dz(1) = 1;
    dz(2) = z(4);
    dz(1:2) = z(3:4);
    dz(3:4) = qdd;
end

function db = dynamicsB(t,z,b, p, p_B)
    global hit;
    % % Compute the contact force (used for problem 2)
    Fc = contact_force(z,p, b, p_B); % this is a 1 by 2 -- think before it was a
    m3 = p_B(1); 
    
    F = Fc(:, 1);
    accel = F/m3 + [0; p_B(2)];

    if hit == false
        accel = [0; 0];
    end 

    db = 0*b; 

    db(1) = b(3);

    if Fc(:, 2) == [0,0] 
        db(2) = -b(4);
    else 
        db(2) = b(4);
    end 

    db(3) = accel(1); 
    db(4) = accel(2); 

end

% this function plots the ball 
function h = circle(x,y,p,n)
    % hold on
    th = 0:pi/50:2*pi;
    r = p(end);
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    if n == 1
        h = plot(xunit, yunit, 'g');
    else 
        h = plot(xunit, yunit, 'w');
    end 
end 

function animateSol(tspan, x, b, p, p_B)
    % Prepare plot handles
    h_AB = plot([0],[0],'LineWidth',2);
    h_BC = plot([0],[0],'LineWidth',2);

    yline(0); % Draw a line where the floor would be
    
    xlabel('x'); ylabel('y');
    h_title = title('t=0.0s');
    
    axis equal
    axis([-1.5 0.5 -0.1 0.5]);

    h = circle(b(1, 1), b(2, 1), p_B, 1);  %plot initial ball

    %Step through and update animation
    for i = 1:length(tspan)
        % skip frame.
        if mod(i,10)
            continue;
        end

        t = tspan(i);
        z = x(:,i); 
        keypoints = keypoints_geometry(z,p);

        rB = keypoints(:,1); % Vector to bottom of the wrist 
        rC = keypoints(:,2); % Vector to tip of the hand 

        set(h_title,'String',  sprintf('t=%.2f',t) ); % update title
        
        % Forearm Display 
        set(h_AB, 'XData',[0 rB(1)]);
        set(h_AB,'YData',[0 rB(2)]);
        
        % Hand Display 
        set(h_BC,'XData',[rB(1) rC(1)]);
        set(h_BC,'YData',[rB(2) rC(2)]);

        % Update and plot ball
        % h = circle(b(1, i-1), b(2, i-1), p_B, 0);
        ball.XData = b(1, i);         %change x coordinate of the ball
        ball.YData = b(2, i);         %change y coordinate of the ball
        h = circle(ball.XData, ball.YData, p_B, 1);
        drawnow
        pause(0.05)  %control speed, if desired

        % plots COM of the ball 
        % plot(ball.XData, ball.YData, 'b.')

        % turns off the plot of the ball 
        h = circle(ball.XData, ball.YData, p_B, 0);

    end
end
% Main script for rigid body dynamics with separate plots and quaternion-based animation
clear all; close all;

% Create animations directory if it doesn't exist
if ~exist('animations', 'dir')
    mkdir('animations');
end

%% Part 1: Define brick dimensions and calculate inertia matrix
x = 0.05;  % 5 cm
y = 0.10;  % 10 cm
z = 0.20;  % 20 cm
rho = 2000; % 2 kg/m^3

% Calculate mass
m = rho * x * y * z;

% Calculate inertia matrix
Jxx = (1/12) * m * (y^2 + z^2);
Jyy = (1/12) * m * (x^2 + z^2);
Jzz = (1/12) * m * (x^2 + y^2);
J = diag([Jxx, Jyy, Jzz]);
disp('Inertia matrix:');
disp(J);

%% Part 2: Calculate characteristic angular velocities
w3 = 2*pi;  % angular velocity for T = 1s
K0 = 0.5 * Jzz * w3^2;  % kinetic energy
w1 = sqrt(2*K0/Jxx);
w2 = sqrt(2*K0/Jyy);

% Scale factor for visualization
scale = 1.0;  % Increased scale for better visibility

%% Part 3: Calculate initial conditions for three cases
% Case 1: Near x-axis
w1_case1 = w1 * 0.99;
w3_case1 = w3 * 0.1;
omega1 = [w1_case1; 0; w3_case1];

% Case 2: Near y-axis
w2_case2 = w2 * 0.99;
w3_case2 = w3 * 0.1;
omega2 = [0; w2_case2; w3_case2];

% Case 3: Near z-axis
w1_case3 = w1 * 0.1;
w3_case3 = w3 * 0.99;
omega3 = [w1_case3; 0; w3_case3];

%% Part 4: Solve Euler equations
% Time settings
t_end = 10;
dt = 0.01;
t = 0:dt:t_end;

% Solve with tighter tolerances
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t1, w1] = ode45(@(t,w) euler_eqs(t, w, J), [0 t_end], omega1, options);
[t2, w2] = ode45(@(t,w) euler_eqs(t, w, J), [0 t_end], omega2, options);
[t3, w3] = ode45(@(t,w) euler_eqs(t, w, J), [0 t_end], omega3, options);

%% Part 5: Create plots for each case
% Case 1: Near x-axis
figure('Position', [50, 500, 1200, 400]);

% Plot ellipsoids
subplot(1, 2, 1);
[ke_x1, ke_y1, ke_z1, am_x1, am_y1, am_z1] = calculate_ellipsoids(omega1, J, scale);
plot_ellipsoids(ke_x1, ke_y1, ke_z1, am_x1, am_y1, am_z1, scale, 'Case 1: Near X-axis - Ellipsoids');

% Plot trajectory
subplot(1, 2, 2);
plot_trajectory(w1, scale, 'Case 1: Near X-axis - Trajectory');

% Case 2: Near y-axis
figure('Position', [50, 300, 1200, 400]);

% Plot ellipsoids
subplot(1, 2, 1);
[ke_x2, ke_y2, ke_z2, am_x2, am_y2, am_z2] = calculate_ellipsoids(omega2, J, scale);
plot_ellipsoids(ke_x2, ke_y2, ke_z2, am_x2, am_y2, am_z2, scale, 'Case 2: Near Y-axis - Ellipsoids');

% Plot trajectory
subplot(1, 2, 2);
plot_trajectory(w2, scale, 'Case 2: Near Y-axis - Trajectory');

% Case 3: Near z-axis
figure('Position', [50, 100, 1200, 400]);

% Plot ellipsoids
subplot(1, 2, 1);
[ke_x3, ke_y3, ke_z3, am_x3, am_y3, am_z3] = calculate_ellipsoids(omega3, J, scale);
plot_ellipsoids(ke_x3, ke_y3, ke_z3, am_x3, am_y3, am_z3, scale, 'Case 3: Near Z-axis - Ellipsoids');

% Plot trajectory
subplot(1, 2, 2);
plot_trajectory(w3, scale, 'Case 3: Near Z-axis - Trajectory');

%% Part 6: Create animations
animate_case(t1, w1, x, y, z, 'Case 1 - Near X-axis');
animate_case(t2, w2, x, y, z, 'Case 2 - Near Y-axis');
animate_case(t3, w3, x, y, z, 'Case 3 - Near Z-axis');

%% Helper Functions

% Euler equations
function dw = euler_eqs(t, w, J)
    dw = zeros(3,1);
    dw(1) = ((J(2,2) - J(3,3)) * w(2) * w(3)) / J(1,1);
    dw(2) = ((J(3,3) - J(1,1)) * w(3) * w(1)) / J(2,2);
    dw(3) = ((J(1,1) - J(2,2)) * w(1) * w(2)) / J(3,3);
end

% Calculate ellipsoids
function [ke_x, ke_y, ke_z, am_x, am_y, am_z] = calculate_ellipsoids(w_init, J, scale)
    % Calculate energy and momentum
    K = 0.5 * w_init' * J * w_init;
    H = J * w_init;
    H0 = norm(H);
    
    % Generate sphere points
    [X, Y, Z] = sphere(50);
    
    % Calculate ellipsoid surfaces
    ke_x = scale * X * sqrt(2*K/J(1,1));
    ke_y = scale * Y * sqrt(2*K/J(2,2));
    ke_z = scale * Z * sqrt(2*K/J(3,3));
    
    am_x = scale * X * H0/J(1,1);
    am_y = scale * Y * H0/J(2,2);
    am_z = scale * Z * H0/J(3,3);
    
    % Print diagnostic information
    fprintf('\nCase Analysis:\n');
    fprintf('Initial ω = [%.2f, %.2f, %.2f]\n', w_init(1), w_init(2), w_init(3));
    fprintf('Kinetic Energy: %.4f\n', K);
    fprintf('Angular Momentum: %.4f\n', H0);
    fprintf('Principal Moments: %.4f, %.4f, %.4f\n', J(1,1), J(2,2), J(3,3));
end

% Plot ellipsoids
function plot_ellipsoids(ke_x, ke_y, ke_z, am_x, am_y, am_z, scale, titleStr)
    hold on;
    
    % Plot surfaces
    surf(ke_x, ke_y, ke_z, 'FaceColor', 'b', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.1);
    surf(am_x, am_y, am_z, 'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeAlpha', 0.1);
    
    % Set axis limits
    max_val = max([max(max(abs(ke_x))), max(max(abs(ke_y))), max(max(abs(ke_z))), ...
                   max(max(abs(am_x))), max(max(abs(am_y))), max(max(abs(am_z)))]);
    axis([-max_val max_val -max_val max_val -max_val max_val]);
    
    % Add coordinate axes
    axisLength = scale;
    quiver3(0, 0, 0, axisLength, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(0, 0, 0, 0, axisLength, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(0, 0, 0, 0, 0, axisLength, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    
    % Labels and formatting
    xlabel('\omega_x', 'FontSize', 12);
    ylabel('\omega_y', 'FontSize', 12);
    zlabel('\omega_z', 'FontSize', 12);
    title(titleStr, 'FontSize', 14);
    legend('Kinetic Energy Ellipsoid', 'Angular Momentum Ellipsoid', 'Location', 'best');
    grid on;
    axis equal;
    view([-45, 30]);
end

% Plot trajectory
function plot_trajectory(w, scale, titleStr)
    hold on;
    
    % Plot path
    plot3(w(:,1), w(:,2), w(:,3), 'k', 'LineWidth', 2);
    plot3(w(1,1), w(1,2), w(1,3), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    plot3(w(end,1), w(end,2), w(end,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
    
    % Add coordinate axes
    axisLength = scale;
    quiver3(0, 0, 0, axisLength, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(0, 0, 0, 0, axisLength, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(0, 0, 0, 0, 0, axisLength, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    
    % Labels and formatting
    xlabel('\omega_x', 'FontSize', 12);
    ylabel('\omega_y', 'FontSize', 12);
    zlabel('\omega_z', 'FontSize', 12);
    title(titleStr, 'FontSize', 14);
    legend('Trajectory', 'Start', 'End', 'Location', 'best');
    grid on;
    axis equal;
    view([-45, 30]);
    
    % Add initial conditions text
    initialw = w(1,:);
    str = sprintf('Initial Angular Velocities:\n\\omega_x = %.2f\n\\omega_y = %.2f\n\\omega_z = %.2f', ...
                  initialw(1), initialw(2), initialw(3));
    annotation('textbox', [0.02, 0.90, 0.15, 0.08], 'String', str, ...
              'EdgeColor', 'none', 'BackgroundColor', 'none', 'FontSize', 8);
end

% Animate case with quaternion-based rotation
function animate_case(t, w, x, y, z, titleStr)
    % Create figure
    fig = figure('Position', [100, 100, 800, 600]);
    
    % Set up video writer
    cleanTitle = strrep(titleStr, ' - ', '_');
    cleanTitle = strrep(cleanTitle, ' ', '_');
    videoFile = fullfile('animations', [cleanTitle '.mp4']);
    
    % Create video writer object
    v = VideoWriter(videoFile, 'MPEG-4');
    v.FrameRate = 60;
    v.Quality = 95;
    open(v);
    
    % Create brick vertices and faces
    vertices = [-x/2 -y/2 -z/2;
                x/2 -y/2 -z/2;
                x/2 y/2 -z/2;
                -x/2 y/2 -z/2;
                -x/2 -y/2 z/2;
                x/2 -y/2 z/2;
                x/2 y/2 z/2;
                -x/2 y/2 z/2];
    
    faces = [1 2 3 4;  % bottom
             5 6 7 8;  % top
             1 2 6 5;  % front
             2 3 7 6;  % right
             3 4 8 7;  % back
             4 1 5 8]; % left
             
    faceColors = [0.8 0.8 1.0;    % bottom - light blue
                  0.8 0.8 1.0;    % top - light blue
                  1.0 0.8 0.8;    % front - light red
                  0.8 1.0 0.8;    % right - light green
                  1.0 0.8 0.8;    % back - light red
                  0.8 1.0 0.8];   % left - light green
    
    % Create patch object
    p = patch('Vertices', vertices, 'Faces', faces, ...
              'FaceColor', 'flat', 'FaceVertexCData', faceColors, ...
              'EdgeColor', 'k', 'LineWidth', 1.5);
    
    % Set up view
    axis equal
    axis([-0.3 0.3 -0.3 0.3 -0.3 0.3])
    grid on
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(titleStr);
    view(3)
    
    % Add coordinate axes
    hold on
    axisLength = 0.15;
    quiver3(0, 0, 0, axisLength, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(0, 0, 0, 0, axisLength, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
    quiver3(0, 0, 0, 0, 0, axisLength, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
   

    % Add axis labels
    text(axisLength*1.1, 0, 0, 'X', 'Color', 'r', 'FontWeight', 'bold');
    text(0, axisLength*1.1, 0, 'Y', 'Color', 'g', 'FontWeight', 'bold');
    text(0, 0, axisLength*1.1, 'Z', 'Color', 'b', 'FontWeight', 'bold');
    
    % Add body-fixed axes
    bodyAxisLength = min([x, y, z]);
    hBodyX = quiver3(0, 0, 0, bodyAxisLength, 0, 0, 'r--', 'LineWidth', 1);
    hBodyY = quiver3(0, 0, 0, 0, bodyAxisLength, 0, 'g--', 'LineWidth', 1);
    hBodyZ = quiver3(0, 0, 0, 0, 0, bodyAxisLength, 'b--', 'LineWidth', 1);
    
    % Add info display
    hTime = text(-0.3, -0.3, 0.3, '', 'FontSize', 10);
    
    % Initialize quaternion
    q = [1; 0; 0; 0];  % Initial quaternion [w; x; y; z]
    
    % Calculate frames
    num_frames = 600;
    t_video = linspace(0, t(end), num_frames);
    
    % Interpolate angular velocities
    wx = interp1(t, w(:,1), t_video, 'pchip');
    wy = interp1(t, w(:,2), t_video, 'pchip');
    wz = interp1(t, w(:,3), t_video, 'pchip');
    
    dt = t(end)/num_frames;
    
    % Animation loop
    for i = 1:num_frames
        if ~ishandle(fig)
            break
        end
        
        % Get current angular velocity
        omega = [wx(i); wy(i); wz(i)];
        
        % Update quaternion with RK4 integration
        q = update_quaternion(q, omega, dt);
        
        % Convert quaternion to rotation matrix
        R = quat2rotm(q);
        
        % Apply rotation to vertices
        rotated_vertices = (R * vertices')';
        
        % Update visualization
        if ishandle(p)
            set(p, 'Vertices', rotated_vertices);
            
            % Update body-fixed axes
            bodyX = R * [bodyAxisLength; 0; 0];
            bodyY = R * [0; bodyAxisLength; 0];
            bodyZ = R * [0; 0; bodyAxisLength];
            
            set(hBodyX, 'UData', bodyX(1), 'VData', bodyX(2), 'WData', bodyX(3));
            set(hBodyY, 'UData', bodyY(1), 'VData', bodyY(2), 'WData', bodyY(3));
            set(hBodyZ, 'UData', bodyZ(1), 'VData', bodyZ(2), 'WData', bodyZ(3));
            
            % Update info display
            current_speed = norm(omega);
            set(hTime, 'String', sprintf('Time: %.2f s\nω = %.2f rad/s', t_video(i), current_speed));
            
            drawnow;
            pause(0.001);
            
            % Capture frame
            frame = getframe(fig);
            writeVideo(v, frame);
        end
    end
    
    % Cleanup
    close(v);
    close(fig);
    fprintf('Animation saved to: %s\n', videoFile);
end

% Quaternion update using RK4
function q_new = update_quaternion(q, omega, dt)
    % RK4 integration steps
    k1 = quaternion_derivative(q, omega);
    k2 = quaternion_derivative(q + 0.5*dt*k1, omega);
    k3 = quaternion_derivative(q + 0.5*dt*k2, omega);
    k4 = quaternion_derivative(q + dt*k3, omega);
    
    % Update quaternion
    q_new = q + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % Normalize to prevent drift
    q_new = q_new / norm(q_new);
end

% Calculate quaternion derivative
function dq = quaternion_derivative(q, omega)
    % Convert angular velocity to quaternion form
    w = [0; omega];
    % Calculate derivative
    dq = 0.5 * quaternion_multiply(q, w);
end

% Quaternion multiplication
function c = quaternion_multiply(a, b)
    % Extract components
    a0 = a(1); a_vec = a(2:4);
    b0 = b(1); b_vec = b(2:4);
    
    % Calculate product
    c0 = a0*b0 - dot(a_vec, b_vec);
    c_vec = a0*b_vec + b0*a_vec + cross(a_vec, b_vec);
    
    c = [c0; c_vec];
end

% Convert quaternion to rotation matrix
function R = quat2rotm(q)
    % Ensure unit quaternion
    q = q / norm(q);
    
    % Extract components
    w = q(1); x = q(2); y = q(3); z = q(4);
    
    % Build rotation matrix
    R = [1-2*y^2-2*z^2,  2*x*y-2*w*z,    2*x*z+2*w*y;
         2*x*y+2*w*z,    1-2*x^2-2*z^2,  2*y*z-2*w*x;
         2*x*z-2*w*y,    2*y*z+2*w*x,    1-2*x^2-2*y^2];
end
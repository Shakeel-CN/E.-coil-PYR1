% Define parameters
D = 200e-12;  % Diffusion coefficient of E. coli (in m^2/s)
v_e = 20e-6;  % Velocity of E. coli (in m/s)
L = 1e-3;  % Length of the channel (in meters)
r = 30e-6;  % Diameter of the channel (in meters)
num_points = 100;  % Number of points in the channel

% Define velocity of water in mm/s
v = 0.54/60;

% Set up grid
dx = L / (num_points - 1);
dy = 2 * r / (num_points - 1);
[x, y] = meshgrid(-L/2:dx:L/2, -r:dy:r);

% Set initial conditions
c = ones(size(x));
c(x < 0) = 0;

% Define the number of iterations and time step
num_iterations = 100; % Number of time steps
dt = 0.5; % Time step in seconds

% Initialize area covered with and without advection
area_no_advection = zeros(num_iterations, 1);
area_with_advection = zeros(num_iterations, 1);

% Initialize speed with and without advection
speed_no_advection = zeros(num_iterations, 1);
speed_with_advection = zeros(num_iterations, 1);

% Create Figure 1 for concentration plot
figure(1)

% Solve advection-diffusion equation
for i = 1:num_iterations
    scaling_factor = (i - 1) / (num_iterations - 1); % Linear scaling factor from 0 to 1
    current_time = (i - 1) * dt; % Calculate current time in seconds
    
    v_term = v * (x / (L/2));
    dc = diff(c, 1, 2);
    dc = [dc, zeros(size(dc, 1), 1)]; % Pad the output of the diff function with zeros
    c = c + (D * (dx^-2) * (dy^-2) * (circshift(c, [0, 1]) + circshift(c, [0, -1]) + circshift(c, [1, 0]) + circshift(c, [-1, 0]) - 4 * c) - (v / v_e) * dc - v_term .* (circshift(c, [0, 1]) - circshift(c, [0, -1])));
    
    % Set boundary conditions
    c(:, 1) = c(:, 2);
    c(:, end) = c(:, end-1);
    
    % Calculate area covered with and without advection
    area_no_advection(i) = sum(sum(c >= 0.5)) * dx * dy;
    area_with_advection(i) = sum(sum((c .* scaling_factor) >= 0.5)) * dx * dy;
    
    % Calculate speed with and without advection
    speed_with_advection(i) = area_no_advection(i) / (pi * r^2 * current_time);
    speed_no_advection(i) = area_with_advection(i) / (pi * r^2 * current_time);
    
    % Plot concentration at current time with scaling factor
    figure(1)
    surf(x, y, c * scaling_factor)
    shading interp
    colormap(jet)
    axis([-L/2 L/2 -r r 0 1])
    xlabel('Distance (meters)')
    ylabel('Radius (meters)')
    zlabel('Concentration')
    title(['Time = ' num2str(current_time) ' s'])
    
    % Plot speed with and without advection as a line graph
    figure(2)
    plot((1:num_iterations) * dt, speed_no_advection, 'b-', 'LineWidth', 2)
    hold on
    plot((1:num_iterations) * dt, speed_with_advection, 'r-', 'LineWidth', 2)
    xlabel('Time (s)')
    ylabel('Speed (micrometers/s)')
    title('Speed with and without advection')
    legend('No advection', 'With advection')
    grid on
    
    drawnow
end

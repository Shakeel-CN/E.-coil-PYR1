% Parameters
D_a = 1.2141e-10; % attractant absolute diffusion coefficient (m^2/s)
Q = 1e-11; % attractant production rate (mol/s)
L = 120e-6; % domain size (m)
x_range = linspace(0, L, 121); % 121 points for 60 um
t_range = linspace(0, 20, 201); % time steps
ecoli_speed = 0.866e-7; % E. coli speed (m/s)
running_time = 1.25; % seconds
tumbling_time = 0.95; % seconds

% Solve PDE to get attractant concentration
sol = pdepe(0, @(x, t, u, DuDx) pdefun(x, t, u, DuDx, Q), @(x) icfun(x, L), @bcfun, x_range, t_range);

% Find maximum concentration along with time
max_concentration = max(sol(:, :, 1), [], 1);
max_concentration_time = t_range(find(max_concentration == max(max_concentration), 1));

% Calculate approximate distance covered and time taken to reach higher concentration
higher_concentration_threshold = max(max_concentration) * 0.8; % Define the threshold
higher_concentration_distances = [];
higher_concentration_times = [];
total_distance_covered = 0; % Initialize total distance covered
for i = 1:length(t_range)
    time = t_range(i);
    if mod(time, tumbling_time + running_time) < tumbling_time
        ecoli_movement = randn(size(x_range)) * ecoli_speed * 0.9; % Tumbling phase
    else
        ecoli_movement = sol(i, :, 1) * ecoli_speed; % Running phase
    end

    higher_concentration_index = find(sol(i, :, 1) >= higher_concentration_threshold, 1);
    if ~isempty(higher_concentration_index)
        higher_concentration_distances = [higher_concentration_distances, x_range(higher_concentration_index) + ecoli_movement(higher_concentration_index) * time];
        higher_concentration_times = [higher_concentration_times, time];
        total_distance_covered = total_distance_covered + ecoli_speed * time; % Accumulate distance covered
    end
end

% Display results
disp(['For Attractant Production Rate: ' num2str(Q) ' mol/s']);
disp(['Maximum Attractant Concentration: ' num2str(max(max_concentration)) ' mol/m^3']);
disp(['Time of Maximum Concentration: ' num2str(max_concentration_time) ' s']);

% Calculate approximate time taken to reach higher concentration from the edge
initial_position = x_range(1); % Position at the edge of the domain
time_to_reach_higher_concentration = (initial_position + higher_concentration_distances(1)) / ecoli_speed - tumbling_time;

if ~isempty(higher_concentration_distances)
    disp(['Absolute Distance Covered by E. coli to Reach Higher Concentration: ' num2str(higher_concentration_distances(1) * 1e6) ' micrometer']);
    disp(['Approximate Distance Covered by E. coli: ' num2str(total_distance_covered) ' m']);
    disp(['Approximate Time Taken by E. coli to Reach Higher Concentration from the Edge: ' num2str(time_to_reach_higher_concentration) ' s']);
else
    disp('Bacteria did not reach higher concentration');
end

% Parameters for Z movement as a sinusoidal function of X and T
amplitude = 5e-6; % 5 micrometers, change as needed
frequency = 2 * pi / (10e-6); % 10 micrometers wavelength, change as needed

% Ecoli movement in X direction over time
ecoli_position = zeros(length(t_range), 1);
z_movement = zeros(length(t_range), 1);

% Initial position of E. coli
ecoli_position(1) = x_range(1);

for i = 2:length(t_range)
    time = t_range(i);
    
    % Calculate concentration gradient
    x_idx = find(x_range <= ecoli_position(i-1), 1, 'last');
    if x_idx == 1
        grad_concentration = (sol(i, x_idx + 1, 1) - sol(i, x_idx, 1)) / (x_range(x_idx + 1) - x_range(x_idx));
    else
        grad_concentration = (sol(i, x_idx, 1) - sol(i, x_idx - 1, 1)) / (x_range(x_idx) - x_range(x_idx - 1));
    end
    
    % Move E. coli based on concentration gradient
    if mod(time, tumbling_time + running_time) < tumbling_time
        ecoli_position(i) = ecoli_position(i-1); % No movement during tumbling phase
    else
        ecoli_position(i) = ecoli_position(i-1) + ecoli_speed * grad_concentration;
    end
    
    % Z movement as sinusoidal function of X and T
    z_movement(i) = amplitude * sin(frequency * ecoli_position(i) + time);
end

% Plot: E. coli movement in 3D (X, T, Z)
figure;
plot3(ecoli_position, t_range, z_movement, '-o');
xlabel('Position (\mum)');
ylabel('Time (s)');
zlabel('Movement (m)');
title('E. coli Movement over Time and Space');
grid on;
colormap(jet);

% Create a new figure for E. coli movement in relation to concentration
figure;

% Initialize E. coli's position
ecoli_position = zeros(length(t_range), 1);

% Initial position of E. coli (e.g., at one end of the domain)
ecoli_position(1) = x_range(1);

% Loop over time to simulate E. coli movement
for i = 2:length(t_range)
    time = t_range(i);
    dt = t_range(i) - t_range(i - 1);
    
    if mod(time, tumbling_time + running_time) < tumbling_time
        % Tumbling phase (random movement)
        ecoli_position(i) = ecoli_position(i - 1) + randn() * ecoli_speed * dt * 0.1;
    else
        % Running phase (movement towards higher concentration)
        concentration_gradient = interp1(x_range, sol(i, :, 1), ecoli_position(i - 1), 'linear', 'extrap');
        ecoli_position(i) = ecoli_position(i - 1) + concentration_gradient * ecoli_speed * dt;
    end
    
    % Keep E. coli within the domain bounds
    ecoli_position(i) = min(max(ecoli_position(i), x_range(1)), x_range(end));
    
    % Get the concentration at E. coli's current position
    current_concentration = interp1(x_range, sol(i, :, 1), ecoli_position(i));
    
    % Plot the E. coli's position at current time, color-coded by concentration
    scatter(ecoli_position(i), time, 50, current_concentration, 'filled'); hold on;
end

% Customize the plot
xlabel('Position (\mum)');
ylabel('Time (s)');
colorbar;
clim([min(sol(:)) max(sol(:))]); % Set color limits to concentration range
colormap('jet');
title('E. coli Movement in Relation to Attractant Concentration');

% Hold off to stop adding to the current plot
hold off;

% Define the icfun function
function u0 = icfun(x, L)
    u0 = zeros(2, 1);
    if abs(x - L/2) <= 20e-6
        u0(1) = 1;
    end
end

% Define the pdefun, icfun, and bcfun functions
function [c, f, s] = pdefun(~, ~, u, DuDx, Q)
    D_a = 1.2141e-10;
    c = [1; 1];
    f = [D_a * DuDx(1); 0];
    s = [-Q; 0];
end

function [pl, ql, pr, qr] = bcfun(~, ~, ~, ~, ~)
    pl = [0; 0];
    ql = [1; 1];
    pr = [0; 0];
    qr = [1; 1];
end
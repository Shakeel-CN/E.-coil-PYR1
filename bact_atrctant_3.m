% Parameters
D_a = 1.2141e-10; % attractant absolute diffusion coefficient (m^2/s)
Q = 1e-11; % attractant production rate (mol/s)
L = 60e-6; % domain size (m)
x_range = linspace(0, L, 121); % 121 points for 60 um
t_range = linspace(0, 10, 101); % time steps
ecoli_speed = 0.866e-6; % E. coli speed (m/s)
running_time = 1.25; % seconds
tumbling_time = 0.75; % seconds

% Solve PDE to get attractant concentration
sol = pdepe(0, @(x, t, u, DuDx) pdefun(x, t, u, DuDx, Q), @(x) icfun(x, L), @bcfun, x_range, t_range);

% Calculate approximate distance covered and time taken to reach higher concentration
higher_concentration_threshold = max(max(sol(:, :, 1))) * 0.8; % Define the threshold
higher_concentration_distances = [];
higher_concentration_times = [];

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
    end
end

% Display results
disp(['For Attractant Production Rate: ' num2str(Q) ' mol/s']);
disp(['Maximum Attractant Concentration: ' num2str(max(max(sol(:, :, 1)))) ' mol/m^3']);

if ~isempty(higher_concentration_distances)
    disp(['Approximate Distance Covered to Reach Higher Concentration: ' num2str(higher_concentration_distances(1) * 1e6) ' µm']);
    disp(['Approximate Time of E. coli to Reach at Maximum Concentration Point: ' num2str(max(t_range)) ' s']);
else
    disp('Bacteria did not reach higher concentration');
end

% Plot 1: Attractant concentration in 3D over time and space
figure;
[x_mesh, t_mesh] = meshgrid(x_range, t_range);
mesh(x_mesh, t_mesh, sol(:, :, 1)); % Removed transpose here
xlabel('Position (\mum)');
ylabel('Time (s)');
zlabel('Attractant Concentration (mol/m^3)');
title('Attractant Concentration over Time and Space');
colormap(jet);
view(3);

% Plot 2: Bacterial movement in 3D over time and space
figure;
ecoli_movement = zeros(length(t_range), length(x_range));
for i = 1:length(t_range)
    time = t_range(i);
    if mod(time, tumbling_time + running_time) < tumbling_time
        ecoli_movement(i, :) = randn(size(x_range)) * ecoli_speed * 0.9; % Tumbling phase
    else
        ecoli_movement(i, :) = sol(i, :, 1) * ecoli_speed; % Running phase
    end
end
[x_mesh, t_mesh] = meshgrid(x_range, t_range);
mesh(x_mesh, t_mesh, ecoli_movement);
xlabel('Position (\mum)');
ylabel('Time (s)');
zlabel('Bacterial Movement (\mum/s)');
title('Bacterial Movement over Time and Space');
colormap(jet);
view(3);

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
function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
    pl = [ul(1); 0]; % Concentration remains the same at left boundary
    ql = [0; 1];
    pr = [ur(1); 0]; % Concentration remains the same at right boundary
    qr = [0; 1];
end



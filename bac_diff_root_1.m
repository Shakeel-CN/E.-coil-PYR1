% Simulation parameters 
diffusion_coefficient = 10e-6;
cell_radius = 1e-5;
surface_area = 2e-6;
growth_rate = 0.39 / 3600;
temperature_range_C = linspace(30, 50, 100);
num_time_steps = 10000;
time_step = 1;

% Initialize data arrays for plotting
time_array = zeros(num_time_steps, 1);
position_array = zeros(num_time_steps, 1);
temperature_array = zeros(num_time_steps, 1);
concentration_array = zeros(num_time_steps, 1);

% Initialize the initial position and concentration
bacteria_position = 0;
bacteria_concentration = 1;

% Simulation loop 
for t = 1:num_time_steps
    bacteria_position = bacteria_position + sqrt(2 * diffusion_coefficient * time_step) * randn();
    bacteria_concentration = bacteria_concentration * exp(growth_rate * time_step);
    time_array(t) = t * time_step;
    position_array(t) = bacteria_position;
    selected_temp_index = randi(length(temperature_range_C));
    selected_temperature_C = temperature_range_C(selected_temp_index);
    temperature_array(t) = selected_temperature_C;
    
    % Update concentration using the diffusion equation
    concentration_array(t) = bacteria_concentration * exp(-bacteria_position^2 / (4 * diffusion_coefficient * t * time_step));
end

% Plot bacterial movement and concentration over time 
figure;
subplot(2, 1, 1);
plot(time_array, position_array);
xlabel('Time (seconds)');
ylabel('Bacteria Position (\mumeters)');
title('(a) Bacterial Diffusion Towards Root Cell');
grid on;

subplot(2, 1, 2);
plot(time_array, concentration_array);
xlabel('Time (seconds)');
ylabel('Bacterial Concentration');
title('(b) Bacterial Concentration Over Time');
grid on;

% Create a 2D scatter plot of temperature vs. position 
figure;
scatter(temperature_array, position_array, 10, 'filled');
xlabel('Temperature (°C)');
ylabel('Bacteria Position (\mumeters)');
title('Temperature Variation (°C) and Bacterial Position');
colormap(jet);
grid on;

% Create a 3D scatter plot of time, temperature, and bacterial position
figure;
scatter3(time_array, temperature_array, position_array, 20, position_array, 'filled');
xlabel('Time (seconds)');
ylabel('Temperature (°C)');
zlabel('Bacterial Position (\mumeters)');
title('Time, Temperature, and Bacterial Position');
colormap(jet);
grid on;

% Display the average and total time
average_time = mean(position_array.^2) / (2 * diffusion_coefficient);
total_time = sum(position_array.^2) / (2 * diffusion_coefficient);

fprintf('Average time taken by bacteria to diffuse into root: %.2f seconds\n', average_time);
fprintf('Total time taken by bacteria to diffuse into root: %.2f seconds\n', total_time);

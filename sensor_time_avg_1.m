% Define parameters with stochasticity and temperature sensitivity
temperature = 25; % Assume room temperature (in °C)
Q10 = 2; % Reaction rate increases by a factor of 2 for every 10°C rise
temperature_factor = Q10^((temperature - 25) / 10);
params = struct( ...
    'k_prod', 0.01 * temperature_factor + randn * 0.01, ... % ABA production rate constant (concentration/time)
    'k_deg', 0.001 * temperature_factor + randn * 0.01, ... % ABA degradation rate constant (1/time)
    'k_on', 10 * temperature_factor + randn * 0.1, ... % ABA binding to receptor (M^-1 min^-1)
    'k_off', 0.1 * temperature_factor + randn * 0.01, ... % ABA dissociation from receptor (min^-1)
    'k_half', 3e-3, ... % Half-maximal concentration for the Hill equation
    'n', 2.5, ... % Cooperativity coefficient in ligand binding
    'k_trans_deg', 0.05 * temperature_factor + randn * 0.01, ... % Degradation rate constant for mRNA (1/time)
    'k_trans_syn', 0.05 * temperature_factor + randn * 0.01, ... % Translation rate constant (concentration/time)
    'k_deg_prot', 0.001 * temperature_factor + randn * 0.01, ... % Degradation rate constant for protein (1/time)
    'k_response', 1e-3 * temperature_factor + randn * 0.01, ... % Rate constant for bioluminescence response (1/time)
    'R_total', 1 ... % Total receptor concentration (concentration)
);

% Define initial conditions
C_ABA_init = 1:2:40; % Initial ABA concentration range (0 to 50 µM)
I_threshold = 0:0.5:5; % Threshold bioluminescence intensity range (0 to 5)
C_ABA_R_init = 1; % Initial concentration of ABA-receptor complex
mRNA_init = 0; % Initial concentration of sensor protein mRNA
Protein_init = 1; % Initial concentration of sensor protein (starting with 1)
I_Bioluminescence_init = 0; % Initial bioluminescence response

% Time span for the simulation
tspan = [0 1000];

% Preallocate results for detection time
detection_time = zeros(length(C_ABA_init), length(I_threshold));

% Solve the ODEs for different initial ABA concentrations and thresholds
for i = 1:length(C_ABA_init)
    for j = 1:length(I_threshold)
        % Initial conditions for the current simulation
        Y0 = [C_ABA_init(i), C_ABA_R_init, mRNA_init, Protein_init, I_Bioluminescence_init];
        
        % Solve ODEs
        [T, Y] = ode15s(@(t, Y) aba_biosensor_odes(t, Y, params), tspan, Y0);
        
        % Find the detection time when bioluminescence exceeds the threshold
        idx = find(Y(:, 5) >= I_threshold(j), 1);
        if ~isempty(idx)
            detection_time(i, j) = T(idx)/60;
        else
            detection_time(i, j) = NaN; % If the threshold is not reached
        end
    end
end

% Calculate average detection time
valid_times = detection_time(~isnan(detection_time));
average_detection_time = mean(valid_times);

% Display the result
fprintf('The average total detection time is %.2f minutes.\n', average_detection_time);

% Plot detection time 
[X, Y] = meshgrid(C_ABA_init, I_threshold);
figure;
surf(X, Y, detection_time');
xlabel('Initial ABA Concentration (µM)');
ylabel('Threshold Bioluminescence Intensity');
zlabel('Detection Time (m)');
title('Detection Time for Different Initial ABA Concentrations and Thresholds');
grid on;
colorbar;

% Define the system of ODEs with all modifications
function dYdt = aba_biosensor_odes(t, Y, params)
    % Unpack parameters
    k_prod = params.k_prod; % ABA production rate (already includes stochasticity)
    k_deg = params.k_deg; % ABA degradation rate (already includes stochasticity)
    k_on = params.k_on; % Association rate (already includes stochasticity)
    k_off = params.k_off; % Dissociation rate (already includes stochasticity)
    k_trans = (0.1 / (1 + Y(4)/10)); % Transcription rate with negative feedback
    k_trans_syn = params.k_trans_syn; % Translation rate (already includes stochasticity)
    k_deg_prot = params.k_deg_prot; % Protein degradation rate (already includes stochasticity)
    k_response = params.k_response; % Bioluminescence response rate (already includes stochasticity)
    
    % Unpack state variables
    C_ABA = Y(1);
    C_ABA_R = Y(2);
    mRNA = Y(3);
    Protein = Y(4);
    I_Bioluminescence = Y(5);

    % Define ODEs with updated dependencies
    dC_ABA_dt = k_prod - k_deg * C_ABA - k_on * C_ABA * (params.R_total - C_ABA_R) + k_off * C_ABA_R;
    dC_ABA_R_dt = k_on * C_ABA * (params.R_total - C_ABA_R) - k_off * C_ABA_R;
    dmRNA_dt = k_trans * (C_ABA^params.n / (params.k_half^params.n + C_ABA^params.n)) - params.k_trans_deg * mRNA;
    dProtein_dt = k_trans_syn * mRNA - k_deg_prot * Protein;
    dI_Bioluminescence_dt = k_response * Protein; % Bioluminescence now depends on protein
    
    % Pack derivatives
    dYdt = [dC_ABA_dt; dC_ABA_R_dt; dmRNA_dt; dProtein_dt; dI_Bioluminescence_dt];
end

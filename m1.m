clc; clear; close all;
data = readmatrix("cmj.xlsx");
g = 9.81;
time_ms = data(:,1); 
time_s = time_ms * 0.001;
time_s = time_s(1261:end,1);
Fz = data(1261:end,4); 
Fz = Fz * g;
BW = 65 * g;
mass = 65;

force_threshold = BW * 0.95;
takeoff_threshold = 5;
landing_threshold = 5;

window_size = 50; 
static_end_idx = 1;
for i = 1:(length(Fz)-window_size)
    if std(Fz(i:i+window_size-1)) < 5 && abs(mean(Fz(i:i+window_size-1)) - BW) < BW*0.1
        static_end_idx = i + window_size - 1;
    else
        break;
    end
end

unloading_start_idx = find(Fz(static_end_idx:end) < force_threshold, 1, 'first') + static_end_idx - 1;
takeoff_idx = find(Fz(unloading_start_idx:end) < takeoff_threshold, 1, 'first') + unloading_start_idx - 1;
takeoff_time = time_s(takeoff_idx);
landing_idx = find(Fz(takeoff_idx:end) > landing_threshold, 1, 'first') + takeoff_idx - 1;
landing_time = time_s(landing_idx);
flight_time = landing_time - takeoff_time;

acceleration = (Fz - BW) / mass;
velocity = cumtrapz(time_s, acceleration);
displacement = cumtrapz(time_s, velocity);

[~, lowest_COM_idx] = min(displacement(unloading_start_idx:takeoff_idx));
lowest_COM_idx = lowest_COM_idx + unloading_start_idx - 1;
breaking_phase = unloading_start_idx:lowest_COM_idx;
propulsive_phase = lowest_COM_idx:takeoff_idx;
landing_phase = landing_idx:length(Fz);

dt_flight = flight_time;
JH_flight = (g * (dt_flight)^2) / 8;
takeoff_V = velocity(takeoff_idx);
impulse_propulsive = trapz(time_s(propulsive_phase), Fz(propulsive_phase) - BW);
JH_impulse = (impulse_propulsive^2) / (2 * mass^2 * g);
peak_V_propulsive = max(velocity(propulsive_phase));
impulse_landing = trapz(time_s(landing_phase), Fz(landing_phase) - BW);
dt_propulsive = time_s(takeoff_idx) - time_s(lowest_COM_idx);
RSI_modified = JH_flight / dt_propulsive;
power = Fz .* velocity;
peak_P_propulsive = max(power(propulsive_phase));
COM_takeoff = displacement(takeoff_idx);
mean_F_propulsive = mean(Fz(propulsive_phase));
peak_F_propulsive = max(Fz(propulsive_phase));
peak_F_breaking = max(Fz(breaking_phase));
mean_P_breaking = mean(power(breaking_phase));
peak_V_negative = min(velocity(breaking_phase));
impulse_breaking = trapz(time_s(breaking_phase), Fz(breaking_phase) - BW);
impulse_unloading = trapz(time_s(unloading_start_idx:lowest_COM_idx), Fz(unloading_start_idx:lowest_COM_idx) - BW);
mean_F_breaking = mean(Fz(breaking_phase));
dt_propulsive = time_s(takeoff_idx) - time_s(lowest_COM_idx);
delta_y_breaking = displacement(lowest_COM_idx) - displacement(unloading_start_idx);
leg_stiffness = peak_F_breaking / delta_y_breaking;
COM_depth = displacement(unloading_start_idx) - displacement(lowest_COM_idx);
dt_jump = time_s(landing_idx) - time_s(unloading_start_idx);
ratio_flight_jump = dt_flight / dt_jump;

results = [dt_flight, JH_flight, takeoff_V, JH_impulse, peak_V_propulsive, impulse_landing,...
    RSI_modified, peak_P_propulsive, COM_takeoff, mean_F_propulsive, peak_F_propulsive,...
    peak_F_breaking, mean_P_breaking, peak_V_negative, impulse_breaking, impulse_unloading,...
    mean_F_breaking, dt_propulsive, leg_stiffness, COM_depth, ratio_flight_jump];

T = array2table(results, 'VariableNames', {'dt_flight', 'JH_flight', 'takeoff_V', 'JH_impulse',...
    'peak_V_propulsive', 'impulse_landing', 'RSI_modified', 'peak_P_propulsive', 'COM_takeoff',...
    'mean_F_propulsive', 'peak_F_propulsive', 'peak_F_breaking', 'mean_P_breaking',...
    'peak_V_negative', 'impulse_breaking', 'impulse_unloading', 'mean_F_breaking',...
    'dt_propulsive', 'leg_stiffness', 'COM_depth', 'ratio_flight_jump'});
writetable(T, 'CMJ_Results.xlsx');

figure;
yyaxis left;
plot(time_s, Fz, 'b', 'LineWidth', 1.5);
ylabel('Force (N)', 'FontWeight', 'bold');
yyaxis right;
plot(time_s, velocity, 'r', 'LineWidth', 1.5);
ylabel('Velocity (m/s)', 'FontWeight', 'bold');
xlabel('Time (s)', 'FontWeight', 'bold');
title('Force-Velocity Profile with Phase Annotation', 'FontSize', 12);
grid on;

xline(time_s(unloading_start_idx), '--k', 'Unloading Start');
xline(time_s(lowest_COM_idx), '--m', 'Breaking End');
xline(time_s(takeoff_idx), '--g', 'Take-off');
xline(time_s(landing_idx), '--c', 'Landing');
legend('Force', 'Velocity', 'Unloading Start', 'Breaking End', 'Take-off', 'Landing',...
    'Location', 'southeast');

phase_labels = {'Unloading', 'Breaking', 'Propulsive', 'Flight', 'Landing'};
phase_colors = [0.9 0.9 0.9; 0.8 0.8 1; 0.7 1 0.7; 1 0.8 0.8; 0.8 0.7 0.9];
x_regions = [time_s(1), time_s(unloading_start_idx);...
    time_s(unloading_start_idx), time_s(lowest_COM_idx);...
    time_s(lowest_COM_idx), time_s(takeoff_idx);...
    time_s(takeoff_idx), time_s(landing_idx);...
    time_s(landing_idx), time_s(end)];


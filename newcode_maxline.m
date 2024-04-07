% Load data
data = readtable('trajectory_26.txt');
N = size(data, 1);

% Parameters
w = 2;
det_t = 0.02;
t_min = w * det_t;
t_max = (N - w) * det_t;

% Time vector
t = t_min:det_t:t_max;
tau_range = (-(N - w) * det_t):det_t:(N * det_t);

% Extracting and computing velocities
time = table2array(data(:, 1));
ax = table2array(data(:, 2));
ay = table2array(data(:, 3));
bx = table2array(data(:, 4));
by = table2array(data(:, 5));

delta_time = diff(time);
vx_a = diff(ax) ./ delta_time;
vy_a = diff(ay) ./ delta_time;
vx_b = diff(bx) ./ delta_time;
vy_b = diff(by) ./ delta_time;

% Magnitude and unit vectors
v_mag_a = sqrt(vx_a.^2 + vy_a.^2);
v_mag_b = sqrt(vx_b.^2 + vy_b.^2);
unit_vx_a = vx_a ./ v_mag_a;
unit_vy_a = vy_a ./ v_mag_a;
unit_vx_b = vx_b ./ v_mag_b;
unit_vy_b = vy_b ./ v_mag_b;

% Initialize the C_ij matrix for the computed values
C_ij = NaN(length(t), length(tau_range));

% Convert the time array to a column vector for broadcasting
time_col = time(:);

for idx_t = 1:length(t)
    current_t = t(idx_t);
    
    for idx_tau = 1:length(tau_range)
        current_tau = tau_range(idx_tau);
        
        sum_v_ij = 0;
        valid_counts = 2*w+1;
        
       
        for k = -w:w
            tk = current_t + k * det_t;
            % Find the indices for a and b at times t and t + tau
            idx_a = find(time_col <= tk, 1, 'last');
            idx_b = find(time_col <= tk + current_tau, 1, 'last');
            
            % Check if both indices are valid and within the range
            if ~isempty(idx_a) && ~isempty(idx_b) && idx_a < length(unit_vx_a) && idx_b < length(unit_vx_b)
                % Compute the dot product for a and b at the given times
                dot_product = unit_vx_a(idx_a) * unit_vx_b(idx_b) + unit_vy_a(idx_a) * unit_vy_b(idx_b);
                sum_v_ij = sum_v_ij + dot_product;
            end
        end
           C_ij(idx_t, idx_tau) = sum_v_ij / valid_counts;
        if C_ij(idx_t,idx_tau) > 0 && C_ij(idx_t,idx_tau) < 0
            C_ij(idx_t, idx_tau) = C_ij(idx_t, idx_tau);
        else
            C_ij(C_ij == 0) = NaN;
        end
    end
end


% Find the maximum delay at each time point
max_delay_per_time = NaN(length(t), 1);
for idx_t = 1:length(t)
    [max_val, idx_max_tau] = max(C_ij(idx_t, :));
    if ~isnan(max_val)
        max_delay_per_time(idx_t) = tau_range(idx_max_tau);
    end
end

% Calculate the average maximum delay
average_max_delay = mean(max_delay_per_time, 'omitnan');

% Plotting
figure;
contourf_handle = contourf(t, tau_range, C_ij.', 20); 
colormap(jet);
colorbar;
hold on;

% Plot the average maximum delay line
line([t(1), t(end)], [average_max_delay, average_max_delay], 'Color', 'white', 'LineStyle', '--', 'LineWidth', 2);

xlabel('Time t (s)'); 
ylabel('Delay \tau (s)');  
title('Delayed correlation C_{ij} with Average Maximum Delay Line'); 

hold off;

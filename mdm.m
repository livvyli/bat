%the latest start at 412
%ihad continue what u are doing. but i dont delete yours, mine until row 138
%find total length of the data
data = readtable('trajectory_21.txt');
sizeOfData = size(data);
N = sizeOfData(1);


% define t and tau
w = 2;
det_t = 0.02;
t_min = w * det_t;
t_max = (N - w)* det_t;
t = t_min: 0.02 : t_max;

% tau_min = -t_min; 
% tau_max = N * det_t;
% %tau = tau_min: 0.02 :tau_max;
% tau=-extime:0.02:tau_max

%define vi and vj
time = table2array(data(:, 1));
ax = table2array(data(:, 2));
ay = table2array(data(:, 3));
bx = table2array(data(:, 4));
by = table2array(data(:, 5));
%difference between consecutive position measurements
delta_ax = diff(ax);
delta_ay = diff(ay);
delta_bx = diff(bx);
delta_by = diff(by);
delta_time = diff(time);
%velocity of a and b
vx_a = delta_ax ./ delta_time;
vy_a = delta_ay ./ delta_time;
vx_b = delta_bx ./ delta_time;
vy_b = delta_by ./ delta_time;

v_mag_a = sqrt(vx_a.^2 + vy_a.^2);
v_mag_b = sqrt(vx_b.^2 + vy_b.^2);
%unit vector for a and b
v_unit_ax = vx_a ./ v_mag_a;
v_unit_ay = vy_a ./ v_mag_a;
v_unit_bx = vx_b ./ v_mag_b;
v_unit_by = vy_b ./ v_mag_b;

% Create lookup tables for velocity components and magnitudes
vx_a_lookup = containers.Map(time(1:end-1), vx_a);
vy_a_lookup = containers.Map(time(1:end-1), vy_a);
vx_b_lookup = containers.Map(time(1:end-1), vx_b);
vy_b_lookup = containers.Map(time(1:end-1), vy_b);
v_mag_a_lookup = containers.Map(time(1:end-1), v_mag_a);
v_mag_b_lookup = containers.Map(time(1:end-1), v_mag_b);

% v_unit_ax_lookup = containers.Map(time(1:end-1), v_unit_ax);
% v_unit_ay_lookup = containers.Map(time(1:end-1), v_unit_ay);
% v_unit_bx_lookup = containers.Map(time(1:end-1), v_unit_bx);
% v_unit_by_lookup = containers.Map(time(1:end-1), v_unit_by);

% define the dot product
total_tau = (-(N - w) * det_t):0.02:(N * det_t);
% Initialize dot_products with NaNs
dot_products = cell(length(t), 1);


for idx_t = 1:length(t)
    tau_min = -(t(idx_t));
    tau_max = N * det_t;
    tau = tau_min:0.02:tau_max;

    % Initialize the cell array for the current row
    dot_products{idx_t} = NaN(1, length(tau));

    for idx_tau = 1:length(tau)
        % Look up velocity components and magnitudes based on time
        if isKey(vx_a_lookup, t(idx_t)) && isKey(vx_b_lookup, t(idx_t) + tau(idx_tau))
            % Get the components of the vector at time idx_t
            v1_x = vx_a_lookup(t(idx_t));
            v1_y = vy_a_lookup(t(idx_t));
            % Normalize the vector
            v1_mag = sqrt(v1_x^2 + v1_y^2);
            v1_x_normalized = v1_x / v1_mag;
            v1_y_normalized = v1_y / v1_mag;
            
            % Get the components of the vector at time idx_t + idx_tau
            v2_x = vx_b_lookup(t(idx_t) + tau(idx_tau));
            v2_y = vy_b_lookup(t(idx_t) + tau(idx_tau));
            % Normalize the vector
            v2_mag = sqrt(v2_x^2 + v2_y^2);
            v2_x_normalized = v2_x / v2_mag;
            v2_y_normalized = v2_y / v2_mag;
            
            % Calculate dot product
            v_ab = v1_x_normalized * v2_x_normalized + v1_y_normalized * v2_y_normalized;

            % Store v_ab in dot_products at the correct index
            dot_products{idx_t}(idx_tau) = v_ab;
        else
            % If one of the keys is missing, store NaN
            dot_products{idx_t}(idx_tau) = NaN;
        end
    end
end

% Combine values starting from the right to the left in each row
for idx_t = 1:length(t)
    % Get the current row
    row_values = dot_products{idx_t};
    
    % Find the last non-NaN index in the row
    last_non_nan_idx = find(~isnan(row_values), 1, 'last');
    
    % Iterate from the last non-NaN index towards the first index
    for idx_tau = last_non_nan_idx-1:-1:1
        % Check if the current value and the previous value are not NaN
        if ~isnan(row_values(idx_tau)) && ~isnan(row_values(idx_tau + 1))
            % Combine the current value with the previous value
            row_values(idx_tau) = row_values(idx_tau) * row_values(idx_tau + 1);
        end
    end
    
    % Update the row in dot_products with the combined values
    dot_products{idx_t} = row_values;
end

% Determine the maximum number of columns needed for any row
max_columns = max(cellfun(@numel, dot_products));

% Create a matrix with NaN values
combined_matrix = NaN(length(t), max_columns);

% Iterate over each row of dot_products
for idx_t = 1:length(t)
    % Determine the starting column index for the current row
    start_column = max_columns - numel(dot_products{idx_t}) + 1;
    
    % Fill in the values of dot_products into the appropriate positions
    combined_matrix(idx_t, start_column:end) = dot_products{idx_t};
end


%below is what u code
%find total length of the data
data = readtable('trajectory_21.txt');
sizeOfData = size(data);
N = sizeOfData(1);

% define t and tau
w = 2;
det_t = 0.02;
t_min = w * det_t;
t_max = (N - w)* det_t;
t = t_min: 0.02 : t_max;

tau_min = -t_max; 
tau_max = N * det_t;
tau = tau_min: 0.02 :tau_max;

%define vi and vj
time = table2array(data(:, 1));
ax = table2array(data(:, 2));
ay = table2array(data(:, 3));
bx = table2array(data(:, 4));
by = table2array(data(:, 5));
%difference between consecutive position measurements
delta_ax = diff(ax);
delta_ay = diff(ay);
delta_bx = diff(bx);
delta_by = diff(by);
delta_time = diff(time);
%velocity of a and b
vx_a = delta_ax ./ delta_time;
vy_a = delta_ay ./ delta_time;
vx_b = delta_bx ./ delta_time;
vy_b = delta_by ./ delta_time;

v_mag_a = sqrt(vx_a.^2 + vy_a.^2);
v_mag_b = sqrt(vx_b.^2 + vy_b.^2);
%unit vector for a and b
v_unit_ax = vx_a ./ v_mag_a;
v_unit_ay = vy_a ./ v_mag_a;
v_unit_bx = vx_b ./ v_mag_b;
v_unit_by = vy_b ./ v_mag_b;

%define the dot product
%%not sure why t_tau_index is idx_t + w + idex_tau -1

dot_products = zeros(length(time)-1, length(tau));
for idx_t = 1:length(time)-1
    for idx_tau = 1:length(tau)
        t_tau_index = idx_t + w + idx_tau - 1;
        %t_tau_index = idx_t + idx_tau - 1;
        if t_tau_index <= length(time)-1
            dot_products(idx_t, idx_tau) = ...
                v_unit_ax(idx_t) * v_unit_bx(t_tau_index) + ... 
                v_unit_ay(idx_t) * v_unit_by(t_tau_index);
        end
    end
end


%%%% second edition
% I don't know whether this is correct or not, but I already work on it for five hours :(
%hope this is helpful... still not sure about this contour graph, if you have any progress feel free to message me, Thanks!

% Load the data from the file
data = readtable('trajectory_24.txt');
N = height(data);

% Define parameters
w = 2;                 
det_t = 0.02;           
t_k = (0:N-1) * det_t;   

t_min = w * det_t;
t_max = (N - w) * det_t;
t = (t_min:det_t:t_max);
% Define the tau vector
tau_min = -w * det_t;
tau_max = (N-w-1)  * det_t;
tau = (tau_min:det_t:tau_max)'; 

% Initialize variables for positions and velocities
ax = data{:, 2};
ay = data{:, 3};
bx = data{:, 4};
by = data{:, 5};

% Calculate velocity
delta_ax = diff(ax);
delta_ay = diff(ay);
delta_bx = diff(bx);
delta_by = diff(by);

% Calculate velocity magnitudes and unit vectors
v_mag_a = sqrt(delta_ax.^2 + delta_ay.^2);
v_mag_b = sqrt(delta_bx.^2 + delta_by.^2);
v_unit_ax = delta_ax ./ v_mag_a;
v_unit_ay = delta_ay ./ v_mag_a;
v_unit_bx = delta_bx ./ v_mag_b;
v_unit_by = delta_by ./ v_mag_b;

% Pad velocities with zeros to maintain the array size
v_unit_ax = [v_unit_ax; 0];
v_unit_ay = [v_unit_ay; 0];
v_unit_bx = [v_unit_bx; 0];
v_unit_by = [v_unit_by; 0];

% Calculate the dot products for each t and tau
C_ij = NaN(length(tau), length(t)); 

% Corrected loop for calculating the dot products
for idx_t = 1:length(t)
    for idx_tau = 1:length(tau)
        sum_vij = 0;
        valid_count = 2*w+1;

        for k = -w:w
        %don't know why everytime I try  t_idx = idx_t + t_k is not work...
            t_idx = idx_t + k;
            tau_idx = idx_tau; 

            % Ensure indices are within bounds and integers
            if t_idx > 0 && t_idx <= length(t) && tau_idx > 0 && tau_idx <= length(tau)
                dot_product = v_unit_ax(t_idx) * v_unit_bx(tau_idx) + ...
                              v_unit_ay(t_idx) * v_unit_by(tau_idx);
                sum_vij = sum_vij + dot_product;
            end
        end

        if valid_count > 0
            C_ij(idx_tau, idx_t) = sum_vij / valid_count; 
        end
    end
end

disp(size(t));    % [1, length(t)]
disp(size(tau));  % [length(tau), 1]
disp(size(C_ij)); % [length(tau), length(t)]

%relative position
data.rx = bx - ax;
data.ry = by - ay;
data.relative_distance = sqrt((data.rx).^2 + (data.ry).^2);
%heading angle
heading_angle_bat_a = atan2d([NaN; diff(ay)], [NaN; diff(ax)]);
heading_angle_bat_b = atan2d([NaN; diff(by)], [NaN; diff(bx)]);
data.relative_heading_angles = abs(heading_angle_bat_b - heading_angle_bat_a);

% Define a suitable angle and distance threshold for chase 
chase_distance_threshold = 5; 
chase_angle_threshold = 20; 
% Define a suitable angle and distance threshold for coordinated flight
coordinated_distance_threshold = 10; 
coordinated_angle_difference_threshold = 90; 
% Initialize interaction type column
data.interaction_type = repmat("no interaction", height(data), 1);
% Identify chases
chase_indices = data.relative_distance < chase_distance_threshold & data.relative_heading_angles < chase_angle_threshold;
data.interaction_type(chase_indices) = "chase";
% Identify coordinated flights
coordinated_indices = data.relative_distance < coordinated_distance_threshold & abs(data.relative_heading_angles) < coordinated_angle_difference_threshold & ~chase_indices; % To ensure no overlap with chases
data.interaction_type(coordinated_indices) = "coordinated flight";
% Manually summarize interaction types
interaction_types = unique(data.interaction_type);
summary_counts = zeros(size(interaction_types));

for i = 1:length(interaction_types)
    summary_counts(i) = sum(data.interaction_type == interaction_types(i));
end

% Create a table for the summary
summary_table = table(interaction_types, summary_counts, 'VariableNames', {'InteractionType', 'Count'});
% Display the summary table
disp(summary_table);

% figure;
% contourf(t, tau, C_ij, contour_levels, 'LineStyle', 'none'); % Filled contour plot
% colormap('jet'); % Use jet colormap for better visibility
% colorbar; % Add a colorbar
% 
% hold on;
% % [~, h_contour] = contour(t, tau, C_ij', contour_levels, 'LineColor', 'k');
% % clabel(h_contour); % Label contour lines
% 
% % Add labels and title
% xlabel('Time t (s)');
% ylabel('Delay τ (s)');
% title('Delayed Correlation C_{ij}(t, τ)');
% 
% hold off; 

% figure;
% 
% Trajectory plot
subplot(1, 2, 1); % Divide the figure window and select the first panel
hold on;
plot(ax, ay, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', 'Bat A');
plot(bx, by, 'r-x', 'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', 'Bat B');
xlabel('X Position (m)');
ylabel('Y Position (m)');
title('Movement Paths of Paired Bats');
legend('show');
grid on;
hold off;
% 
% % Assuming C_ij, t, and tau have been calculated as in your contour plot script
% % Contour plot
% subplot(1, 2, 2); % Select the second panel
% % Add your contour plotting code here. Example:
% contourf(t, tau, C_ij, 20, 'LineStyle', 'none'); % Adjust parameters as needed
% colorbar;
% xlabel('Time t (s)');
% ylabel('Delay τ (s)');
% title('Delayed Correlation C_{ij}(t, τ)');
% 
% % Adjust the layout for better visibility
% sgtitle('Bat Trajectories and Temporal Correlation Analysis'); % Super title for the figure

% Plot trajectories and contour plot side by side
% figure;
% 
% % Trajectory plot
% subplot(1, 2, 1); % Divide the figure window and select the first panel
% hold on;
% % Plot trajectories for bat A and bat B with time markers
% plot(ax, ay, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', 'Bat A');
% plot(bx, by, 'r-x', 'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', 'Bat B');
% 
% % Annotate specific time points on the trajectories, if necessary
% % Example: plot(ax(10), ay(10), 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
% 
% xlabel('X Position (m)');
% ylabel('Y Position (m)');
% title('Movement Paths of Paired Bats');
% legend('show');
% grid on;
% hold off;
% Load your data and preprocess it as you've done before.
% Assume C_ij is already calculated as the delayed directional correlation matrix
% For each time point, find the delay that maximizes the correlation



% Contour plot
subplot(1, 2, 2); % Select the second panel
[~, h_contourf] = contourf(t, tau, C_ij, 50, 'LineStyle', 'none'); % Adjust the number of levels as needed
colorbar;
hold on;
% Overlay contour lines at a specific correlation threshold, e.g., 0.95
%[~, h_contour] = contour(t, tau, C_ij', [0.95, 0.95], 'LineColor', 'k');
%clabel(h_contour);

xlabel('Time t (s)');
ylabel('Delay τ (s)');
title('Delayed Correlation C_{ij}(t, τ)');
hold off;

% Adjust the layout for better visibility
sgtitle('Bat Trajectories and Temporal Correlation Analysis'); % Super title for the figure

% Example of separate figure for the histogram
figure;
histogram(data.relative_distance, 'BinWidth', 0.1, 'FaceColor', 'b', 'DisplayName', 'Coord - Flight');
hold on;
histogram(data.relative_heading_angles, 'BinWidth', 1, 'FaceColor', 'r', 'DisplayName', 'Chase - Flight');
xlabel('Delay τ (sec)');
ylabel('Count');
title('Histogram of Interaction Types');
legend('show');
hold off;

%ainnewcode
%find total length of the data
data = readtable('trajectory_2.txt');
sizeOfData = size(data);
N = sizeOfData(1);


% define t and tau
w = 2;
det_t = 0.02;
t_min = w * det_t;
t_max = (N - w)* det_t;
t = t_min: 0.02 : t_max;

% tau_min = -t_min; 
% tau_max = N * det_t;
% %tau = tau_min: 0.02 :tau_max;
% tau=-extime:0.02:tau_max

%define vi and vj
time = table2array(data(:, 1));
ax = table2array(data(:, 2));
ay = table2array(data(:, 3));
bx = table2array(data(:, 4));
by = table2array(data(:, 5));
%difference between consecutive position measurements
delta_ax = diff(ax);
delta_ay = diff(ay);
delta_bx = diff(bx);
delta_by = diff(by);
delta_time = diff(time);
%velocity of a and b
vx_a = delta_ax ./ delta_time;
vy_a = delta_ay ./ delta_time;
vx_b = delta_bx ./ delta_time;
vy_b = delta_by ./ delta_time;

v_mag_a = sqrt(vx_a.^2 + vy_a.^2);
v_mag_b = sqrt(vx_b.^2 + vy_b.^2);
%unit vector for a and b
v_unit_ax = vx_a ./ v_mag_a;
v_unit_ay = vy_a ./ v_mag_a;
v_unit_bx = vx_b ./ v_mag_b;
v_unit_by = vy_b ./ v_mag_b;

% Create lookup tables for velocity components and magnitudes
vx_a_lookup = containers.Map(time(1:end-1), vx_a);
vy_a_lookup = containers.Map(time(1:end-1), vy_a);
vx_b_lookup = containers.Map(time(1:end-1), vx_b);
vy_b_lookup = containers.Map(time(1:end-1), vy_b);
v_mag_a_lookup = containers.Map(time(1:end-1), v_mag_a);
v_mag_b_lookup = containers.Map(time(1:end-1), v_mag_b);



% define the dot product
total_tau = (-(N - w) * det_t):0.02:(N * det_t);
% Initialize dot_products with NaNs
dot_products = cell(length(t), 1);


for idx_t = 1:length(t)
    tau_min = -(t(idx_t));
    tau_max = N * det_t;
    tau = tau_min:0.02:tau_max;

    % Initialize the cell array for the current row
    dot_products{idx_t} = NaN(1, length(tau));

    for idx_tau = 1:length(tau)
        % Look up velocity components and magnitudes based on time
        if isKey(vx_a_lookup, t(idx_t)) && isKey(vx_b_lookup, t(idx_t) + tau(idx_tau))
            % Get the components of the vector at time idx_t
            v1_x = vx_a_lookup(t(idx_t));
            v1_y = vy_a_lookup(t(idx_t));
            % Normalize the vector
            v1_mag = sqrt(v1_x^2 + v1_y^2);
            v1_x_normalized = v1_x / v1_mag;
            v1_y_normalized = v1_y / v1_mag;
            
            % Get the components of the vector at time idx_t + idx_tau
            v2_x = vx_b_lookup(t(idx_t) + tau(idx_tau));
            v2_y = vy_b_lookup(t(idx_t) + tau(idx_tau));
            % Normalize the vector
            v2_mag = sqrt(v2_x^2 + v2_y^2);
            v2_x_normalized = v2_x / v2_mag;
            v2_y_normalized = v2_y / v2_mag;
            
            % Calculate dot product
            v_ab = v1_x_normalized * v2_x_normalized + v1_y_normalized * v2_y_normalized;

            % Store v_ab in dot_products at the correct index
            dot_products{idx_t}(idx_tau) = v_ab;
        else
            % If one of the keys is missing, store NaN
            dot_products{idx_t}(idx_tau) = NaN;
        end
    end
end

% Combine values starting from the right to the left in each row
for idx_t = 1:length(t)
    % Get the current row
    row_values = dot_products{idx_t};
    
    % Find the last non-NaN index in the row
    last_non_nan_idx = find(~isnan(row_values), 1, 'last');
    
    % Iterate from the last non-NaN index towards the first index
    for idx_tau = last_non_nan_idx-1:-1:1
        % Check if the current value and the previous value are not NaN
        if ~isnan(row_values(idx_tau)) && ~isnan(row_values(idx_tau + 1))
            % Combine the current value with the previous value
            row_values(idx_tau) = row_values(idx_tau) * row_values(idx_tau + 1);
        end
    end
    
    % Update the row in dot_products with the combined values
    dot_products{idx_t} = row_values;
end

% Determine the maximum number of columns needed for any row
max_columns = max(cellfun(@numel, dot_products));

% Create a matrix with NaN values
combined_matrix = NaN(length(t), max_columns);

% Iterate over each row of dot_products
for idx_t = 1:length(t)
     % Get the current row values
    row_values = dot_products{idx_t};
    % Determine the starting column index for the current row
    start_column = max_columns - numel(dot_products{idx_t}) + 1;
    
    % Fill in the values of dot_products into the appropriate positions
    combined_matrix(idx_t, start_column:end) = dot_products{idx_t};


end

tau_max = N * det_t;
rownames=t_min: 0.02 : t_max;
columnnames=-t_max:0.02:tau_max;
% Convert row and column names to cell arrays of strings
rownames_cell = arrayfun(@num2str, rownames, 'UniformOutput', false);
columnnames_cell = arrayfun(@num2str, columnnames, 'UniformOutput', false);

% Convert combined_matrix to a table
T = array2table(combined_matrix, 'VariableNames', columnnames_cell, 'RowNames', rownames_cell);

%next

% end


% Initialize tdcc as a table with the same size as T
tdcc = array2table(NaN(size(T)), 'VariableNames',columnnames_cell, 'RowNames', rownames_cell);


for ttime=t_min: 0.02 : t_max
    for ttau=-t_max:0.02:tau_max
        summ=[];
        for k=-2:0.02:2
            det_t=0.02;
            t_k=k.*det_t;
            x_pos=ttime+t_k;
            y_pos=ttau;
            % Determine the row index in the table T corresponding to the desired row
            row_idx = find(strcmp(rownames_cell, num2str(x_pos)));
            % Determine the column index in the table T corresponding to the desired column
            col_idx = find(strcmp(columnnames_cell, num2str(y_pos)));
            % Access the value at the specified row and column
            v_ab=T{row_idx, col_idx};

            if ~isnan(v_ab)
                summ = [summ, v_ab];

            end

        end
        %Calculate c_ij only if summ is not empty
        if ~isempty(summ)
            c_ij = 1 / (2 * w + 1) * sum(summ);
        else
            c_ij = NaN;
        end
        
      %tdcc(ttime,ttau)=c_ij
        tdcc{rownames_cell(row_idx), columnnames_cell(col_idx)} = c_ij;
       
    end
end

data.rx = bx - ax;
data.ry = by - ay;
data.relative_distance = sqrt((data.rx).^2 + (data.ry).^2);
%heading angle
heading_angle_bat_a = atan2d([NaN; diff(ay)], [NaN; diff(ax)]);
heading_angle_bat_b = atan2d([NaN; diff(by)], [NaN; diff(bx)]);
data.relative_heading_angles = abs(heading_angle_bat_b - heading_angle_bat_a);

% Define a suitable angle and distance threshold for chase 
chase_distance_threshold = 5; 
chase_angle_threshold = 20; 
% Define a suitable angle and distance threshold for coordinated flight
coordinated_distance_threshold = 10; 
coordinated_angle_difference_threshold = 90; 
% Initialize interaction type column
data.interaction_type = repmat("no interaction", height(data), 1);
% Identify chases
chase_indices = data.relative_distance < chase_distance_threshold & data.relative_heading_angles < chase_angle_threshold;
data.interaction_type(chase_indices) = "chase";
% Identify coordinated flights
coordinated_indices = data.relative_distance < coordinated_distance_threshold & abs(data.relative_heading_angles) < coordinated_angle_difference_threshold & ~chase_indices; % To ensure no overlap with chases
data.interaction_type(coordinated_indices) = "coordinated flight";
% Manually summarize interaction types
interaction_types = unique(data.interaction_type);
summary_counts = zeros(size(interaction_types));

for i = 1:length(interaction_types)
    summary_counts(i) = sum(data.interaction_type == interaction_types(i));
end

% Create a table for the summary
summary_table = table(interaction_types, summary_counts, 'VariableNames', {'InteractionType', 'Count'});
% Display the summary table
disp(summary_table);

rownames=t_min: 0.02 : t_max;
columnnames=-t_max:0.02:tau_max;

% Convert tdcc to a numeric matrix
tdcc_matrix = table2array(tdcc);

% Replace NaN values with a value within the range of your data
tdcc_no_nan = fillmissing(tdcc_matrix, 'constant', 0);

% Transpose tdcc_no_nan to match the dimensions expected by heatmap
tdcc_no_nan_transposed = tdcc_matrix';
  
% Plot the heatmap with custom colormap
figure;
contourf(rownames, columnnames, tdcc_no_nan_transposed);
colormap(jet)
%'Colormap', jet
%'ColorMethod', 'none
xlabel('t');
ylabel('tau');
title('Heatmap of tdcc');
colorbar;

% Example of separate figure for the histogram
figure;
histogram(data.relative_distance, 'BinWidth', 0.1, 'FaceColor', 'b', 'DisplayName', 'Coord - Flight');
hold on;
histogram(data.relative_heading_angles, 'BinWidth', 1, 'FaceColor', 'r', 'DisplayName', 'Chase - Flight');
xlabel('Delay τ (sec)');
ylabel('Count');
title('Histogram of Interaction Types');
legend('show');
hold off;

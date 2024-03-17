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



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



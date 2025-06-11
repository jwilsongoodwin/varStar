clear; clc; close all;

%% === Load Light Curve ===
[file, path] = uigetfile({'*.csv;*.txt;*.tbl', 'Light Curve Files'}, 'Select a Variable Star Light Curve');
if isequal(file, 0), return; end
filename = fullfile(path, file);
disp(['Loading: ', filename]);

[~, ~, ext] = fileparts(filename);
raw = readmatrix(filename);

% === Auto-detect time and magnitude columns ===
if size(raw,2) >= 6 && max(raw(:,6)) > 0
    time = raw(:,2);    % Specific format case
    mag = raw(:,6);
else
    time = raw(:,1);    % General format
    mag = raw(:,2);     % Magnitude
end

% === Clean and Deduplicate ===
valid = isfinite(time) & isfinite(mag);
time = time(valid); mag = mag(valid);
[time, idx] = sort(time); mag = mag(idx);
[time, unique_idx] = unique(time, 'stable');
mag = mag(unique_idx);

% === Clean and Deduplicate ===
valid = isfinite(time) & isfinite(mag);
time = time(valid); mag = mag(valid);
[time, idx] = sort(time); mag = mag(idx);
[time, unique_idx] = unique(time, 'stable');
mag = mag(unique_idx);

% === Convert JD to days since first observation
fprintf('Original time range: %.2f to %.2f\n', min(time), max(time));

% Step 1: Convert JD to MJD (if needed)
if min(time) > 2400000
    time = time - 2400000.5;
    fprintf('Subtracted JD epoch: Now in MJD.\n');
end

% Step 2: Shift to zero
time = time - min(time);

% Debug print
fprintf('Adjusted time range: %.2f to %.2f days\n', min(time), max(time));
%% === Interpolate to Uniform Time Grid ===
dt = median(diff(time));
t_uniform = time(1):dt:time(end);
mag_uniform_raw = interp1(time, mag, t_uniform, 'linear');   % For plotting and phase-folding
mag_uniform = detrend(mag_uniform_raw);                      % For FFT only

%% === FFT Analysis ===
N = length(mag_uniform);
Fs = 1 / dt;
F = fft(mag_uniform);
f = (0:N-1)*(Fs/N);
power = abs(F).^2;
periods = 1 ./ f;

valid = periods > 0.2 & periods < min(300, (time(end)-time(1))/2);
periods = periods(valid);
power = power(valid);

[~, idx_max] = max(power);
dominant_period = periods(idx_max);
fprintf('Detected dominant period: %.4f days\n', dominant_period);

%% === Plot Original Light Curve ===
figure;
plot(time, mag, 'k.');
xlabel('Time (days)'); ylabel('Magnitude');
title('Original Light Curve');
set(gca, 'YDir', 'reverse');  % Brighter magnitudes higher

%% === Plot FFT Power Spectrum ===
figure;
semilogx(periods, power, 'b');
xlabel('Period (days)'); ylabel('Power');
title('FFT Power Spectrum');
xline(dominant_period, 'r--', sprintf(' %.4f d', dominant_period));
grid on;

%% === Zero-Padded FFT and User Selection ===
F_padded = fft([mag_uniform, zeros(1, 5*N)]);
f_padded = (0:length(F_padded)-1)*(Fs/length(F_padded));
periods_padded = 1 ./ f_padded;
power_padded = abs(F_padded).^2;

valid_padded = periods_padded > 0.01 & periods_padded < 0.3;
fprintf('Click on the peak period in the zero-padded FFT plot...\n');
fprintf('Or press Enter without clicking to skip Delta Scuti refinement.\n');

figure;
plot(periods_padded(valid_padded), power_padded(valid_padded), 'b');
xlabel('Period (days)'); ylabel('Power');
title('Zero-Padded FFT – Click to Select Peak Period');
grid on;

[x_click, ~, button] = ginput(1);
isDeltaScuti = false;
if isempty(x_click) || (exist('button', 'var') && button == 13)
    fprintf('No peak selected. Continuing with automatic classification.\n');
    used_period = dominant_period;
else
    periods_in_range = periods_padded(valid_padded);
    [~, idx_nearest] = min(abs(periods_in_range - x_click));
    refined_period = periods_in_range(idx_nearest);
    fprintf('You selected: %.5f days\n', x_click);
    fprintf('Snapped to FFT peak at: %.5f days\n', refined_period);
    use_it = input('Use this period for refined phase-folded light curve? (y/n): ', 's');
    isDeltaScuti = strcmpi(use_it, 'y');
    used_period = isDeltaScuti * refined_period + ~isDeltaScuti * dominant_period;
end

%% === Update Classification ===
if used_period < 0.2
    class = 'Delta Scuti';
elseif used_period < 1.0
    class = 'RR Lyrae';
elseif used_period < 100
    class = 'Cepheid';
else
    class = 'Mira';
end
fprintf('Updated classification: %s\n', class);

%% === Phase-Folded Light Curve ===
phase = mod(t_uniform, used_period) / used_period;
phase(phase > 0.5) = phase(phase > 0.5) - 1;
[phase_sorted, idx] = sort(phase);
mag_sorted = mag_uniform_raw(idx);  % Use raw magnitudes for folding

nbins = 100;
edges = linspace(-0.5, 0.5, nbins+1);  % Corrected phase range
centers = (edges(1:end-1) + edges(2:end)) / 2;
mean_mag = zeros(1, nbins);

for i = 1:nbins
    in_bin = phase_sorted >= edges(i) & phase_sorted < edges(i+1);
    if any(in_bin)
        mean_mag(i) = mean(mag_sorted(in_bin));
    else
        mean_mag(i) = NaN;
    end
end

figure;
plot(centers, mean_mag, '+k');
xlabel('Phase'); ylabel('Mean Magnitude');
title(sprintf('Refined Phase-Folded Light Curve\nP = %.5f d (%s)', used_period, class));
set(gca, 'YDir', 'reverse');  % Magnitude: smaller = brighter
xlim([-0.5, 0.5]);             % Keep it consistent
grid on;

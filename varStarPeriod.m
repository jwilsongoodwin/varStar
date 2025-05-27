clear; clc; close all;

%% === Load Light Curve ===
[file, path] = uigetfile({'*.csv;*.txt;*.tbl', 'Light Curve Files'}, 'Select a Variable Star Light Curve');
if isequal(file, 0), return; end
filename = fullfile(path, file);
disp(['Loading: ', filename]);

[~, ~, ext] = fileparts(filename);
raw = readmatrix(filename);

% === Auto-detect time and flux columns ===
if size(raw,2) >= 6 && max(raw(:,6)) > 0
    time = raw(:,2);    % For V445 Lyr
    flux = raw(:,6);
else
    time = raw(:,1);    % General case
    flux = raw(:,2);
end

% === Clean and Deduplicate ===
valid = isfinite(time) & isfinite(flux);
time = time(valid); flux = flux(valid);
[time, idx] = sort(time); flux = flux(idx);
[time, unique_idx] = unique(time, 'stable');
flux = flux(unique_idx);

if numel(time) < 100
    error('Too few valid time points.');
end

%% === Interpolate to Uniform Time Grid ===
dt = median(diff(time));
t_uniform = time(1):dt:time(end);
flux_uniform = interp1(time, flux, t_uniform, 'linear');
flux_uniform = detrend(flux_uniform);

%% === FFT Analysis ===
N = length(flux_uniform);
Fs = 1 / dt;
F = fft(flux_uniform);
f = (0:N-1)*(Fs/N);
power = abs(F).^2;
periods = 1 ./ f;

% Valid period range for variable stars
valid = periods > 0.2 & periods < min(300, (time(end)-time(1))/2);
periods = periods(valid);
power = power(valid);

% Detect dominant period
[~, idx_max] = max(power);
dominant_period = periods(idx_max);
fprintf('Detected dominant period: %.4f days\n', dominant_period);

%% === Classification by Period Range (initial guess only) ===
if dominant_period < 0.2
    class = 'Delta Scuti (auto)';
elseif dominant_period < 1.0
    class = 'RR Lyrae';
elseif dominant_period < 100
    class = 'Cepheid';
else
    class = 'Mira';
end
fprintf('Auto-classified as: %s\n', class);

%% === Plot Original Light Curve ===
figure;
plot(time, flux, 'k.');
xlabel('Time (days)');
ylabel('Flux');
title('Original Light Curve');

%% === Plot FFT Power Spectrum ===
figure;
semilogx(periods, power, 'b');
xlabel('Period (days)');
ylabel('Power');
title('FFT Power Spectrum');
xline(dominant_period, 'r--', sprintf(' %.4f d', dominant_period));
grid on;

%% === Zero-Padded FFT and User Selection ===
F_padded = fft([flux_uniform, zeros(1, 5*N)]);
f_padded = (0:length(F_padded)-1)*(Fs/length(F_padded));
periods_padded = 1 ./ f_padded;
power_padded = abs(F_padded).^2;

% Limit to high-frequency region (Delta Scuti zone)
valid_padded = periods_padded > 0.01 & periods_padded < 0.3;
fprintf('Click on the peak period in the zero-padded FFT plot...\n');
fprintf('Or press Enter without clicking to skip Delta Scuti refinement.\n');

figure;
plot(periods_padded(valid_padded), power_padded(valid_padded), 'b');
xlabel('Period (days)'); ylabel('Power');
title('Zero-Padded FFT – Click to Select Peak Period');
grid on;

[x_click, ~, button] = ginput(1);  % Wait for click or Enter

if isempty(x_click) || (exist('button', 'var') && button == 13)
    fprintf('No peak selected. Continuing with automatic classification.\n');
    isDeltaScuti = false;
else
    refined_period = x_click;
    fprintf('You selected period: %.5f days\n', refined_period);
    use_it = input('Use this period to create refined phase-folded plot? (y/n): ', 's');
    isDeltaScuti = strcmpi(use_it, 'y');
end


%% === Phase-Folded Light Curve Using FFT Period ===
phase = mod(t_uniform, dominant_period) / dominant_period;
[phase_sorted, idx] = sort(phase);
flux_sorted = flux_uniform(idx);

% Downsample folded curve
figure;
step = max(1, floor(length(phase_sorted)/2000));
plot(phase_sorted(1:step:end), flux_sorted(1:step:end), '.k');
xlabel('Phase'); ylabel('Flux');
title(sprintf('Phase-Folded Light Curve\nP = %.4f d (%s)', dominant_period, class));
grid on;

% Smoothed phase-folded curve
window_size = 100;
smoothed_flux = movmean(flux_sorted, window_size);
figure;
plot(phase_sorted, smoothed_flux, '-k');
xlabel('Phase'); ylabel('Smoothed Flux');
title(sprintf('Smoothed Phase-Folded Light Curve\nP = %.4f d (%s)', dominant_period, class));
grid on;

%% === If Confirmed, Use Clicked Period for Final Plot ===
if isDeltaScuti
    fprintf('\n=== Refining With Selected Period ===\n');
    
    % Phase-fold using selected period
    phase_final = mod(t_uniform, refined_period) / refined_period;
    [phase_sorted_final, idx_final] = sort(phase_final);
    flux_final_sorted = flux_uniform(idx_final);

    % Smooth the full sorted flux
    window_size = 100;
    smoothed_final = movmean(flux_final_sorted, window_size);

    % Downsample both phase and smoothed flux for plotting
    step = max(1, floor(length(phase_sorted_final)/2000));
    phase_ds = phase_sorted_final(1:step:end);
    flux_ds = smoothed_final(1:step:end);

    % Plot clean, smoothed, downsampled curve
    figure;
    plot(phase_ds, flux_ds, '-k', 'LineWidth', 1.5);
    xlabel('Phase');
    ylabel('Smoothed Flux');
    title(sprintf('Refined Phase-Folded Light Curve (Smoothed)\nP = %.5f d (Selected)', refined_period));
    grid on;
end

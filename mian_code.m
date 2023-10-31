% Author: Syed Talha Abid Ali
% Date: Oct 31, 2023
% Paper Implementation: R-peak detection algorithm for ECG using double difference and RR interval processing 

%% To Select Signal from windows directory MIT-BIH Arrythmia Dataset of ECG from physionet (Download all signals in .mat format)
clc
clear all
close all

[filename,pathname] = uigetfile('*.*','Select the ECG Signal'); 
filewithpath = strcat(pathname,filename);
%Fs = input('Enter Sampling Rate: ');
Fs = 360;
ecg = load(filename);
ecgsig = (ecg.val)./200; 
t = 1:length(ecgsig);
N = length(ecgsig);
tx = t./Fs;
plot(ecgsig)

%% Removing DC Component

ecgsig_01 = ecgsig-mean(ecgsig);

%% Band Pass Filtering using Wavelet Decomposition (Removing Levels) [My own way to implement band pass filtering]

lev=4;
wt = modwt(ecgsig,4,'sym4'); 
wtrec = zeros(size(wt)); 
wtrec(3:4,:) = wt(3:4,:); 
y = imodwt(wtrec,'sym4'); 

%% Step 1: Calculate squared double differences

n = length(y);
d1 = diff(y);
d2 = diff(d1);
d = d2.^2 

for i = 1:n-1
    d1(i) = y(i+1) - y(i);
end

for j = 1:n-2
    d2(j) = d1(j+1) - d1(j);
    d(j) = d2(j)^2;
end

%% Step 2: Sorting 'd(j)' in descending order .
[d_sorted, d_indices] = sort(d, 'descend');
%% Step 3: Selecting peaks above a threshold
threshold = 0.1 * max(d); % 3% of the maximum value (I got good results on 0.1)
selected_peaks = d_indices(d_sorted > threshold);
%% Step 4: Eliminate peaks within ±75 ms of each other (considering them as closely spaced)

max_duration = 150; % Maximum duration of QRS regions in ms
time_interval = 75; % ±75 ms

final_peaks = selected_peaks(1); % Initialize with the first peak

for i = 2:length(selected_peaks)
    if (selected_peaks(i) - selected_peaks(i-1)) > time_interval
        final_peaks = [final_peaks, selected_peaks(i)];
    end
end
%% Step 5: QRS regions are within ±75 ms of each final peak (Creates Red Bounding boxes around potential QRS)

qrs_regions = [];
for i = 1:length(final_peaks)
    qrs_start = max(1, final_peaks(i) - time_interval);
    qrs_end = min(n, final_peaks(i) + time_interval);
    qrs_regions = [qrs_regions; qrs_start, qrs_end];
end

% Display the QRS regions or perform further analysis.
disp('QRS Regions:');
disp(qrs_regions);

% Plot the original ECG signal
figure;
plot(y);
title('Original ECG Signal');
xlabel('Sample Number');
ylabel('Amplitude');
hold on;

% Plot QRS regions as red rectangles
for i = 1:size(qrs_regions, 1)
    qrs_start = qrs_regions(i, 1);
    qrs_end = qrs_regions(i, 2);
    qrs_duration = qrs_end - qrs_start + 1;
    
    % Highlight QRS regions by drawing rectangles
    rectangle('Position', [qrs_start, min(y), qrs_duration, max(y) - min(y)], 'EdgeColor', 'r', 'LineWidth', 1);
end

legend('ECG Signal', 'QRS Regions');
hold off;


%% Step 6: Initializing arrays to store R-peak locations
r_peak_locations = [];
%% Step 7: Defining a window size for relative magnitude comparison
window_size = 70; % Adjust as needed
%% Step 8: Loop through each QRS region to find R Peaks
for i = 1:size(qrs_regions, 1)
    qrs_start = qrs_regions(i, 1);
    qrs_end = qrs_regions(i, 2);
    
    % Ensure QRS window does not extend beyond the signal
    qrs_start = max(qrs_start, 1);
    qrs_end = min(qrs_end, length(y));
    
    % Extract the QRS window
    qrs_window = y(qrs_start:qrs_end);
    
    % Calculate the maximum and minimum amplitudes
    max_amplitude = max(qrs_window);
    min_amplitude = min(qrs_window);
    
    % Calculate the mean of maximum and minimum amplitudes
    mean_amplitude = (max_amplitude + min_amplitude) / 2;
    
    % Calculate relative magnitudes by subtracting the mean from all data points
    relative_magnitudes = qrs_window - mean_amplitude;
    
    % Find the position of the maximum relative magnitude within the window
    [~, max_rel_magnitude_index] = max(relative_magnitudes);
    
    % Convert the relative index to the absolute index
    r_peak_location = qrs_start - 1 + max_rel_magnitude_index;
    
    % Add the R-peak location to the list
    r_peak_locations = [r_peak_locations, r_peak_location];
end
%% Step 9: Plot the original ECG signal with detected R peaks (Make Red Circle around the peaks)

figure;
plot(y);
title('Original ECG Signal with Detected R Peaks');
xlabel('Sample Number');
ylabel('Amplitude');
hold on;
plot(r_peak_locations, y(r_peak_locations), 'ro', 'MarkerSize', 10);
legend('ECG Signal', 'Detected R Peaks');
hold off;
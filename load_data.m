clearvars
clc
close all

addpath('input')

Test1 = struct2cell(load('Test1.mat'));
Test2 = struct2cell(load('Test2.mat'));

EMG_Test = {Test1{1}.matrix(:,1) Test1{1}.matrix(:,2) Test2{1}.matrix(:,1)};
EMG_Acceleration = {Test1{1}.matrix(:,3:5) Test1{1}.matrix(:,5:end) Test2{1}.matrix(:,2:end)};

% Parameters
Fs    = 2000; % Sampling Frequency 
Ts    = 1/Fs; % Sampling Interval 
dFs   = 500; % Down sampled Frequency
Fnyq  = Fs/2; % Nyquist Frequency
Freq1 = 30; % Passband Frequency, Hz (Lower)
Freq2 = 450; % Passband Frequency, Hz (Upper)
Fenv  = 4; % Envelope Frequency

% Filters Note:
% When specifying frequencies for digital filters in Matlab, the frequencies
% should be normalized by the Nyquist frequency. 
% The cutoff frequency (for envelope ) has to adjusted upward by 25% in 2nd
% order filter because the filter will be applied twice (forward and backward). 
% The adjustment assures that the actual -3dB frequency after two passes 
% will be the desired fco specified above. This 25% adjustment factor is 
% correct for a 2nd order Butterworth; 
% For a 4th order Butterworth used twice so multiply by 1.116.

W = (1 / Fnyq) * [Freq1, Freq2];
% [B, A]  = butter(4, 1.116 * W, 'Bandpass'); % Bandpass filter
[B, A] = fir1(300, W);

figure('units','normalized','outerposition',[0 0 1 1])
title('Bode Plot of the Filter');
freqz(B, A);
% freqz(B, A); % Filter Bode Plot

for i = 1:size(EMG_Test, 2)
    EMG_data = EMG_Test{i} - mean(EMG_Test{i}); % EMG Data, removing offset
    EMG_Filtered{i} = filtfilt(B, A, EMG_data); % Filter EMG
    rectified{i} = abs(EMG_Filtered{i}); % Rectify signal
    [b, a] = butter(4, 1.116 * (Fenv/ Fnyq), 'low');   
    EMG{i} = filtfilt(b, a, rectified{i}); % Smooth by LPF: 4th order
    
    %EMG_Acceleration{i} = EMG_Acceleration{i} - mean(EMG_Acceleration{i});
    EMG_ResultantAcc{i} = sqrt(sum(EMG_Acceleration{i}.^2')');
    
    figure('units','normalized','outerposition',[0 0 1 1])
    title('EMG Signal, Filtered and Rectified EMG Signal, Linear Envelope');
    grid;
   
    t = (1:max(size(EMG{i})))/Fs;
    lim = (max(size(t)) - 1000)/Fs;
    
    % Top two plots and bottom one
    subplot(2,2,1); 
    plot(t, EMG_data, 'b')
    hold on
    plot(t, EMG_Filtered{i} , 'r')
    hold off
    xlabel('Time (s)'); 
    xlim([0 lim]);
    ylabel('EMG (\muV)');
    legend('EMG signal', 'EMG filtered');
        
    subplot(2,2,2); 
    plot(t, rectified{i}, 'g')
    hold on
    plot(t, EMG{i}, 'r');
    hold off
    xlabel('Time (s)'); 
    xlim([0 lim]);
    ylabel('EMG (\muV)');
    legend('EMG rectified', 'Linear envelope');   
    
    dsFactor = round(Fs / dFs); % Down sampling factor
    EMG_downsampled{i} = downsample(EMG{i},dsFactor); 
    
    subplot(2,2,[3 4]); 
    yyaxis left
    xlabel('Time (s)'); 
    plot(downsample(t, dsFactor), EMG_downsampled{i}, 'LineWidth', 1.5);
    ylim([-50 130])
    ylabel('EMG (\muV)');
    hold on
    yyaxis right
    ylabel('Acceleration (m/s^2)')
    xlabel('Time(s)'); 
    xlim([0 lim]);
    plot(t, EMG_ResultantAcc{i}, 'b');
    ylim([0 4])
    hold off
    legend('Linear envelope', 'Acceleration');    
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% - - - - - - - -  - - - - - - QUESTION A - - - - - - - - - - - - - - - - %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% The downsamplig phase needs to be after the filtering all the noise will
% alias back in, raising the floor noise of the system. Even if there is no
% noise, downsample before filtering will lead to aliasing from frequency
% above the Nyquist rate.

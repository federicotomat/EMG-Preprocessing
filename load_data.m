clearvars
clc
close all

addpath('input')
if(~exist('results', 'dir' )), mkdir('results'), end

Test1 = struct2cell(load('Test1.mat'));
Test2 = struct2cell(load('Test2.mat'));

EMG.test   = {Test1{1}.matrix(:,1) Test1{1}.matrix(:,2) Test2{1}.matrix(:,1)};
EMG.labels = {'Biceps' 'Triceps' 'Deltoid'};
EMG.acc    = {Test1{1}.matrix(:,3:5) Test1{1}.matrix(:,5:end) Test2{1}.matrix(:,2:end)};

% Parameters
Fs    = 2000; % Sampling Frequency 
Ts    = 1/Fs; % Sampling Interval 
dFs   = 200;  % Down sampled Frequency
Fnyq  = Fs/2; % Nyquist Frequency
Freq1 = 30;   % Passband Frequency, Hz (Lower)
Freq2 = 450;  % Passband Frequency, Hz (Upper)
Fenv  = 4;    % Envelope Frequency
progr = 0;    % Variable used for the display percentage

dsFactor = round(Fs / dFs); % Down sampling factor

% Filters Note:
% When specifying frequencies for digital filters in Matlab, the
% frequencies should be normalized by the Nyquist frequency. 
% The cutoff frequency (for envelope ) has to adjusted upward by 25% in 2nd
% order filter because the filter will be applied twice (forward and 
% backward). 
% The adjustment assures that the actual -3dB frequency after two passes 
% will be the desired fco specified above. This 25% adjustment factor is 
% correct for a 2nd order Butterworth; 
% For a 4th order Butterworth used twice so multiply by 1.116.

% Use the function freqz to display the Bode plot of the filter.

W      = (1 / Fnyq) * [Freq1, Freq2]; % Bandpass frequency
[B, A] = fir1(300, W);                % Bandpass filter

for i = 1:size(EMG.test, 2)
    EMG.data         = EMG.test{i} - mean(EMG.test{i});        % EMG Data, removing offset
    EMG.filt{i}      = filtfilt(B, A, EMG.data);               % Filter EMG
    EMG.rectified{i} = abs(EMG.filt{i});                       % Rectify signal
    [b, a]           = butter(4, 1.116 * (Fenv/ Fnyq), 'low'); % Envelope filter   
    EMG.smooth{i}    = filtfilt(b, a, EMG.rectified{i});       % Smooth by LPF: 4th order
    %EMG.acc{i}       = EMG.acc{i} - mean(EMG.acc{i});
    EMG.acc{i}       = sgolayfilt((100*EMG.acc{i}), 3, 657);   % Filter the single component (and convert it in cm/s^2)
    EMG.resAcc{i}    = sqrt(sum(EMG.acc{i}.^2')');

 
    f = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    grid;
   
    % Parametes for plotting time
    t   = (1:max(size(EMG.smooth{i})))/Fs;
    lim = (max(size(t)) - 1000)/Fs;
    
    % Top two plots and bottom one for the EMG preprocessing
    subplot(2,2,1); 
    plot(t, EMG.data, 'b')
    hold on
    plot(t, EMG.filt{i} , 'r')
    hold off
    xlabel('Time (s)'); 
    xlim([0 lim]);
    ylabel('EMG (\muV)');
    legend('EMG signal', 'EMG filtered');
    title([EMG.labels{i} ': EMG recorded and filtered']);
        
    subplot(2,2,2); 
    plot(t, EMG.rectified{i}, 'g')
    hold on
    plot(t, EMG.smooth{i}, 'r');
    hold off
    xlabel('Time (s)'); 
    xlim([0 lim]);
    ylabel('EMG (\muV)');
    legend('EMG rectified', 'Linear envelope');   
    title([EMG.labels{i}, ': EMG rectified and envelope']);
    
    % Downsample the signal
    EMG.downsampled{i} = downsample(EMG.smooth{i}, dsFactor); 
    
    subplot(2,2,[3 4]); 
    yyaxis left
    xlabel('Time (s)'); 
    plot(downsample(t, dsFactor), EMG.downsampled{i}, 'LineWidth', 1.05);
    ylim([-50 130])
    ylabel('EMG (\muV)');
    hold on
    yyaxis right
    ylabel('Acceleration (cm/s^2)')
    xlabel('Time(s)'); 
    xlim([0 lim]);
    plot(t, EMG.resAcc{i}, 'LineWidth', 1.05);
    hold off
    legend('Linear envelope', 'Acceleration');   
    title([EMG.labels{i}, ': resultant acceleration and downsampled envelope']);
    
    saveas(f, ['results/' EMG.labels{i}, '_preprocessing.jpg']);
   
    % Acceleration component plot
    f = figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    grid;
    title('Components Acceleration');
    
    subplot(3,1,1); 
    plot(t, EMG.acc{i}(:,1));
    xlim([0 lim]);
    xlabel('Time (s)');
    ylabel('x-Acceleration (m/s^2)')
    title([EMG.labels{i}, ': acceleration along x']);
    
    subplot(3,1,2); 
    plot(t, EMG.acc{i}(:,2));
    xlim([0 lim]);
    xlabel('Time (s)');
    ylabel('y-Acceleration (m/s^2)')
    title([EMG.labels{i}, ': acceleration along y']);
    
    subplot(3,1,3); 
    plot(t, EMG.acc{i}(:,3));
    xlim([0 lim]);
    xlabel('Time (s)');
    ylabel('z-Acceleration (m/s^2)')
    title([EMG.labels{i}, ': acceleration along z']);
    
    saveas(f, ['results/' EMG.labels{i}, '_acceleration.jpg']);
    
    if i == 1 
        fprintf(1,'   EMG Processing: %3d%%', progr);
    end
    progr = floor((100*(i/(size(EMG.test, 2)))));
    fprintf(1,'\b\b\b\b%3.0f%%', progr);
end

% Save processed data into a .mat file
EMGprocessed.labels = EMG.labels;
EMGprocessed.matrix = EMG.downsampled;
save('results/Output.mat', 'EMGprocessed');

fprintf('\n   The results of the pre-processing are available in the *results* folder\n');

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %
% - - - - - - - -  - - - - - - QUESTION A - - - - - - - - - - - - - - - - %
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - %

% The downsamplig phase needs to be after the filtering all the noise will
% alias back in, raising the floor noise of the system. Even if there is no
% noise, downsample before filtering will lead to aliasing from frequency
% above the Nyquist rate.

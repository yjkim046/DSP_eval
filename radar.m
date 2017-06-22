% Given Constraints
tgt_pd = 0.9;            % Probability of detection
tgt_pfa = 1e-6;          % Probability of false alarm
max_range = 5000;    % Maximum unambiguous range
range_res = 50;      % Required range resolution
rcs = 1;             % Required target radar cross section

% EM Parameters
c = 299792458; % light speed
fc = 100e6; lambda = c/fc;
Pt = 5000; % Tx power of 5kW after Tx Gain

% Pulse Parameters
repetition = 10;
pulse_width = (2*range_res)/c;
pulse_interval = 2*max_range/c;
pulse_bw = 1/pulse_width;
Fs = pulse_bw; Ts = 1/Fs; % Sampling Frequency

% Rx Parameters
detection_threshold = 0.2e-13; % Needs to be adjusted depending on SNR

num_pulses = 100000;
pwr_detected = zeros(1,num_pulses);
% PD
for ll = 1:num_pulses
    % Tx
    x = 1; % 1 for Probability of Detection, 0 for PFA

    % Propagation Channel and noise addition
    dist = max_range;
    Aeff = lambda^2/4/pi;
    Pr = Pt/(4*pi*dist^2)*rcs/(4*pi*dist^2)*Aeff;
    Pr_dBm = 10*log10(Pr/1e-3);
    rxWaveform = sqrt(Pr)*x;

    % Noise addition
    noise_dBm = -174 + 10*log10(Fs);
    N = 10^(noise_dBm/10)*1e-3/repetition;
    v = sqrt(N/2)*complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
    rxWaveform = rxWaveform + v;

    % Rx
    pwr_detected(ll) = abs(rxWaveform).^2;
end
Pwr1 = mean(pwr_detected);
PD = sum(pwr_detected>detection_threshold)/num_pulses;

% PFA
for ll = 1:num_pulses
    % Tx
    x = 0; % 1 for Probability of Detection, 0 for PFA

    % Propagation Channel and noise addition
    dist = max_range;
    Aeff = lambda^2/4/pi;
    Pr = Pt/(4*pi*dist^2)*rcs/(4*pi*dist^2)*Aeff;
    Pr_dBm = 10*log10(Pr/1e-3);
    rxWaveform = sqrt(Pr)*x;

    % Noise addition
    noise_dBm = -174 + 10*log10(Fs);
    N = 10^(noise_dBm/10)*1e-3/repetition;
    v = sqrt(N/2)*complex(randn(size(rxWaveform)), randn(size(rxWaveform)));
    rxWaveform = rxWaveform + v;

    % Rx
    pwr_detected(ll) = abs(rxWaveform).^2;
end
Pwr0 = mean(pwr_detected);
PFA = sum(pwr_detected>detection_threshold)/num_pulses;

% Print the results
SNR_dB = Pr_dBm - noise_dBm;
fprintf('At SNR = %4.2f (dB)\n',SNR_dB);

fprintf('Probability of Detection = %4.2f',PD);
if PD>tgt_pd
    fprintf('.....PASSED\n');
else
    fprintf('.....FAILED\n');
end

fprintf('Probability of False Alarm = %4.2f',PFA);
if PFA<tgt_pfa
    fprintf('.....PASSED\n');
else
    fprintf('.....FAILED\n');
end

%% precursor

%files and info identities
StudRegNum = 202113405;
s1 = ('speaker1.wav');
s2 = ('speaker2.wav');

%testing audio file
u1 = audioread(s1); % read in a single channel signal
u2 = filter([0 0 1],1,u1); % delay by two sampling periods
fs = 8000; % sampling frequency in Hertz
sound([u1 u2],fs); % play back stereo
%sound([u2 u1],fs); % reverse stereo channels

[x1,x2] = AssignmentScenario(StudRegNum);
sound(x1,8000);

%% Question 1

%from plot this is a discrete time signal. therefore period equals time between 2 points
x1F200Samples = x1(1:200);
t = (1:200)/fs;
fo = (2*pi)/Ts; % signal angular freq = 50,265.482
foeqn = (2*pi*fo/fs);
%time is 1/freq = 2.25ms from plot
sampleTime = 0.00225;


figure (1);
plot(t,x1F200Samples,'r--');
title("Question 1");
xlabel("Sample time (200 Samples)");
ylabel("x1 from audio");

%angular frequency using 2*pi*f  = 5.026548245743669e+04
%angular frequency using 2*pi*8000/angularfreq(above) = 39.4784 Rad/s
%fo from 1/Ts in plot = 1/Tgraph = 1/0.00225 = 444.44hz - fo from time
%domain

%plot generated

%% Question 2

fs = 8000; % sampling frequency in Hertz

%FIRST 200 Samples
%f = (0:length(x1F200Samples)-1)/length(x1F200Samples)*fs; % frequency scale from first 200 sampies
%plot(f,abs(fft(x1F200Samples))); % discrete Fourier transform of x
%estimated frequency = 7560 - 440 =7120 (discrete fourier signal so f = distance between 2 peaks)
%COMMENT - massive difference of 7080.5216 (7120 - 39.4784)

%ALL SAMPLES
figure(2);
f = (0:length(x1)-1)/length(x1)*fs; % frequency scale
plot(f,abs(fft(x1))); % discrete Fourier transform of x
xlabel('frequency f / [Hz]'); % label axes
ylabel('magnitude of discrete fourier transform');
%plot generated, peaks at 440Hz and 7120Hz
%estimated frequency of peaks in plot generated
%estimated frequency = 7555.56 - 444.44 = 7111.12 = 124.112rad/s (discrete fourier signal so f = distance between 2 peaks)

%COMMENT - massive difference of 84.6336Rad/s (124.112 - 39.4784) or
%4849.15Hz when comparing foeqn with fo found in plot 2
%unsure as to why this is the case tbh something to do with fourier
%alignment

%% Question 3

%Over the range of sample indices 19500 < n ≤ 21548,
%speaker 2 is relatively quiet. From the time domain signals x1[n] and x2[n], try to estimate the relative delay between
%these two signals.

%plot for x1
%set x1 and x2 to all sample constraints in Q3
x1fsamplesQ3 = x1(19500:21548);
x2fsamplesQ3 = x2(19500:21548);

%setting time to constraint parameters to accomodate for samples in in Q3
%so with parameters you are dealing with x1 and x2 as a whole so set yur
%time parameter to 0.000225 * 35999(number of x1 and samples)
%time delay for full of x1 and x2  = 4.499875
%tQ3 = (0:Ts:0.256);
tQ3 = (19500:21548)/fs;
%time delay for full of x1 and x2  = );

figure(3);
plot(tQ3,x1fsamplesQ3);
title("analysis q3")
xlabel("time(s)"); % label axes
ylabel("sample parameters");
hold on
plot(tQ3,x2fsamplesQ3);
title("analysis q3")
xlabel("time(s)"); % label axes
ylabel("sample parameters");
legend('x1fsample', 'x2fsample');
%xlim([19501 19508]) - delay = 0.26ms

%%attempted to estimate 2 similar peaks then figure out delay between them.


% Calculate cross-correlation
cross_corr = xcorr(x2fsamplesQ3, x1fsamplesQ3);
% Find the index of the maximum correlation
[max_corr, max_index] = max(cross_corr);
% Estimate the relative delay in terms of sample indices
relative_delay = max_index - length(x1fsamplesQ3);
fprintf('Relative delay between x1 and x2: %d samples\n', relative_delay);
%comes out to 1 sample


%% Q4
%For the segments in Question 3, let Xm(ejΩ),
%m = 1, 2, be their discrete time Fourier transforms. Assuming that X1(ejΩ) 6= 0∀Ω, what would you expect for the
%quantities

% Calculate Fourier transforms
X1 = fft(x1fsamplesQ3);
X2 = fft(x2fsamplesQ3);

% Calculate Gain (G) and Phase (A)
G = abs(X2 ./ X1);
A = angle(X2 ./ X1);
fx1 = (0:(length(X1)-1))/length(X1)*fs;
fx2 = (0:(length(X2)-1))/length(X2)*fs;

% Plot the results
figure;
subplot(2, 1, 1);
plot(fx1, G);
title('Gain vs Frequency');
xlabel('Frequency (\Omega)');
ylabel('Gain');

subplot(2, 1, 2);
plot(fx1, A);
title('Phase vs Frequency');
xlabel('Frequency (\Omega)');
ylabel('Phase');

%Display the results for a specific frequency (e.g., at the peak)
peak_frequency_index = find(G == max(G)); %finding sample number where frequency is highest
fprintf('Gain at peak frequency: %.4f\n', G(peak_frequency_index)); %gain for sample number 
fprintf('Phase at peak frequency: %.4f radians\n', A(peak_frequency_index)); %angle for sample number
%Also could be guessed from graph estimation but this is still the same
%principle justified - periodic angle and discrete magnitude

%% Q5

fs = 8000; % sampling frequency in Hertz
Range = (19501:21548); % window of samples to investigate
N = length(Range); % number of samples
X1 = fft(x1(Range));% Fourier coefficients via fft()
X2 = fft(x2(Range));
f = (0:(N-1))/N*fs; % frequency scale
G_Q5 = abs(X2) ./ abs(X1); % Fourier coefficients via fft()

%Mag plot
figure();
plot(f,abs(X1)); % plot magnitude of coefficient
xlabel('frequency f / [Hz]'); % label axes
ylabel('magnitude');
title('Mag diff');
hold on;
plot(f,abs(X2)); % plot magnitude of coefficient
legend('X1', 'X2');

figure();
plot(f,abs(X1)); % plot magnitude of coefficient
xlabel('frequency f / [Hz]'); % label axes
ylabel('magnitude');
hold on;
plot(f,abs(X2)); % plot magnitude of coefficient
title('Mag diff');
legend('X1', 'X2');
xlim([199.218 199.220]); %mag diff taken from plot is 0.9
%xlim([415.300 415.320]);

%Gain Plot
figure();
plot(f,G_Q5); % plot magnitude of coefficient
xlabel('frequency f / [Hz]'); % label axes
ylabel('Gain');

%The gain and delay difference plots will provide insights into how the signals
%x1[n] and x2[n]​ differ in terms of magnitude and phase across different frequencies.
%from mag plot gain diff = 0.5528, taken the peaks and checked y axis
%difference
%from angle plot delay estimate = 0.3 radians, 17.1887 degrees


%Delay
A3Q5 = angle(X2./X1);
figure();
plot(f,A3Q5); % plot magnitude of coefficient
xlabel('frequency f / [Hz]'); % label axes
ylabel('Angle of delay Rad/s');
hold on;
mq5 = -0.0016207;
yq5= mq5 *f;
plot(f,yq5); % plot magnitude of coefficient
xlabel('frequency f / [Hz]'); % label axes
ylabel('Angle of delay Rad/s');



%part 2
AQ5 = angle(X1); %estimated angle found - 30 deg
figure();
plot(f, AQ5);
title('Phase Angle vs Frequency');
xlabel('Frequency (Hz)');
ylabel('Phase Angle (radians)');

%spikes between 3 and -3 so find 3 radians in degrees = +/- 171.887 degrees
%could maybe take one angle away from other to understand the phase delay,
%maybe inverse the fourier transform to the time domain and then estimate
%delay from a plot, gain is a little different as you can evauluate from
%plots when dealing with magnitude, find gain of q5 and q6 then evaluate
%the peaks?


%% Q6  (Estimation for 2nd speaker.)
% The first speaker is quiet during the period 24800 < n ≤ 25824.
% Use this to estimate the gain and delay difference between x1[n] and x2[n] for speaker 2.
% W.r.t. Fig. 1, what is the estimated angle ϑ2?

Range_Q6 = (24800:25824);      % window of samples to investigate - 1025
N_Q6 = length(Range_Q6);       % number of samples
X1_Q6 = fft(x1(Range_Q6));     % Fourier coefficients via fft()
X2_Q6 = fft(x2(Range_Q6));
f_Q6 = (0:(N_Q6-1))/N_Q6*fs;   % Frequency scale
t_Q6 = (0:(N_Q6-1))/N_Q6*(1025*sampleTime);    % Time scale for inverse 
theta_Q6=angle(X2_Q6);

delta_Q6= angle(X2_Q6)-angle(X1_Q6);
G_Q6=abs(X2_Q6)./abs(X1_Q6);  % Gain between between 2 signals in range

figure();
subplot(2,1,1);
plot(f_Q6,abs(X2_Q6)); % plot magnitude of coefficient
title("fourier magnitude of X2");
xlabel("frequency f / [Hz]"); % label axes
ylabel("magnitude");

subplot(2,1,2);
plot(f_Q6,abs(X1_Q6)); % plot magnitude of coefficient
title("fourier magnitude of X1");
xlabel("frequency f / [Hz]"); % label axes
ylabel("magnitude");

%mag diff from plots = 0.1569 - is this gain difference? may need to find
%the difference in magnitude between both transforms although mag is not
%equal to gain


%Gain Diff - dont think this is necessary

figure();
plot(f_Q6,G_Q6); % plot magnitude of coefficient
title("Gain of both sample ranges");
xlabel("frequency f / [Hz]"); % label axes
ylabel("Gain");
hold on;
plot(f,G_Q5);
figure();
plot(f_Q6,G_Q6); % plot magnitude of coefficient
title("Gain but zoomed in");
xlabel("frequency f / [Hz]"); % label axes
ylabel("Gain");
hold on;
plot(f,G_Q5);
legend('GQ6', 'GQ5');
xlim([1800 1815]); % at similar peak GQ6 - 1810.73 & GQ5 - 1804.69 gain diff = 6.04

%plotting angles
figure();
subplot(2,1,1);
plot(f_Q6,angle(X1_Q6)); % plot magnitude of coefficient
title("fourier angle of X1");
xlabel("frequency f / [Hz]"); % label axes
ylabel("angle");

subplot(2,1,2);
plot(f_Q6,theta_Q6); % plot magnitude of coefficient
title("fourier angle of X2");
xlabel("frequency f / [Hz]"); % label axes
ylabel("angle");

%delay difference angle
phaseQ6 = angle(X2_Q6) - angle(X1_Q6);
AQ6 = unwrap(phaseQ6);
figure();

plot(f_Q6, AQ6);
title('Delay Difference vs Frequency');
xlabel('Frequency (Hz)');
ylabel('Delay Difference (radians)');
xlim([0 8000]); %delay = gradient change in 1/8th of frequency so 1000hz - 2.34798RAD/s or 0.3737ms



%% Q7 (Estimation of relative transfer functions.) 
% Using the estimated relative gains and delays for the two speakers
% determine the matrix of relative transfer functions H (z) •—◦ H[n]

%Fourier transform function= G(X2)/G(X1)
GainH_Q7=abs(G_Q6/G_Q5);
PhaseH_Q7=angle(theta_Q6/thetaQ5);
Hz=[GainH_Q7, PhaseH_Q7];
%[GainH_Q7 PhaseH_Q7]

%Transfer Function
% Given gain and delay differences for both speakers use estimations
tau_Q5 = 171.887;   % Replace with your actual delay difference for speaker 1
tau_Q6 = 171.887; % Replace with your actual delay difference for speaker 2

% Angular frequency vector
omega = 2 * pi * f; % f is the frequency vector from your analysis

% Calculate the relative transfer functions
H_Q5 = G_Q5 * exp(-1i * omega * tau_Q5);
H_Q6 = G_Q6 * exp(-1i * omega * tau_Q6);

% Matrix of relative transfer functions
H = [H_Q5; H_Q6];




%% Question 8 (Construction of separation filters.)
% Based on your estimated H (z), construct and implement the unmixing matrix G(z) = inv(H(z))
% What do the output signals of G(z) sound like?

%Gz=inv(Hz); no idea,, idk if need to implement below code from
%Appendix A or nah. 


%Look at Q10 in advance its asking to adjust the rho value, therefore maybe
%need to adjust it.

%Matrix inversion
%Hinv = pinv(H);
%sound(abs(Hinv));

[Hinv] = AssignmentScenario(StudRegNum);
sound(x1,fs);

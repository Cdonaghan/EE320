%% Semester 1 EE320 Assignment
%%Precursor
Studentregno=202113405;
u1=audioread("speaker1.wav");
u2=audioread("speaker2.wav");
fs = 8000; % sampling frequency 

%sound(u2,fs);
%sound(u1,fs);
%sound([s1 s2],fs); 
u2=filter([0 0 1],1,u1);
%sound([u1 u2],fs);
[x1,x2]=AssignmentScenario(Studentregno);
%sound(x1,fs);
%sound(x2,fs);

%% Q1 (Time domain analysis of interference frequency.) 

% range first 200 elements
x1_first_200 = x1(1:200);

tQ1=(1:200)/fs;  

figure(1);
plot(tQ1,x1_first_200);
xlim ([0 0.01])
title("Q1. 200 samples x[1]")
xlabel("Time(s)");
ylabel("Amplitude");

%fo from 1/Ts in plot = 1/Tgraph = 1/0.00225 = 444.44hz - fo from time
%domain
%in angular domain 444.44hz = 15.51 Rad/s
%plot generated



%% Q2 (Frequency domain analysis of interference frequency.) 
%  Estimate the interference frequency by evaluating the discrete-time Fourier transform

f = (0:length(x1)-1)/length(x1)*fs;
plot(f,abs(fft(x1)));                      
title("Q2. Frequency Domain analysis")
xlabel("Frequency (Hz)");          
ylabel("Magnitude of Fourier Transform");


%% Q3 (Frequency domain analysis of interference frequency.) 
%  Estimate the interference frequency by evaluating the discrete-time Fourier transform
%Over the range of sample indices 19500 < n ≤ 21548,
%speaker 2 is relatively quiet. From the time domain signals x1[n] and x2[n], try to estimate the relative delay between
%these two signals.

f = (0:length(x1)-1)/length(x1)*fs;
plot(f,abs(fft(x1)));                      
title("Q2. Frequency Domain analysis")
xlabel("Frequency (Hz)");          
ylabel("Magnitude of Fourier Transform");

%plot for x1
%set x1 and x2 to all sample constraints in Q3
x1fsamplesQ3 = x1(19500:21548);
x2fsamplesQ3 = x2(19500:21548);

%setting time to constraint parameters to accomodate for samples in in Q3
%so with parameters you are dealing with x1 and x2 as a whole so set yur
%tQ3 = (0:Ts:0.256);
tQ3 = (19500:21548)/fs; % - this works better

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
%xlim([19501 19508]) - delay = 0.26ms - 1 sample
%attempted to estimate 2 similar peaks then figure out delay between them.

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
%fx2 = (0:(length(X2)-1))/length(X2)*fs;

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

%%taken from MATLAB as another way to justify Gain and Angle
%Display the results for a specific frequency (e.g., at the peak)
peak_frequency_number = find(G == max(G)); %finding sample number where frequency is highest
fprintf('Gain at peak frequency: %.4f\n', G(peak_frequency_number)); %gain for sample number 
fprintf('Phase at peak frequency: %.4f radians\n', A(peak_frequency_number)); %angle for sample number
%Also could be taken from graph estimation but this is still the same
%Handwritten analysis for justification - shown in document file

%% Q5 (Delay and gain estimation in the Fourier domain.)

Range = (19501:21548);
N = length(Range);     
X1 = fft(x1(Range));   
X2 = fft(x2(Range));   
f_Q5 = (0:(N-1))/N*fs; 

%Gain and angle
G_Q5=abs((X2)./(X1));
phase1= angle(X2./X1);

figure();
plot(f_Q5,abs(X1)); hold on;
plot(f_Q5,abs(X2)); hold off;
avgGain=mean(G_Q5);

legend ('X1','X2');
%xlim([1144.4 1144.7]);
title("Q5. Fourier magnitude plots");
xlabel("Frequency (Hz)");
ylabel("Magnitude");


figure();
plot(f_Q5,G_Q5); hold on;
yline(1/avgGain,"r--","LineWidth",2); hold off;
ylim([-5 35])
title("Q5. Gain Plot");
xlabel("Frequency (Hz)"); 
ylabel("Gain");

figure();
m1=-0.0016207;
y=m1*f_Q5;
plot(f_Q5,phase1); hold on;
plot(f_Q5,y,"r--");
title("Q5. Getting angle");
xlabel("Frequency (Hz)"); 
ylabel("Delay Difference (radians)");

%gain = 1/1.336 = 0.8821
%delay = 2 sample - 0.0005
%angle = 31.027

%% Q6  (Estimation for 2nd speaker.)

Range_Q6 = (24800:25824);      
N_Q6 = length(Range_Q6);       
X1_Q6 = fft(x1(Range_Q6));     
X2_Q6 = fft(x2(Range_Q6));
f_Q6 = (0:(N_Q6-1))/N_Q6*fs;  

% Gain and angle Calculations
G_Q6=abs(X2_Q6)./abs(X1_Q6);
phase2= angle(X2_Q6./X1_Q6);
avgGainQ6=mean(G_Q6);
G_Q6=abs(X2_Q6)./abs(X1_Q6);  

%plots
figure();
plot(f_Q6,abs(X1_Q6)); hold on;
plot(f_Q6,abs(X2_Q6)); hold off;
legend ('X1','X2');
title("Q6. Fourier magnitude plots");
xlabel("Frequency (Hz)");
ylabel("Magnitude");


figure();
plot(f_Q6,G_Q6); hold on;
yline(1/avgGainQ6,"r--","LineWidth",2); hold off;
ylim([-5 35])
title("Q6. Gain Plot");
xlabel("Frequency (Hz)"); 
ylabel("Gain");


figure();
m2=2.33405e-3;
y2=m2*f_Q6;
plot(f_Q6,phase2); hold on;
plot(f_Q6,y2); hold off;
title("Q6. Getting angle");
xlabel("Frequency (Hz)"); 
ylabel("Delay Difference (radians)");


% Gain: Avg_Gain=1.1731.  1/1.1731= 0.85244
% Delay Difference: 3.46e-4 : 3 samples
% Angle= 47.87

%% Q7 (Estimation of relative transfer functions.) 

% 'H' Transfer function delays from q5 and q6
tau_Q5 = 0.00025;    
tau_Q6 = 0.000375;   

% Angular frequency vector
omega = 2 * pi * f;  % f is the frequency vector from your analysis

% Calculate the relative transfer functions
H_Q5 = G_Q5 * exp(-1i * omega * tau_Q5);
H_Q6 = G_Q6 * exp(-1i * omega * tau_Q6);

% Matrix of relative transfer functions
H = [H_Q5; H_Q6];
%Handwritten analysis to justify this implementation, in document


%% Question 8 (Construction of separation filters.)

rho1= 1/1.133646820460025;
rho2= 1/1.173089387767823;
rho3= rho1 * rho2 ;

% delay1= 2 samples
% delay2= 3 samples

y1prime = filter(1,[1 0 0 0 -rho3],x2) + filter([0 -rho1],[1 0 0 0 -rho3],x1);
y2prime = filter(-rho2, [1 0 0 -rho3],x2) + filter(1,[1 0 0 0 -rho3],x1);

%sound(y1prime);
%sound(y2prime);
%it inverses the signal priority

%% Question 9 (Bonus – efficient implementation.) 
%%woodbury  formula

%testmatrix
omega1 = 2 * pi * fs;
H1Q9 = 0.8821 * exp(-1i * omega1 * tau_Q5);
H2Q9 = 0.85244 * exp(-1i * omega1 * tau_Q6);
H3Q9 = [H1Q9; H2Q9]; %gives [0.882100000000000 + 0.000000000000000i;
%0.852440000000000 + 0.000000000000001i]

%inversion
W9 = 0.8821;
X9 = exp(-1i * omega1 * tau_Q5);
Y9 = 0.85244;
Z9 = exp(-1i * omega1 * tau_Q6);
inverseE = inv(W9);
inverseG = inv(Y9);
inversenormal = inv(W9+X9*Y9*Z9);
inversewoodbury = inverseE - inverseE*X9*inv(inverseG + Z9*inverseE*X9)*Z9*...
inverseE;
%both inmversewoodbury 1 and 2 come out to the same value of 0.5765 + 0.0000i 
%thus proving that this theory can impliment efficient matrix inversion 

%% Question 10 (Notch filter anaysis.)
% Determine the denominator and numerator coefficients ai and bi, i =0, 1, 2
%should maybe be a function to extract it??

%plot 1
Omega0 = 444.444; rho4 = 0.9; %frequency from q1
a = [1 -2*rho4*cos(Omega0) rho4^2];
b = [1 -2*cos(Omega0) 1];
Omega = (0:8000)/8000*2*pi; % freq. scale
Q = (b(1) + b(2)*exp(-1i*Omega) + b(3)*exp(-1i*2*Omega)) ./...
(a(1) + a(2)*exp(-1i*Omega) + a(3)*exp(-1i*2*Omega));
figure();
subplot(2, 1, 1);
plot(Omega/2*pi*8000,abs(Q));
title("Notch Filter r = 0.9");
xlabel('freq. f/[Hz]'); ylabel('magnitude');

%plot 2
rho5 = 0.99;
r = [1 -2*rho5*cos(Omega0) rho5^2];
l = [1 -2*cos(Omega0) 1];
P = (l(1) + l(2)*exp(-1i*Omega) + l(3)*exp(-1i*2*Omega)) ./...
(r(1) + r(2)*exp(-1i*Omega) + r(3)*exp(-1i*2*Omega));
subplot(2, 1, 2);
plot(Omega/2*pi*8000,abs(P));
title("Notch Filter r = 0.99");
xlabel('freq. f/[Hz]'); ylabel('magnitude');

%% Question 11 (Notch filter implementation and application.)
% Implement the notch filter in Fig 2. How do the output signals y1[n] and y2[n] 
% sound when using the notch filters with different γ = 0.9 and γ = 0.99?
[x1,x2] = AssignmentScenario(202113405);

% Implementing Notch Filter
 r = 0.99; r2 = 0.9; omega0 = 2*pi/5; 
% normalised frequency of sinusoidal interference

% coefficents of denominator 
a = [1 -2*r* cos(omega0) r^2]; 
a2 = [1 -2*r2* cos(omega0) 1, r2^2]; 

% coefficents, of numerator
b = [1 -2* cos(omega0) 1]; 

omega = (0:1000)/1000* 2* pi;

Q = zeros(length(omega)); 
Q2 = zeros(length(omega)); 

for k = 1:length(omega)

Q(k) = (b(1) + b(2)* exp(-1i*omega(k)) + b(3)* exp(-1i* 2*omega(k))) /...
(a(1) + a(2)* exp(-1i*omega(k)) + a(3)*exp(-1*2* omega(k)));

Q2(k) = (b(1) + b(2)* exp(-11*omega (k)) + b(3)* exp(-1i*2*omega(k))) / ...
(a2(1) + a2(2)*exp(-1i*omega(k)) + a2(3)*exp(-1i*2*omega(k))); 

end
figure();
subplot(2,1,1) 
plot(omega/(2* pi)* 8000, abs(Q)); 
title("Notch Filter r = 0.99"); 
xlabel("frequency [Hz]"); 
ylabel("magnitude"); 
subplot(2,1,2)
plot(omega/(2* pi)* 8000,abs(Q2)); 
title("Notch Filter r = 0.9");
xlabel("Frequency (Hz]"); 
ylabel("magnitude");

% seperated signals from Q7
rho1 = 1/1.133646820460025;
rho2 = 1/1.173089387767823;

rho3 = rho1* rho2;% multiplication of coefficients used because it is an inverse matrix
y1_prime = filter (1, [1 0 0 0 -rho3],x2) + filter ([0 0 0 -rho1], [1 0 0 0 rho3],x1);
y2_prime = filter(-rho2,[1 0 0 -rho3],x2) + filter(1,[1 0 0 0 -rho3],x1);

% Filtering interference from signals

y1 = filter(b,a,y1_prime); 
y2 = filter (b,a,y2_prime);
%sound(y1);
%sound(y2);

%% Question 12 (Filter sequence.) 
%answer is explained in document



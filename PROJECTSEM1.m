%% precursor

%files and info identities
StudRegNum = 202113405;
s1 = ('speaker1.wav');
s2 = ('speaker2.wav');

%testing audio file
u1 = audioread(s1); % read in a single channel signal
u2 = filter([0 0 0 1],1,u1); % delay by two sampling periods
fs = 8000; % sampling frequency in Hertz
sound([u1 u2],fs); % play back stereo
%sound([u2 u1],fs); % reverse stereo channels

[x1,x2] = AssignmentScenario(StudRegNum);
sound(x1,8000);

%% Question 1
%time is 1/freq = 0.125ms

Ts = 1/fs; % sampling period [s]
t = (0:Ts:0.00025); % sampled time scale over 200 samples
fo = (2*pi*fs); % signal angular freq
stem(t,fo);


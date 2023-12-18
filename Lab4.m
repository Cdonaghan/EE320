%% time domain 1

t = (0:0.01:10);
h = StepFct(t,1) - StepFct(t,2);
x = ImpulseFct(t,0);
y = FilterFct(h,x,t);
plot(t,y);

%% time domain 2


t = (0:0.01:10);
h = StepFct(t,0) - StepFct(t,2);
x = ImpulseFct(t,0);
y1 = FilterFct(h,x,t);
y2 = FilterFct(h,y1,t);
plot(t,y2);
y3 = FilterFct(h,y2,t);

%% fourier domain 1
[H,f] = FourierTransform(h,t);

%returns the Fourier transform in H and a frequency scale (in Hertz) in f.
%You can display e.g. magitude
figure(1);
plot(f,abs(H));
figure(2);
%and phase
plot(f,angle(H));

%changing stepfunction 
figure(3);
[H1,f] = FourierTransform(h1,t);
h1 = StepFct(t,1) - StepFct(t,2);
plot(f,abs(H1));
figure(4)
plot(f,angle(H1));


%changing stepfunction again
figure(5);
h2 = StepFct(t,2) - StepFct(t,3);
[H2,f] = FourierTransform(h2,t);
plot(f,abs(H2));
figure(4)
plot(f,angle(H2));




%% Fourier Domain 2 
%functions are defined in the time domain:

t = (0:0.01:10);
f1 = StepFct(t,0) - StepFct(t,2);
f2 = t.*(StepFct(t,0) - StepFct(t,2));
%◮ display the two functions, and their Fourier transforms
[F1,f] = FourierTransform(f1,t);
[F2,f] = FourierTransform(f2,t);
%◮ compute the convolution, e.g. by
g = FilterFct(f1,f2,t);

[G,f] = FourierTransform(g,t);
fourierconvolution = F1.*F2;

figure(1)
subplot(211);
plot(f,abs(fourierconvolution));
hold on;
plot(f,abs(G), 'r--');
subplot(212);
plot(f,angle(fourierconvolution));
hold on;
plot(f,angle(G), 'r--');







function [X,f] = FourierTransform(x,t)
% X = FourierTransform(x,t);
%
% Fourier transform implementation --- this is a discrete approximation.
%  
% S. Weiss, 12/11/2013

Ts = t(2)-t(1);          % sampling period
X = fft(x)*Ts;           % Fourier transform ---scaled for continuous approx.
f = (0:(length(t)-1))/length(t)/Ts;     % freq. scale

L = ceil(length(x)/2);
X = [X(L+1:end) X(1:L)];
f = f - f(L);

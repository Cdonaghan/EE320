function y = FilterFct(h,x,t)
% s = ImpulseFct(t,tau);
%
% Approximation of a continuous convolution.
% 
% S. Weiss, 12/11/2013

y = filter(h,1,x) * (t(2)-t(1));

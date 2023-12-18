function s = ImpulseFct(t,tau)
% s = ImpulseFct(t,tau);
%
% Implementation of the discrete version of an impulse delta(t-tau).
% 
% S. Weiss, 12/11/2013

s = (t==tau);  % the logic comparison returns 1 if true,
              % and 0 otherwise.

s = s / (t(2)-t(1));  % height of impulse = 1/ width

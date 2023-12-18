function s = StepFct(t,tau)
% s = StepFct(t,tau);
%
% Implementation of the step function s(t-tau).
% 
% S. Weiss, 27/10/2013
 
s = (t>tau);  % the logic comparison returns 1 if true,
              % and 0 otherwise.

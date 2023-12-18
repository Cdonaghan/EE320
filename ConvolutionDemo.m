function ConvolutionDemo()
% ConvolutionDemo();
%
% Visualisation of the convolution of two simulated
% continuous-time functions f1 and f2.
% To generate different functions, please edit f1
% and f2 below as appropriate.
%
% S Weiss, 27/10/2013

% time scale
t = (-1:0.01:10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ARGUMENT 1 TO CONVOLUTION
%
%   insert your function f1 here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=2*(StepFct(t,0) - StepFct(t,.5));
f4 = StepFct(t,0) - StepFct(t,1);
f3 = t.*(StepFct(t,0) - StepFct(t,1));
f1=f3;
f2 = f4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ARGUMENT 2 TO CONVOLUTION
%
%   insert your function f2 here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = StepFct(t,0) -StepFct(t,1.5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
subplot(311);
plot(t,f1,'b-','linewidth',2);
hold on;
plot(t,f2,'r--','linewidth',2);
legend('f_1(t)','f_2(t)');
ylabel('amplitude'); 
Range = axis();
axis([-1 10 Range(3) Range(4)]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NUMERICAL SOLUTION TO CONVOLUTION
%   AND DISPLAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2=fliplr([zeros(1,900) f2]);

y = []; l=1; L=5;
for T = 0:0.01:10,
   product= f1.*f2(901:end);    
   f2=[0 f2(1:end-1)];
   yT = sum(product)*0.01;
   y = [y yT];

   % refresh plot every L-th step
   if mod(l-1,L)==0,
      subplot(312);
      hold off;
      plot(t,f1,'b-','linewidth',2);
      hold on;
      plot(t,f2(901:end),'r--','linewidth',2);
      plot(t,product,'k-.','linewidth',2);
      area(t,product,'FaceColor','y');
      plot(t,f1,'b-','linewidth',2);
      plot(t,f2(901:end),'r--','linewidth',2);
      plot(t,product,'k-.','linewidth',2);
      plot(T,0,'ko','linewidth',3);
      Range=axis();
      axis([-1 10 Range(3) Range(4)]);
      grid on;
      legend('f_1(t)',sprintf('f_2(%.1f-t)',T),sprintf('f_1(t) f_2(%.1f-t)',T));
      ylabel('amplitude');
      subplot(313);
      hold off;
      plot((0:0.01:T),y,'k','linewidth',2);
      hold on;
      plot(T,yT,'ko','linewidth',3);
      Range=axis();
      axis([-1 10 Range(3) Range(4)]);
      grid on;
      legend('y(t)');
      ylabel('amplitude');
      xlabel('time t/s');
      drawnow;
   end;
   l = l+1;
end;


function s = StepFct(t,tau)
% Implements the step function s(t-tau)
s = (t>tau);  

function d = ImpulseFct(t,tau);
% Implememts a Dirac Impulse
d = (t==tau)/(t(2)-t(1));

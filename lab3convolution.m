
t = [-1:0.01:10];

f1 = StepFct(t,1);




f2 = StepFct(t,0) - StepFct(t,1);




f3 = t.*(StepFct(t,0) - StepFct(t,1));



f4 = (t-2).*(StepFct(t,2) - StepFct(t,3));



f5 = exp(-t).*StepFct(t,0);



f6 = (3-t).*(StepFct(t,2) - StepFct(t,3));


f7 = StepFct(t,0) - StepFct(t,1) - StepFct(t,1) + StepFct(t,2);

f8 = t.*(StepFct(t,0) - StepFct(t,1)) + StepFct(t,1) - StepFct(t,2);


f9 = (1-t).* (StepFct(t,0) - StepFct(t,2));

f10 = (2*t).* (StepFct(t,0) - StepFct(t,0.5)) + StepFct(t,0.5) - StepFct(t,1.5) + (4-2*t).*(StepFct(t,1.5)- StepFct(t,2));

f11 = (StepFct(t,0) - 2*(StepFct(t,0.5)) + 2*(StepFct(t,1.5)) - StepFct(t,2));

f12 = t.*(StepFct(t,0) - StepFct(t,1)) - t.*(StepFct(t,1) - StepFct(t,2));
    plot(t,f12);

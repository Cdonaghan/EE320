
%%state

a = [2,2/3, 1/3; 1, 5/3, 1/3;  0, 2/3, 7/3];

b = [2/3, 0; -1/3, -1; 2/3, 0];

c = [3 1 2];




[T newa] = eig(a);

newb = inv(T)*b;

%%obs - statE 2 is not observable

newC = c*T;

%%controllability - state 2 is not controllable

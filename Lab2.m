
%%state

a = [5, 1, -2, 0; 1, 5, 0, -2; -2, 0, 5, 1; 0, -2, 1, 5];

b = [2, 1, 1; 2, -1, -3; 0, -3, -1; 0, -1, 3];

c = [0, 0, 1, 0; 0, 0, 0, 1];




[T newa] = eig(a);

newb = inv(T)*b;

%%obs - every state is observable

newC = c*T;

%%controllability - state one is not controllable





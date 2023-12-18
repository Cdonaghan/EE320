t = (0:0.01:10);
x1 = sin(2*pi*2*t);
plot(t,x1);

x2 = 0.5*sin(2*pi*4*t);
y = x1 + x2;
plot(t,y)
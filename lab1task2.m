t = 0:0.1:10;
f = 0.5;
w=figure;

%%x =  sin(2*pi*f*t);
%%plot(t,x);
x6=zeros(1,length(t));


for v = 1:+1:5
  
x4 = 2*v-1;
x2 = 1/x4;
x3 = sin(2*pi*x4*f*t);
x6 = x6+(x2*x3);
figure;plot(t,x6,'r--o','Linewidth',2);


%%stem(t,x,'ro');
fo='Fourier no .%d ';
grid on;
str=sprintf(fo,v);
title(str);

    



end



%%last part
p = 0.2*randn(size(t));
b = p + x;
plot(t,b);

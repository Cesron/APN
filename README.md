syms y t;
disp('Método Rugen Kuttan 4');
f=input('Ingrese la ecuación diferencial dy/dt: ');
intervalo=input('Ingrese el intervalo [a,b]');
y0=input('Ingrese el valor inicial y0: ');
F=input('Ingrese la solución de la ecuación diferencial: ');
h=input('Ingrese el valor de h: ');
a=intervalo(1);
b=intervalo(2);
n=int16((b-a)/h);
T=[a:h:b];
Y(1)=y0;
fprintf(' ti\t\t\t\t || k1\t\t\t\t || k2\t\t\t\t || k3\t\t\t\t || k4\t\t\t\t || Yi\t\t\t\t || F(ti)\t\t\t\t || error\n')
fprintf(' %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %e\n',T(1),0,0,0,0,double(Y(1)),double(subs(F,t,T(1))),0)
a1=1/6;
a2=1/3;
a3=1/3;
a4=1/6;
b1=1/2;
b2=1/2;
b3=1;
c11=1/2;
c21=0;
c22=1/2;
c31=0;
c32=0;
c33=1;
for i=1:n
   k1=subs(f,{t,y},{T(i),Y(i)});
   k2=subs(f,{t,y},{T(i)+b1*h,Y(i)+c11*h*k1});
   k3=subs(f,{t,y},{T(i)+b2*h,Y(i)+c21*h*k1+c22*h*k2});
   k4=subs(f,{t,y},{T(i)+b3*h,Y(i)+c31*h*k1+c32*h*k2+c33*h*k3});
   Y(i+1)=Y(i)+h*(a1*k1+a2*k2+a3*k3+a4*k4);
   exacta=subs(F,T(i+1));
   error=abs(exacta-Y(i+1));
   fprintf(' %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %.15f\t\t\t\t || %e\n',T(i+1),double(k1),double(k2),double(k3),double(k4),double(Y(i+1)),double(exacta),double(error))
end
fprintf('Valor aproximado Y(%.15f)= %.15f\n',b,double(Y(n+1)));
fprintf('Valor exacto F(%.15f)= %.15f\n',b,double(exacta));
fprintf('Error= %e\n', double(error));


% Sistema de ecuaciones
syms t x y;
A=24*exp(-2*t)-27*exp(-t)+4*cos(t)-5*sin(t);
B=-2*exp(-2*t)+18*exp(-t)+7*cos(t)-8*sin(t);
Dx=(-4*diff(B)-7*B+7*diff(A)-9*A-101*x)/-88;
Dy=(7*diff(A)-4*diff(B)-8*A+5*B+101*y)/88;
[X,Y]=dsolve('Dx=(101*x)/88 - (189*exp(-t))/44 + (277*exp(-2*t))/44 + cos(t) - (101*sin(t))/88', 'Dy=(101*y)/88 + (567*exp(-t))/88 - (277*exp(-2*t))/44','x(0)=0','y(0)=-1');
T=[0:1/200:1/40];
[T' double(subs(X,T))' double(subs(Y,T))']


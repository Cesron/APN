% Cuadratura Gaussiana
syms x
f=(pi/24*x+pi/8)^2*sin(pi/24*x+pi/8)*pi/24
p=x^6-(15/11)*x^4+5/11*x^2-5/231
n=6
X=double(solve(p))
Y=double(subs(f,X))
h=diff(p)
for i=1:n
d=int(p/(x-X(i)))
m=subs(d,1)-subs(d,-1)
W(i)=double(norm(1/subs(h,X(i))*m))
end
double(W)
integral=0
for i=1:n
integral=integral+W(i)*Y(i)
end
vpa(integral)
clear
% Fin cuadratura gaussiana

% Richardson
syms x;
f=x*log10(x+2)-x^2*acos(exp(2*x))
d=diff(f);
h=1/100
h2=h;
c=-0.65
m=6
fprintf('\n');
N=zeros(m);
for i=1:m
    N(i,1)=(subs(f,c+h)-subs(f,c-h))/(2*h);
    fprintf('(1/2*%.7f)*(f(%.7f+(%.7f))-f(%.7f-(%.7f)))\n',double(h),double(c),double(h),double(c),double(h));
    h=h/2;
end
for j=2:m
    for i=1:m+1-j       
        N(i,j)=(4^(j-1)*N(i+1,j-1)-N(i,j-1))/(4^(j-1)-1);
        fprintf('(%d*N%d%d - N%d%d)/%d\n',double(4^(j-1)),double(i+1),double(j-1),double(i),double(j-1),double(4^(j-1)-1));
    end
end
a=1;
for i=1:m  
    for j=1:m+1-i
        fprintf('N(%d)(h/%d)= %.15f   ',j,a,N(i,j));
    end
    fprintf('\n');
    a=a*2;
end
fprintf('\nValor aproximado: f�(c)= %.15f',N(1,m));
exacto=double(subs(d,c));
fprintf('\nValor exacto: f�(c)= %.15f',exacto);
fprintf('\nError: %e\n',abs(exacto-N(1,m)));
% Fin de Richardson

% Romberg
clear
clc
syms x;
f=sin(x)
a=0
b=pi/4
tam=3
h1=b-a;
h(1)=h1;
for i=2:tam
h(i)=h(i-1)/2;
end
R=zeros(tam);
R(1,1)=h(1)/2*(subs(f,a)+subs(f,b));
for i=2:tam
k=0;
for n=1:2^(i-2)
k=k+subs(f,a+(2*n-1)*h(i));
end
R(i,1)=1/2*(R(i-1,1)+h(i-1)*k);
end
for i=2:tam
for j=2:i
R(i,j)=(4^(j-1)*R(i,j-1)-R(i-1,j-1))/(4^(j-1)-1);
end
end
for i=1:tam
for j=1:i
fprintf('R(%d,%d)= %.15f',i,j,R(i,j));
fprintf('\n');
end
fprintf('\n');
end
fprintf('\n');
fprintf('\nEl valor aproximado de la funci�n f(x) es: %.15f',R(tam,tam));
exacto=double(int(f,a,b));
fprintf('\nEl valor exacto de la integral de f(x) es: %.15f',exacto);
error=abs(exacto-R(tam,tam));
fprintf('\nEl error en la aproximacion es: %e\n',error);
% Fin Romberg

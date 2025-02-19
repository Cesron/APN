disp('------------- Prueba de Bisección ------------- '); 
syms x;
f=input('Introduzca la función f(x): ');
a=input('Introduzca el punto a: ');
b=input('Introduzca el punto b: ');
tol=input('Introduzca el margen de error: 10^-');
tol=10^-tol;
fa=subs(f,a);
fb=subs(f,b);
if fa*fb<0;
    count=1;
    c=(a+b)/2;
    fc=subs(f,c);
    error=abs(fc);
    fprintf('n  || a\t\t\t\t\t|| b\t\t\t\t || c\t\t\t\t  || error\n');
    fprintf('%02d || %.15f || %.15f || %.15f || %e\n', count, double(a), double(b), double(c), double(error));
    while error>tol;
        count=count+1;
        if fa*fc<0;
            b=c;
            c=(a+b)/2;
            error=abs(b-c);
        else
            a=c;
            c=(a+b)/2;
            error=abs(a-c);
        end
        fa=subs(f,a);
        fb=subs(f,b);
        fc=subs(f,c);
        fprintf('%02d || %.15f || %.15f || %.15f || %e\n', count, a, b, c, error);
    end
    fprintf('\nEl valor aproximado de x es: %.15f\n', c);
end






disp('--------------- Prueba punto fijo ---------------');
syms x;
g=input('Introduzca la función g(x): ');
x0=input('Introduzca el valor de x0: ');
tol=input('Introduzca el error: 10^-');
tol=10^-tol;
x1=subs(g,x0);
error=abs(x0-x1);
count=1;
fprintf('n  || x0\t\t\t\t|| x1\t\t\t\t || error\n')
fprintf('%02d || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(error))
while error>tol;
    count=count+1;
    x0=x1;
    x1=subs(g,x0);
    error=abs(x0-x1);
    fprintf('%02d || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(error))
end
fprintf('\nEl valor aproximado de x es: %.15f\n', double(x1))





disp('--------------- Prueba Newton Raphson ---------------');
syms x;
f=input('Introduzca la función f(x): ');
x0=input('Introduzca el primer de x0: ');
tol=input('Introduzca el margen de error: 10^-');
tol=10^-tol;
df=diff(f);
x1=x0-subs(f,x0)/subs(df,x0);
error=abs(x0-x1);
count=1;
fprintf('n  || x0\t\t\t  || x1\t\t\t || error\n');
fprintf('%02d || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(error));
while error>tol;
    count=count+1;
    x0=x1;
    x1=x0-subs(f,x0)/subs(df,x0);
    error=abs(x0-x1);
    fprintf('%02d || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(error));
end
fprintf('\nEl valor aproximado de x es: %.15f\n', double(x1));






disp('--------------- Prueba Secante ---------------');
syms x;
f=input('Introduzca la función f(x): ');
x0=input('Introduzca el valor de x0: ');
x1=input('Introduzca el valor de x1: ');
tol=input('Introduzca el margen de error: 10^-');
tol=10^-tol;
x2=x1-(subs(f,x1)*(x1-x0))/(subs(f,x1)-subs(f,x0));
error=abs(x2-x1);
count=1;
fprintf('n  || x0\t\t\t\t|| x1\t\t\t\t || x2\t\t\t\t  || error\n');
fprintf('%02d || %.15f || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(x2), double(error));
while error>tol
    count=count+1;
    x0=x1;
    x1=x2;
    x2=x1-(subs(f,x1)*(x1-x0))/(subs(f,x1)-subs(f,x0));
    error=abs(x2-x1);
    fprintf('%02d || %.15f || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(x2), double(error));
end
fprintf('\nEl valor aproximado de x es: %.15f\n', double(x2));






disp('----------------------- Método de Posición Falsa -----------------------');
syms x;
f=input('Introduzca la función f(x): ');
x0=input('Introduca el valor de x0: ');
x1=input('Introduzca el valor de x1: ');
tol=input('Introduzca el margen de error: 10^-');
tol=10^-tol;
x2=x1-(subs(f,x1)*(x1-x0))/(subs(f,x1)-subs(f,x0));
error=abs(x2-x1);
count=1;
fprintf('n  || x0\t\t\t\t || x1\t\t\t\t || x2\t\t\t\t || error\n');
fprintf('%02d || %.15f || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(x2), double(error));
while error>tol;
    count=count+1;
    if subs(f,x0)*subs(f,x2) < 0;
        x1=x2;
        x2=x1-(subs(f,x1)*(x1-x0))/(subs(f,x1)-subs(f,x0));
        error=abs(x1-x2);
    else
        x0=x2;
        x2=x1-(subs(f,x1)*(x1-x0))/(subs(f,x1)-subs(f,x0));
        error=abs(x0-x2);
    end
    fprintf('%02d || %.15f || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(x2), double(error));
end
fprintf('\nEl valor aproximado de x es: %.15f\n', double(x2));



disp('-------------------- Método de Steffensen --------------------');
syms x;
g=input('Introduzca la función g(x): ');
x0=input('Introduzca el valor de x0: ');
tol=input('Introduzca el margen de error: 10^-');
tol=10^-tol;
x1=subs(g,x0);
x2=subs(g,x1);
x3=x0-((x1-x0)^2)/(x2-2*x1+x0);
error=abs(x3-x0);
count=1;
fprintf('n  || x0\t\t\t\t|| x1\t\t\t\t || x2\t\t\t\t  || x3\t\t\t\t   || error\n');
fprintf('%02d || %.15f || %.15f || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(x2), double(x3), double(error));
while error>tol;
    count=count+1;
    x0=x3;
    x1=subs(g,x0);
    x2=subs(g,x1);
    x3=x0-((x1-x0)^2)/(x2-2*x1+x0);
    error=abs(x3-x0);
    fprintf('%02d || %.15f || %.15f || %.15f || %.15f || %e\n', count, double(x0), double(x1), double(x2), double(x3), double(error));
end
fprintf('\nEl valor aproximado de x es: %.15f\n', double(x3));




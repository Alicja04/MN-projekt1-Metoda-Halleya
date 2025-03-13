%Poniżej znajduje się kod generujący tabele przydatne do analizy zbieżności
%ciekawych wielomianów
%Przykład 1
ciekawe_x0 = [0, -sqrt(3/2), sqrt(3/2),-sqrt(1/2), sqrt(1/2), 0.7, 0.9, 1, 10];
a = [ 8     0   -12     0];
alfa = [0 sqrt(3/2) -sqrt(3/2)];
eps = 2e-16;
T1 = convergace_table(a, alfa, ciekawe_x0, eps);
%Przykład 2
alfa = [1 2 3 4 -1 -2 -3 -4];
ciekawe_x0 = [alfa 0 1.5026 2.5903 3.6787 0.77068 2.1180 3.3154 0.7 0.9 10];
a = [1 0 -30 0 273 0 -820 0 576];
T2 = convergace_table(a, alfa, ciekawe_x0, eps);

%Przykład 3
a = [1 0 -2 2];
real_alfa = -(nthroot(9 - sqrt(57), 3)/3^(2/3)) - 2/nthroot(3*(9-sqrt(57)), 3);
A = roots(a);
alfa = [real_alfa, A(2), A(3)];
ciekawe_x0 = [alfa -sqrt(2/3) sqrt(2/3) 0, real(A(2)) real(A(3)) abs(A(2)) abs(A(3))];
T3 = convergace_table(a, alfa, ciekawe_x0, eps);

%Przykład 4
a = [1 -3.003 3.006 -1.003];
alfa = [1, 1.001, 1.002];
ciekawe_x0 = [alfa 1.00042 1.00158 10 8 2 4 1.003 -9.6];
T4 = convergace_table(a, alfa, ciekawe_x0, eps);

%Przykład 5
a = [1 2 -11 -12 36];
alfa = [-3, 2];
xc1 = -0.5 -(5/(2*sqrt(3)));
xc2 = (5/(2*sqrt(3))) -0.5;
ciekawe_x0 = [alfa -1/2 xc1 xc2 4 -4 10 2.0001 -3.0001];
T5 = convergace_table(a, alfa, ciekawe_x0, eps);
%Przykład 6
a = [1 -6 11 -6];
alfa = [1 2 3];
xc1 = 2 -1/sqrt(3);
xc2 = 2 +1/sqrt(3);
ciekawe_x0 = [alfa xc1 xc2 0 0.5 10 -2 -3 -0.001 1.9];
T6 = convergace_table(a, alfa, ciekawe_x0, eps);
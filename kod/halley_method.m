function [xn, iter, convergence] = halley_method(a, x0, max_iter, eps)
%halley_method Implementacja metody Halley'a do znajdowania pierwiastków wielomianu
%
%   Dane wejściowe:
%   a - Wektor współczynników wielomianu (np. dla x^3 + 2x^2 - 5 -> a = [1, 2, 0, -5])
%   x0 -  Punkt początkowy iteracji
%   max_iter - Maksymalna liczba iteracji
%   eps - Tolerancja błędu między kolejnymi przybliżeniami (opcjonalne, domyślnie 2e-16)
%
%   Dane wyjściowe:
%   xn- Przybliżony pierwiastek uzyskany metodą Halley'a
%   iter - Liczba iteracji potrzebnych do znalezienia rozwiązania
%   convergence - Flaga logiczna (true/false), czy metoda zbiegła do rozwiązania

    if nargin < 4
        eps=2e-16;
    end
    if length(x0) > 1
        return
    end
    %inicjalizacja zmiennych
    convergence=false;
    xn=x0;
    iter=0;
    while iter<max_iter
        F=horner(xn,a);
        f=F(:,2);
        f1=F(:,3);
        f2=F(:,4);
        %sytuacja dzielenia przez 0
        if abs(f1) < 2e-16
            break;
        end
        div=f1 - f2*(f/(2*f1));
        if abs(div) < 2e-16
            break;
        end
        %następne przybliżenie miejsca zerowego
        x_next=xn-f/div;   
        iter=iter+1;
        if abs(x_next-xn)<eps
            convergence=true;
            xn = x_next;
            break;
        end
        xn=x_next;
    end
end
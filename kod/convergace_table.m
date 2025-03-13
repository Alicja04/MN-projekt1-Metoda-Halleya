function [T] = convergace_table(a, alfa, x0, eps)
%convergace_table Zwraca tabelę pozwalającą na szeroką analizę dokładności
%                 wyznaczianych pierwiastków wielomianów za pomocą metody Halley'a implementowanej funkcją 'halley_method' 
%
%   Dane wejściowe:
%   alfa      - Wektor pierwiastków, do których porównywane są obliczone pierwiastki
%   a         - Wektor współczynników wielomianu przekazywany do funkcji 'halley_method2'
%   x0        - Wektor punktów początkowych dla funkcji 'halley_method'
%   eps       - Tolerancja błędu między kolejnymi przybliżeniami przekazywana do funkcji 'halley_method' (opcjonalne, domyślnie 2e-16)
%
%   Dane wyjściowe:
%   T         - Tablica, której kolejne wiersze to:
%               * x0 - przybliżenie początkowe
%               * xn - wyznaczone miejsce zerowe
%               * Liczba iteracji zwracana przez funkcję 'halley_method'
%               * Dokładny pierwiastek wielomianu (jeśli wielomian ma kilka pierwiastków, to podoawany ten najbliżej wyznaczonego xn)
%               * Różnica pomiędzy dokładnym pierwiastkiem a xn
%               * Kolejne trzy kolumny to wartość funkcji w punkcie x0 oraz wartość 1 i 2 pochodnej też w tym punkcie
   
    if nargin < 4
        eps=2e-16;
    end
    solution = ones(length(x0),8);
    for i=(1:length(x0))
        [xn, iter, ~] = halley_method(a, x0(i), 100, eps);
        solution(i, 2) = xn;
        solution(i, 3) = iter;
        solution(i, 1) = x0(i);
        [~, idx] = min(abs(alfa-xn));
        solution(i, 4) = alfa(idx);
        solution(i,5) = abs(alfa(idx) - xn);
        F = horner(x0(i), a);
        solution(i, 6) = F(:,2);
        solution(i, 7) = F(:,3);
        solution(i, 8) = F(:,4);
    end
    T = array2table(solution, 'VariableNames', {'x0 - przybliżenie początkowe', 'xn - wyznaczone miejsce zerowe', 'Liczba iteracji', 'Dokładny root', 'Różnica pomiędzy dokładnym root a xn', 'f(x0)', 'f1(x0)', 'f2(x0)'});   
end

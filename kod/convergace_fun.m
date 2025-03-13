function [x0, y] = convergace_fun(alfa, a, x01, x0n, step, max_iter, epsConvergace, epsHalley)
%convergace_fun Oblicza ilość iteracji metody Halley'a dla różnych początkowych przybliżeń
%
%   Dane wejściowe:
%   alfa      - Wektor pierwiastków, do których porównywane są obliczone pierwiastki
%   a         - Wektor współczynników wielomianu przekazywany do funkcji 'halley_method2'
%   x01       - Punkt początkowy dla zakresu początkowych przybliżeń
%   x0n       - Punkt końcowy dla zakresu początkowych przybliżeń
%   step      - Krok między kolejnymi początkowymi przybliżeniami
%   max_iter  - Maksymalna liczba iteracji dla metody Halley'a
%   epsConvergace - Tolerancja błędu między obliczonymi pierwiastkami a wartościami w 'alfa' (opcjonalne, domyślnie 2e-16)
%   epsHalley - Toleracja przekazywana do funkcji 'halley_method' (opcjonalne, domyślnie 2e-16)
%
%   Dane wyjściowe:
%   x0        - Wektor początkowych przybliżeń (od x01 do x0n)
%   y         - Wektor liczb iteracji dla każdego początkowego przybliżenia. 
%               Wartość 'y(i)' to liczba iteracji, jeśli metoda Halley'a zbiegnie do pierwiastka 
%               w 'alfa' w granicach 'eps', lub -1, jeśli nie osiągnięto zbieżności

    if nargin < 7
        epsConvergace=2e-16;
    end
    if nargin < 8
        epsHalley = 2e-16;
    end
    round_point = length (num2str(step - fix(step))) - 2;
    x0 = x01 : step : x0n;
    y = -ones(1, length(x0));
    roots_all = NaN(length(x0), 2);
    convergance_all = false(1, length(x0));
    for i=(1:length(x0))
        [xn, iter, convergance] = halley_method(a, round(x0(i), round_point), max_iter, epsHalley);
        roots_all(i, :) = [xn, iter];
        convergance_all(i) = convergance;
    end

    for i=(1:length(x0))
        [~, idx] = min(abs(alfa - roots_all(i,1)));
        if abs(alfa(idx) - roots_all(i, 1)) <= epsConvergace
            y(i) = roots_all(i, 2);
        end
    end
end
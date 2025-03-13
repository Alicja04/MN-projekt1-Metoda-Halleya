function F = horner(x, a)
%horner : Funkacja wylicza wartości wielomianu oraz dwóch pierwszych pochodnych (f, f', f") dla punktów x
%
%   Dane wejściowe:
%   x - Wektor punktów dla któych liczymy wartości;
%   a - Współczynniki wielomianu (od najwyższej do najniższej potęgi);
%
%   Dane wyjściowe:
%   F- Macierz której kolejne kolumny to x, f, f1, f2 przedstawiające:
%       *x - punkty podane na wejściu
%       *f - wartości wielomianu w kolejnych punktach x
%       *f1 - wartości pierwszej pochodnej w kolejnych punktach x
%       *f2 - wartości drugiej pochodnej w kolejnych punktach x
    x=x(:);
    a=reshape(a',1,[]);
    n=length(a);
    f=a(1);
    f1=f;
    f2=f;
    if n==2
        f2=0*ones(length(x), 1);
        f=0;
        f1=0;
    end
    if n==1
        f1=0;
        f2=0*ones(length(x), 1);
    end
    for i=2:(n-2)
        f=a(i)+x.*f; 
        f1=f+x.*f1;
        f2=f1+x.*f2;
    end
    f=a(n-1)+x.*f;
    f1=f+x.*f1;
    f=a(n)+x.*f;
    f2=2*f2;
    F=[x,f,f1,f2];
end
function [xk] = Homeier (a, x0)
%% Wejście
% a - wektor kolejnych wartości współczynników wielomianu zaczynając od an
% czyli pierwszy wyraz wektora to an, drugi to an-1, a ostatni to a0
% x - punkt, w którym chcemy otrzymać wartość funkcji
%podaje użytkownikowi wielomiany przynajmniej 3 stopnia bo w niższych łatwo
%wyznaczyć zera wielomianów
%% Wyjście
%funkcja zwraca miejsce zerowe jeśli zachodzi warunek stopu lub po 1000
%liczbie iteracji zwraca, że metoda nie zbiega, ponadto zwraca mi liczbę k,
%która pokazuje ile było potrzebnych iteracji do osiągnięcia przybliżenia
%%
k = 0;
horn = zeros(1,3);
hornpom = zeros(1,3);
x = x0;
xk = zeros(1,2);
blad = 1;
while blad > 0,0001;
    if k > 1000
        return
    end
    xold = x;
    horn = Hornerf(a,x);    %horn[f(x), f'(x), f''(x)]
    y = x - horn(1)/horn(2);
    hornpom = Hornerf(a,y);
    x = x - 1/2*horn(1)*(1/(horn(2)) + 1/(hornpom(2)));
    blad = abs(xold - x);
    xk(1) = x;
    k = k+1;
    xk(2) = k;
end




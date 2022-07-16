function [xk] = Halley (a, x0)
%% Wejście
% a - wektor kolejnych wartości współczynników wielomianu zaczynając od an
% czyli pierwszy wyraz wektora to an, drugi to an-1, a ostatni to a0
% x - punkt, w którym chcemy otrzymać wartość funkcji
%podaje użytkownikowi wielomiany przynajmniej 3 stopnia bo w niższych łatwo
%wyznaczyć zera wielomianów
%% Wyjście
%funkcja zwraca miejsce zerowe jeśli zachodzi warunek stopu lub po 1000
%liczbie iteracji zwraca, że metoda nie zbiega, ponadto zwraca ile było
%potrzebnych iteracji do osiągnięcia miejsca zerowego
%%
k = 0;
hall = zeros(1,3);
x = x0;
xk = zeros(1,2);
blad = 1;
while blad > 0,0001;
    if k > 1000
        return
    end
    xold = x;
    hall = Hornerf(a,x);
    x = x - (2*hall(1)*hall(2))/(2*(hall(2)^2) - hall(1)* hall(3));
    blad = abs(x - xold);
    xk(1) = x;    
    k = k+1;
    xk(2) = k;
end







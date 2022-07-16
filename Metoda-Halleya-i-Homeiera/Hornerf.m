function [f] = Hornerf (a, x)
%% Wejście
% a - wektor kolejnych wartości współczynników wielomianu zaczynając od an
% czyli pierwszy wyraz wektora to an, drugi to an-1, a ostatni to a0
% x - punkt, w którym chcemy otrzymać wartość funkcji
%podaje użytkownikowi wielomiany przynajmniej 3 stopnia bo w niższych łatwo
%wyznaczyć zera wielomianów
%% Wyjście
% zwracamy wektor [f, f', f''], z którego będziemy korzystać przy
% wyznaczaniu miejsc zerowych w metodzie Halleya i Homeyera
%%
n = length(a);
w = a(1);
p = w;
r = p;
for k = 2:n-2
    w = a(k) + x*w;
    p = w + x*p;
    r = p + x*r;
end
w = a(n-1) + x*w;
p = w + x*p;
w = a(n) + x*w;
f = zeros(1,3);
f(1) = w;
f(2) = p;
f(3) = 2*r;

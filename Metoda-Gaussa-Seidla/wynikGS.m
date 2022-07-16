%% METODA GAUSSA SEIDLA %%
% warunki wystarczające zbieżności metody G-S dla równań AX = B to:
% 1) A jest macierzą dominującą wierszowo (aii większy niż suma pozostaych
% wyrazów w danym wierszu)
% 2) A jest macierzą dominującą kolumnowo (analogicznie dla kolumn
% z tego powodu często będę uwzględniał takie macierze w testach
% programistycznych

% Najpier pokażę na podsawie przykładu macierzy zbieżnej że mój kod jest
% poprawny
A0 = [7 1 1; 2 10 3; 1 1 8];
B0 = [2 3 4; 7 8 1; 1 3 1];
[Xa, ka] = GaussSeidel(A0,B0,100);
[Xb, kb] = GSodwr(A0,B0,100);
Za = Xa - inv(A0)*B0;          %powinno zwracać macierz wartości bardzo bliskich zeru
Zb = Xb - B0*inv(A0);          % wektory ka oraz kb są zerowe od pewnego momentu

%% Przykłady macierzy w których metoda jest szybko zbieżna, wolno zbieżna i rozbieżna 
% wynik przedstawie za pomocą wykresu wartości błędu od numeru kroku (będę
% brał co 10 krok, aby było widać różnice pomiędzy wartościami
B4 = randi([0,10], [3,3]);

% zbieżna szybko
A4 = [18 1 1; 3 54 1; 2 5 64];
[X4a, k4a] = GaussSeidel(A4,B4,500);
[X4b, k4b] = GSodwr(A4,B4,500);

%zbieżna wolno
A5 = [18 11 7; 53 54 1; 12 51 64];
[X5a, k5a] = GaussSeidel(A5,B4,500);
[X5b, k5b] = GSodwr(A5,B4,500);

%zbieżna bardzo wolno
A6 = [18 11 9; 53 54 51; 15 51 64];
[X6a, k6a] = GaussSeidel(A6,B4,500);
[X6b, k6b] = GSodwr(A6,B4,500);

wyk4a = zeros(100);
wyk5a = zeros(100);
wyk6a = zeros(100);
wyk4b = zeros(100);
wyk5b = zeros(100);
wyk6b = zeros(100);
for i=1:100
   wyk4a(i) = k4a(5*i);
   wyk4b(i) = k4b(5*i);
   wyk5a(i) = k5a(5*i);
   wyk5b(i) = k5b(5*i);
   wyk6a(i) = k6a(5*i);
   wyk6b(i) = k6b(5*i);
end
x= (1:100);
hold on
wykres2 =subplot(2,2,1);
wykres2 = semilogy(x,wyk4a,'-*r', x, wyk4b,'-*g', x, wyk5a, '-+r', x, wyk5b, '-+g', x, wyk6a, '--r', x, wyk6b, '--g');
title('Log wartości błędu dla macierzy zbieżnych');
xlabel('k');
ylabel('Log błędu');
hold off


%rozbieżna powoli
A7 = [18 11 12; 53 54 63.4; 15 51 64];
[X7a, k7a] = GaussSeidel(A7,B4,500);
[X7b, k7b] = GSodwr(A7,B4,500);

%rozbieżna szybko
A8 = [18 11 12; 53 54 86; 15 51 64];
[X8a, k8a] = GaussSeidel(A8,B4,500);
[X8b, k8b] = GSodwr(A8,B4,500);

%rozbieżna bardzo szybko
A9 = [18 11 12; 53 54 786; 15 51 64];
[X9a, k9a] = GaussSeidel(A9,B4,500);
[X9b, k9b] = GSodwr(A9,B4,500);

for i=1:100
   wyk4a(i) = k7a(5*i);
   wyk4b(i) = k7b(5*i);
   wyk5a(i) = k8a(5*i);
   wyk5b(i) = k8b(5*i);
   wyk6a(i) = k9a(5*i);
   wyk6b(i) = k9b(5*i);
end
x= (1:100);
hold on
wykres3 =subplot(2,2,2);
wykres3 = semilogy(x,wyk4a,'-*r', x, wyk4b,'-*g', x, wyk5a, '-+r', x, wyk5b, '-+g', x, wyk6a, '--r', x, wyk6b, '--g');
title('Log wartości błędu dla macierzy rozbieżnych');
xlabel('k');
ylabel('Log błędu');
hold off
% ostatnia macierz jest tak szybko zbieżna że nawet wykres logarytmiczny
% nie pokazuje czytelnego wykresu

%% MACIERZ JEDNOSTKOWA
% Teraz przeprowadzę kilka testów dla macierzy jednostkowej. Jest to
% interesujący przykład, ponieważ wynik działań AX=B oraz XA = B powinien
% być identyczny, jednakże metoda G-S nie zawsze zwróci dobry wynik

%% macierz dominująca rzędowo

A1 = [4 1 3; 22 253 215; 1 189 191];
B1 = [1 0 0; 0 1 0; 0 0 1];
[X1a, k1a] = GaussSeidel(A1,B1,100);
[X1b, k1b] = GSodwr(A1,B1,100);
Z1a = abs(X1a - X1b);
% tutaj wynik jest inny jednak różnica jest rzędu 10^(-10) ponieważ wektor
% na tej macierzy metoda wolno zbiega
[X1c, k1c] = GaussSeidel(A1,B1,1000);
[X1d, k1d] = GSodwr(A1,B1,1000);
Z1c = X1c - X1d;
% przy tysiącu iteracji wynik już ma różnicę rzędu 10^(-16)

%jednakże pokażę przykład że wynik może się dużo bardziej różnić:
%góno prawda xd
% A2 = [ 30 -29 29;69 -70 69; 999 1000 1001];
% B2 = [1 0 0; 0 1 0; 0 0 1];
% [X2a, k2a] = GaussSeidel(A2,B2,1000);
% [X2b, k2b] = GSodwr(A2,B2,1000);
% Z2a = X2a - X2b

%% macierz dominująca kolumnowo ale nie dominująca wierszowo
A2 = [7 150 69; 2 152 1; 4 1 71];
B2 = [1 0 0; 0 1 0; 0 0 1];
[X2a, k2a] = GaussSeidel(A2,B2,1000);
[X2b, k2b] = GSodwr(A2,B2,1000);
Z2a = abs(X2a - X2b);
% tutaj również metoda G-S jest zbieżna w obu przypadkach i wyniki znacząco
% się nie różnią
%% dowolna macierz ale nadal zbieżna w jednym przypadku
A3 = [15 14 10.9; 491 878 697; 1 120 100];
B3 = [1 0 0; 0 1 0; 0 0 1];
[X3a, k3a] = GaussSeidel(A3,B3,10000);
[X3b, k3b] = GSodwr(A3,B3,10000);
Z3a = abs(X3a - X3b);
% tutaj próbowałem manipulować pewnymi liczbami i zauważyłem że jeśli
% metoda jest zbieżna w AX=B to jest również zbieżna w XA=B jednakże są one
% zbieżne w różnym tempie i zwracają różny wynik co widać na macierzy Z3a

%% CZY METODA G-S ZACZYNA BYĆ ROZBIEŻNA DLA RÓWNAĆ AX=B ORAZ XA=B W TYM SAMYM MOMENCIE?
%% CZY MOŻE ONA ZWRÓCIĆ KOMPLETNIE INNE WYNIKI POMIMO ŻE POWINNY ZWRÓCIĆ IDENTYCZNE?
%pokaże manipulując jedną wartością w macierzy, że macierze zaczynają
%rozbiegać w tym samym momencie
katka = zeros(20);
blada = zeros(20);
bladb = zeros(20);
for i=1:20  
    j = 12 - i/10;
    A3 = [15 14 j; 491 878 697; 1 120 100];
    B3 = [1 0 0; 0 1 0; 0 0 1];
    [X3a, k3a] = GaussSeidel(A3,B3,10000);
    [X3b, k3b] = GSodwr(A3,B3,10000);
    dok = inv(A)*B;                   %wywołuje tylko jedną dokładną wartość, bo dla B=I AX=XA
    katka(i) = k3a(1000);
    poma = max(abs(dok - X3a),[], "all");
    pomb = max(abs(dok - X3b),[], "all");
    blada(i) = poma;
    bladb(i) = pomb;
end
y= (1:20);
hold on
wykres1 = subplot(2,2,3);
wykres1 = semilogy(y,katka,'-*','DisplayName', 'dokładność przybliżenia');
wykres1 = semilogy(y, blada, '--r', 'DisplayName','błąd przybliżenia AX');
wykres1 = semilogy(y, bladb, '-g', 'DisplayName','błąd przybliżenia XA');
hold off
% zauważmy, że dla pierwszych 6 iteracji jest tak mały, że matlab wyrzuca 0
% i wykres nie obejmuje tych wartości. Widać również że od wartości j
% równej około 10.6 metoda zaczyna być rozbieżna, jednakże w jednym
% przypadku rozbiega szybciej niż w drugim.
% Zatem wniosek jest taki, że dla macierzy A, dla których metoda Gaussa
% Seidla jest rozbieżna, wyniki równań AX oraz XA są inne

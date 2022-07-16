function [x,k] = GSodwr(A,B, c)
%% Wejście
% A: macierz NxN
% B: macierz NxN
% c: ilość iteracji, które chcemy przeprowadzić
%% Wyjście
% Macierz X o rozmiarach NxN
% wektor k maksymalnych błędów w danej iteracji
%%
while size(A,1) ~= size(B,1)
    input('Blad')
    return
end
while size(A,1) ~= size(A,2)
    input('Blad')
    return
end
while size(A,1) ~= size(B,2)
    input('Blad')
    return
end
n = size(A,1);
x = zeros(n,n); %tworzę pustą macierz rozwiązań aby nadać jej odpowiedni wymiar i uzupełniać ją
k= zeros(1,c);                 %wektor rozwiązań analogicznie jak macierz X
itr = 1;
while itr<c+1
    for m=1:length(k)
        x_old = x;
        for l=1:n
            for i=1:n
                sum=0;
                for j=1:l-1                        %ta i kolejna pętla pozwalają wyznaczyć sumę którą odejmujemy od wyrazu B(i,l)
                     sum=sum+A(j,l)*x(i,j);    
                end
                for j=l+1:n
                    sum=sum+A(j,l)*x_old(i,j);
                end
                x(i,l)=(1/A(l,l))*(B(i,l) - sum); 
            end
        end
        itr = itr+1;
        err = max(abs(x_old - x), [], "all");   
        k(1,m)=max(err);
    end
end

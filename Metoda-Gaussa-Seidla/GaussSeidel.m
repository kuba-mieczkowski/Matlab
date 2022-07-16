function [x,k] = GaussSeidel(A,B,c)
%% Wejście
% A: macierz NxN
% B: macierz NxM
% c: ilość iteracji, które chcemy przeprowadzić
%% Wyjście
% Macierz X o rozmiarach NxM
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
x = zeros(size(A,2),size(B,2)); %tworzę pustą macierz rozwiązań aby nadać jej odpowiedni wymiar i uzupełniać ją
k= zeros(1,c);                 %wektor rozwiązań analogicznie jak macierz X
itr = 1;
while itr< c+1 
    for m=1:length(k)
        x_old = x;
        for l=1:size(B,2)
            for i=1:size(A,2)
                sum=0;
                for j=1:i-1                        %ta i kolejna pętla pozwalają wyznaczyć sumę którą odejmujemy od wyrazu B(i,l)
                     sum=sum+A(i,j)*x(j,l);        
                end
                for j=i+1:size(A,2)
                    sum=sum+A(i,j)*x_old(j,l);
                end
                x(i,l)=(1/A(i,i))*(B(i,l) - sum); 
            end
        end
        itr = itr+1;
        err = max(abs(x_old - x), [], "all");
        k(1,m)=max(err);                        %wyciągam największą różnicę między macierzą x w i-tej i i-1-szej iteracji, która pozwala mi określić zbieżność metody
    end                                         % a następnie zapisuje to na wektorze k, który potem będę wykorzystywał do wykresów zbieżności
end




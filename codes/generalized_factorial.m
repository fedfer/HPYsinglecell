% programma per generare una matrice dei logaritmi dei coefficienti fattoriali
% generalizzati partendo dalla nota relazione di ricorrenza. Attenzione:
% al posto (n,l) ci sar? il logaritmo del coefficiente fattoriale generalizzato
% C(n-1,l-1,s), dato che voglio inserire anche i coefficienti con n=0 e
% l=0. 
% ATTENZIONE: se s=1 il programma non funziona, perch? gli elementi della
% diagonale sono nulli e non ne posso fare il logaritmo: in tal caso
% andrebbe vista una soluzione a parte.

function LogC=generalized_factorial(N,L,s)

% C(N,L,s) ? il massimo coefficente fattoriale che voglio generare, quindi
% le dimensioni della matrice saranno (N+1)x(L+1), visto che inserisco
% anche i coefficienti fattoriali con lo 0. Indico con LogC il logaritmo
% del coefficiente fattoriale generalizzato.
LogC=zeros(N+1,L+1);
% la matrice ? triangolare inferiore, visto che se l>n C(n,l;s)=0
% inizializzo C(0,0;s)
LogC(1,1)=0;
% genero i logaritmi dei coefficienti fattoriali C(n,1,s) e quelli sulla
% diagonale perch? nella formula ricorsiva non voglio fare il logaritmo di
% zero
for n=1:N
    LogC(n+1,n+1)=n*log(s);
end
for n=2:N
    LogC(n+1,2)=log(n-1-s)+LogC(n,2);
end
for n=2:N
    for l=3:n
        % C(n+1,l) ? il coefficiente fattoriale C(n,l-1), n ed l sono
        % riferiti alla matrice!! NOn all'ordine del coefficiente
        % fattoriale
        x=log(s+exp(log(n-1-s*(l-1))+LogC(n,l)-LogC(n,l-1)));
        if x==inf
            % in tal caso non passo all'esponenziale e approssimo x, in
            % realt? non perdo nulla essendo s \in (0,1)
            x=log(n-1-s*(l-1))+LogC(n,l)-LogC(n,l-1);
            disp('Approssimazione di x');
        end
        LogC(n+1,l)= LogC(n,l-1)+x;
    end
end

% quindi se serve il coefficiente fattoriale C(n,l;s), bisogna calcolare la
% matrice C in (n+1,l+1)

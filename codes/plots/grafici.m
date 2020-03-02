clear all;
close all;


load clementina_grafici_FRUTTO1.mat
load clementina_grafici_FRUTTO2.mat

figure(1);
title('posterior histogram of \theta_0');
ist=histc(M_parametri(:,1),1:max(M_parametri(:,1)));
bar(ist/N);

figure(2);
title('posterior histogram of \theta');
ist=histc(M_parametri(:,3),1:max(M_parametri(:,3)));
bar(ist/N);

figure(3);
title('posterior histogram of \sigma_0');
ist=histc(M_parametri(:,2),0:0.001:1);
bar(0:0.001:1,ist/N,0.01);

figure(4);
title('posterior histogram of \sigma');
ist=histc(M_parametri(:,4),0:0.001:1);
bar(0:0.001:1,ist/N,0.01);

fprintf('Stime di sigma = %f sigma_0  = %f  theta  = %f theta_0  = %f \n',sum(M_parametri(:,4))/N,...
    sum(M_parametri(:,2))/N,sum(M_parametri(:,3))/N,sum(M_parametri(:,1))/N );

% FLOWER

fprintf('---------------- FIRST LIBRARY ----------------');


for i=1:length(M1)
        % trovo posteriro e intervalli di confidenza per n=m
    m=M1(i);
    figure(i+4);
    plot(Valori_predittiva1{i},P_posterior1{i},'-');
    legend('Posterior distribution of K_{n}^{(m)}: FIRST LIBRARY');
    fprintf('La media a posteriori delle nuove specie nelle prossime %i osservazioni è: %f  \n',m,media1(i));
    fprintf('Intervallo di confidenza al 0.95 per la media a posteriori: %f in (%f , %f)  \n',media1(i),media_quantili1(i,1),media_quantili1(i,2));
    fprintf('La probabilità di osservare una nuova nelle prossima %i osservazione è: %f  \n',m+n1+1,media_pr_nuova1(i));
    fprintf('Intervallo di confidenza al 0.95 per la discovery: %f in (%f , %f)  \n',media_pr_nuova1(i), ...
        media_pr_nuova_quantili1(i,1), media_pr_nuova_quantili1(i,2));
    fprintf('La media a posteriori dei nuovi geni per X1 ma non per X2 nelle prossime %i osservazioni è: %f  \n',m,media_nuove_X1_non_X2(i));
    fprintf('Intervallo di confidenza al 0.95 per questa media: %f in (%f , %f)  \n',media_nuove_X1_non_X2(i), ...
        media_nuove_X1_non_X2_quantili(i,1),  media_nuove_X1_non_X2_quantili(i,2));
    fprintf('La media a posteriori delle nuove specie distinte per X1 ma non per X2 nelle prossime %i osservazioni è: %f  \n',m,media_distinte_nuove_X1_non_X2(i));
    fprintf('Intervallo di confidenza al 0.95 per questa media: %f in (%f , %f)  \n',media_distinte_nuove_X1_non_X2(i), ...
        media_distinte_nuove_X1_non_X2_quantili(i,1),  media_distinte_nuove_X1_non_X2_quantili(i,2));
end

figure(5+length(M1));
hold on;
plot(M1,media_pr_nuova1,'*-');
plot(M1,media_pr_nuova_quantili1(:,1)',':*r');
plot(M1,media_pr_nuova_quantili1(:,2)',':*g');
title('FIRST LIBRARY');
legend('decadimento di osservare una nuova specie','quantili di ordine inferiore','quantili di ordine superiore');


% ROOT

fprintf('---------------- SECOND LIBRARY ----------------');


for i=1:length(M2)
        % trovo posteriro e intervalli di confidenza per n=m
    m=M2(i);
    figure(i+5+length(M1));
    plot(Valori_predittiva2{i},P_posterior2{i},'-');
    legend('Posterior distribution of K_{n}^{(m)}: SECOND LIBRARY');
    fprintf('La media a posteriori delle nuove specie nelle prossime %i osservazioni è: %f  \n',m,media2(i));
    fprintf('Intervallo di confidenza al 0.95 per la media a posteriori: %f in (%f , %f)  \n',media2(i),media_quantili2(i,1),media_quantili2(i,2));
    fprintf('La probabilità di osservare una nuova nelle prossima %i osservazione è: %f  \n',m+n2+1,media_pr_nuova2(i));
    fprintf('Intervallo di confidenza al 0.95 per la discovery: %f in (%f , %f)  \n',media_pr_nuova2(i), ...
        media_pr_nuova_quantili2(i,1), media_pr_nuova_quantili2(i,2));
    fprintf('La media a posteriori delle nuove osservazioni per X2 ma non per X1 nelle prossime %i osservazioni è: %f  \n',m,media_nuove_X2_non_X1(i));
    fprintf('Intervallo di confidenza al 0.95 per questa media: %f in (%f , %f)  \n',media_nuove_X2_non_X1(i), ...
        media_nuove_X2_non_X1_quantili(i,1),  media_nuove_X2_non_X1_quantili(i,2));
    fprintf('La media a posteriori delle nuove specie distinte per X2 ma non per X1 nelle prossime %i osservazioni è: %f  \n',m,media_distinte_nuove_X2_non_X1(i));
    fprintf('Intervallo di confidenza al 0.95 per questa media: %f in (%f , %f)  \n',media_distinte_nuove_X2_non_X1(i), ...
        media_distinte_nuove_X2_non_X1_quantili(i,1),  media_distinte_nuove_X2_non_X1_quantili(i,2));
end

figure(6+length(M2)+length(M1));
hold on;
plot(M2,media_pr_nuova2,'*-');
plot(M2,media_pr_nuova_quantili2(:,1)',':*r');
plot(M2,media_pr_nuova_quantili2(:,2)',':*g');
title('SECOND LIBRARY');
legend('decadimento di osservare una nuova specie','quantili di ordine inferiore','quantili di ordine superiore');
% read results

Runs=40;
addsample=100;

% Dati finali
M=zeros(Runs,addsample);
DATAfinal_HPY_big = M;
DATAfinal_uniform_big = M;
DATAfinal_Oracle_big = M;
DATAfinal_GoodTulming_big = M;


for III=1:Runs
    
    load_path = strcat('/Users/felpo/Desktop/resultsHPY/results_',string(III),'.mat');
    load(load_path)  
    
    DATAfinal_HPY_big(III,:)=DATAfinal_HPY(III,:);
    DATAfinal_uniform_big(III,:)=DATAfinal_uniform(III,:);
    DATAfinal_Oracle_big(III,:)=DATAfinal_Oracle(III,:);
    DATAfinal_GoodTulming_big(III,:)=DATAfinal_GoodTulming(III,:);
    
end

DATAfinal_HPY = DATAfinal_HPY_big;
DATAfinal_uniform = DATAfinal_uniform_big;
DATAfinal_Oracle = DATAfinal_Oracle_big;
DATAfinal_GoodTulming= DATAfinal_GoodTulming_big;
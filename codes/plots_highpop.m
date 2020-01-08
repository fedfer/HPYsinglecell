plot(1:addsample,sum(DATAfinal.HPY)/Runs,'r');
hold on
plot(1:addsample,sum(DATAfinal.uniform)/Runs,'g');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs,'b');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs,'k');
plot(1:addsample,quantile(DATAfinal.HPY,0.025),':r');
plot(1:addsample,quantile(DATAfinal.HPY,0.975),':r');
plot(1:addsample,quantile(DATAfinal.uniform,0.025),':g');
plot(1:addsample,quantile(DATAfinal.uniform,0.975),':g');
plot(1:addsample,quantile(DATAfinal.Oracle,0.025),':b');
plot(1:addsample,quantile(DATAfinal.Oracle,0.975),':b');
plot(1:addsample,quantile(DATAfinal.GoodTulming,0.025),':k');
plot(1:addsample,quantile(DATAfinal.GoodTulming,0.975),':k');
legend('HPY','Uniform','Oracle','GT','Location','NorthWest');

%plot with bands2 with sd
sd.HPY = std(DATAfinal.HPY);
sd.uniform = std(DATAfinal.uniform);
sd.GT = std(DATAfinal.GoodTulming);
sd.Oracle = std(DATAfinal.Oracle);
plot(1:addsample,sum(DATAfinal.HPY)/Runs,'r');
hold on
plot(1:addsample,sum(DATAfinal.uniform)/Runs,'g');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs,'b');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs,'k');
plot(1:addsample,sum(DATAfinal.HPY)/Runs+sd.HPY,':r');
plot(1:addsample,sum(DATAfinal.HPY)/Runs-sd.HPY,':r');
plot(1:addsample,sum(DATAfinal.uniform)/Runs+sd.uniform,':g');
plot(1:addsample,sum(DATAfinal.uniform)/Runs-sd.uniform,':g');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs+sd.GT,':k');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs-sd.GT,':k');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs+sd.HPY,':b');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs-sd.HPY,':b');
legend('HPY','Uniform','Oracle','GT','Location','NorthWest');


%plot with sd
shadedErrorBar(1:addsample,DATAfinal.HPY, {@mean,@std}, 'lineprops', 'r','transparent',true,'patchSaturation',0.075);
hold on
shadedErrorBar(1:addsample,DATAfinal.Oracle, {@mean,@std}, 'lineprops', 'b','transparent',true,'patchSaturation',0.075);
shadedErrorBar(1:addsample,DATAfinal.uniform, {@mean,@std}, 'lineprops', 'g','transparent',true,'patchSaturation',0.075);
shadedErrorBar(1:addsample,DATAfinal.GoodTulming, {@mean,@std}, 'lineprops', 'k','transparent',true,'patchSaturation',0.075);
legend('HPY','Uniform','Oracle','GT','Location','NorthWest');

M=zeros(1,addsample);
Final=struct('uniform',M,'Oracle',M,'GoodTulming',M,'GoodTuring2',M,'GoodTuring3',M);
Final.uniform=sum(DATAfinal.uniform)/Runs;
Final.HPY=sum(DATAfinal.HPY)/Runs;
Final.Oracle=sum(DATAfinal.Oracle)/Runs;
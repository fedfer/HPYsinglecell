
plot(1:addsample,sum(DATAfinal_HPY)/Runs,'r');
hold on
plot(1:addsample,sum(DATAfinal_uniform)/Runs,'g');
plot(1:addsample,sum(DATAfinal_Oracle)/Runs,'b');
plot(1:addsample,sum(DATAfinal_GoodTulming)/Runs,'k');

legend('HPY','Uniform','Oracle','Good-Tulming','Location','NorthWest');

%plot with bands2 with sd
sd.HPY = std(DATAfinal_HPY);
sd.uniform = std(DATAfinal_uniform);
sd.GT = std(DATAfinal_GoodTulming);
sd.Oracle = std(DATAfinal_Oracle);
plot(1:addsample,sum(DATAfinal_HPY)/Runs,'r');
shadedErrorBar(1:addsample,DATAfinal_HPY, {@mean,@std}, 'lineprops', 'r','transparent',true,'patchSaturation',0.075);
hold on
plot(1:addsample,sum(DATAfinal_uniform)/Runs,'g');
shadedErrorBar(1:addsample,DATAfinal_uniform, {@mean,@std}, 'lineprops', 'g','transparent',true,'patchSaturation',0.075);
plot(1:addsample,sum(DATAfinal_Oracle)/Runs,'b');
shadedErrorBar(1:addsample,DATAfinal_Oracle, {@mean,@std}, 'lineprops', 'b','transparent',true,'patchSaturation',0.075);
plot(1:addsample,sum(DATAfinal_GoodTulming)/Runs,'k');
shadedErrorBar(1:addsample,DATAfinal_GoodTulming, {@mean,@std}, 'lineprops', 'k','transparent',true,'patchSaturation',0.075);

legend('HPY','Uniform','Oracle','GT','Location','NorthWest');
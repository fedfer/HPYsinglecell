%plot with sd
plot(1:addsample,sum(DATAfinal.HPY)/Runs,'r');
hold on
plot(1:addsample,sum(DATAfinal.uniform)/Runs,'g');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs,'b');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs,'k');

figure; hold on,
h_1 = shadedErrorBar(1:addsample,DATAfinal.HPY, {@mean,@std}, 'lineprops', 'r','transparent',true,...
    'patchSaturation',0.075,'LineStyle','--');
title('Delayed Abundance case');
%title('Incidence case');
%title('Simulation Experiment');
xlabel('Runs');
ylabel('Number of Species Discovered');
h_2 = shadedErrorBar(1:addsample,DATAfinal.Oracle, {@mean,@std}, 'lineprops', 'b','transparent',true,...
    'patchSaturation',0.075,'LineStyle','-');
h_3 = shadedErrorBar(1:addsample,DATAfinal.uniform, {@mean,@std}, 'lineprops', 'g','transparent',true,...
    'patchSaturation',0.075,'LineStyle','-.');
h_4 = shadedErrorBar(1:addsample,DATAfinal.GoodTulming, {@mean,@std}, 'lineprops', 'k','transparent',true,...
    'patchSaturation',0.075,'LineStyle',':');
legend([h_2.mainLine h_1.mainLine h_4.mainLine h_3.mainLine], ...
    'Oracle','HPY','GT','Uniform','Location','NorthWest');

%legend('HPY','Uniform','Oracle','GT','Location','NorthWest');


%plot with sd
sd.HPY = std(DATAfinal.HPY);
sd.uniform = std(DATAfinal.uniform);
sd.GT = std(DATAfinal.GoodTulming);
sd.Oracle = std(DATAfinal.Oracle);
plot(1:addsample,sum(DATAfinal.HPY)/Runs,'r');
hold on
plot(1:addsample,sum(DATAfinal.uniform)/Runs,'g');
plot(1:addsample,sum(DATAfinal.uniform)/Runs,'lineProps', 'b','LineStyle','--');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs,'b');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs,'k');
plot(1:addsample,sum(DATAfinal.HPY)/Runs+sd.HPY,'--r');
plot(1:addsample,sum(DATAfinal.HPY)/Runs-sd.HPY,'--r');
plot(1:addsample,sum(DATAfinal.uniform)/Runs+sd.uniform,'--g');
plot(1:addsample,sum(DATAfinal.uniform)/Runs-sd.uniform,'--g');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs+sd.GT,'--k');
plot(1:addsample,sum(DATAfinal.GoodTulming)/Runs-sd.GT,'--k');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs+sd.HPY,'--b');
plot(1:addsample,sum(DATAfinal.Oracle)/Runs-sd.HPY,'--b');
legend('HPY','Uniform','Oracle','GT','Location','NorthWest');
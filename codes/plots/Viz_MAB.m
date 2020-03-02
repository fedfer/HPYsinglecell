% PLOT OF PAPER
% be careful to use the new version of shadedErrorBar with:
% comment off edge in function
% add linewidth parameter

shadedErrorBar(1:addsample,DATAfinal_HPY, {@mean,@std}, 'lineprops', 'r--','transparent',true,'patchSaturation',0.075);
hold on
shadedErrorBar(1:addsample,DATAfinal_uniform, {@mean,@std}, 'lineprops', 'g-.','transparent',true,'patchSaturation',0.075);
shadedErrorBar(1:addsample,DATAfinal_Oracle, {@mean,@std}, 'lineprops', 'b','transparent',true,'patchSaturation',0.075);
% plot(1:addsample,sum(DATAfinal_GoodTulming)/Runs,'k');
shadedErrorBar(1:addsample,DATAfinal_GoodTulming, {@mean,@std}, 'lineprops', 'k:','transparent',true,'patchSaturation',0.075);

legend('HPY','Uniform','Oracle','GT','Location','NorthWest');
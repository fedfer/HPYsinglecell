%plot with sd

shadedErrorBar(1:addsample,DATAfinal.HPY, {@mean,@std}, 'lineprops', 'r--','transparent',true,'patchSaturation',0.075);
title('Delayed Abundance case');
hold on
% shadedErrorBar(1:addsample,DATAfinal_HPY_hyper, {@mean,@std}, 'lineprops', 'r:','transparent',true,'patchSaturation',0.075);
shadedErrorBar(1:addsample,DATAfinal.uniform, {@mean,@std}, 'lineprops', 'g-.','transparent',true,'patchSaturation',0.075);
shadedErrorBar(1:addsample,DATAfinal.Oracle, {@mean,@std}, 'lineprops', 'b','transparent',true,'patchSaturation',0.075);
shadedErrorBar(1:addsample,DATAfinal.GoodTulming, {@mean,@std}, 'lineprops', 'k:','transparent',true,'patchSaturation',0.075);

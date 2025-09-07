%%
addpath('~/CoSMoMVPA/mvpa/')
addpath('~/Repository/CommonFunctions/matplotlib/')

load('../derivatives/results/rsa_stats.mat','stats')

%%
f=figure(4);clf
f.Position(3:4) = [750 900];

%time time correlation young
res=stats.time_time_young;
subplot(3,2,4)  
RP = res.R.*res.thresholded_cluster_map;
imagesc(stats.infant_timevec,stats.adult_timevec,RP',[0 .02])
c=colorbar();
colormap viridis
a=gca;
a.YDir = 'normal';
title('RSA correlation infant-adult (<= 6 months)')
xlabel('time infant (ms)')
ylabel('time adult (ms)')
xlim(minmax(stats.infant_timevec))
ylim(minmax(stats.infant_timevec))
colormap viridis
hold on
axis square
plot(a.XLim,a.YLim,'r')

%time time correlation old
res=stats.time_time_old;
subplot(3,2,6)
RP = res.R.*res.thresholded_cluster_map;
imagesc(stats.infant_timevec,stats.adult_timevec,RP',[0 .02])
c=colorbar();
colormap viridis
a=gca;
a.YDir = 'normal';
title('RSA correlation infant-adult (> 6 months)')
xlabel('time infant (ms)')
ylabel('time adult (ms)')
xlim(minmax(stats.infant_timevec))
ylim(minmax(stats.infant_timevec))
colormap viridis
hold on
axis square
plot(a.XLim,a.YLim,'r')

%adult model correlation
co=tab10();h=[];
mt = {stats.adult_model_animacy,stats.adult_model_category,stats.adult_model_object};
for p=1:3
    t = 4-p;
    a=subplot(3,2,2);hold on
    mu = mt{t}.R;
    h(p) = plot(stats.adult_timevec,mu,'LineWidth',2,'Color',co(t,:));
    plot(stats.adult_timevec,0*stats.adult_timevec,'k--')
    idx = mt{t}.thresholded_cluster_map>0;
    plot(stats.adult_timevec(idx),0*stats.adult_timevec(idx)-t*.008-.02,'.','Color',co(t,:),'MarkerSize',20)
    xlabel('time adult (ms)')
    ylabel('correlation')
    title('RSA model fit (adults)')
    xlim(minmax(stats.adult_timevec))
end
legend(h,fliplr(stats.modelnames))

% infant model correlation (both)
co=tab20();
mtold = {stats.infant_old_model_animacy,stats.infant_old_model_category,stats.infant_old_model_object};
mtyoung = {stats.infant_young_model_animacy,stats.infant_young_model_category,stats.infant_young_model_object};
for t=1:3
    pp=[5 3 1];
    a=subplot(3,2,pp(t));hold on
    plot(stats.infant_timevec,mtold{t}.R,'LineWidth',2,'Color',co(t*2-1,:))
    plot(stats.infant_timevec,mtyoung{t}.R,'LineWidth',2,'Color',co(t*2,:))
    plot(stats.infant_timevec,0*stats.infant_timevec,'k--')

    idx = mtold{t}.thresholded_cluster_map>0;
    plot(stats.infant_timevec(idx),0*stats.infant_timevec(idx)-.015,'.','Color',co(t*2-1,:),'MarkerSize',20)
    
    idx = mtyoung{t}.thresholded_cluster_map>0;
    plot(stats.infant_timevec(idx),0*stats.infant_timevec(idx)-.015,'.','Color',co(t*2,:),'MarkerSize',20)
    
    legend({'> 6 months','<= 6 months'},'Location','NE')
    xlabel('time infant (ms)')
    ylabel('correlation')
    title(sprintf('RSA model fit: %s (infants)',stats.modelnames{t}))
    xlim(minmax(stats.infant_timevec))
    ylim([-.02 .04])
end

%%
fn = './figures/figure3';
tn = tempname;
print(gcf,'-dpng','-r500',tn)
im=imread([tn '.png']);
[i,j]=find(mean(im,3)<255);margin=2;
imwrite(im(min(i-margin):max(i+margin),min(j-margin):max(j+margin),:),[fn '.png'],'png');

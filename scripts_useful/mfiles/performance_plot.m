clear all
close all
num_orig = xlsread("performance_results.xlsx")
num=num_orig/1000;
ylabel_='GStencils/s';
ylims={[0 20],[0 20],[0 20]}

% amd 512,1024,2048
data_512=num(2:4,2:3);
data_1024=num(6:8,2:3);
data_2048=num(10:12,2:3);
plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,'./perf_results/amd_2000.png')

% cascadelake_2000 512,1024,2048
data_512=num(2:4,5:6);
data_1024=num(6:8,5:6);
data_2048=num(10:12,5:6);
plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,'./perf_results/cascadelake_2000.png')

% skylake_2000 512,1024,2048
data_512=num(2:4,8:9);
data_1024=num(6:8,8:9);
data_2048=num(10:12,8:9);
plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,'./perf_results/skylake_2000.png')

% haswell_2000 512,1024,2048
data_512=num(2:4,11:12);
data_1024=num(6:8,11:12);
data_2048=num(10:12,11:12);
plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,'./perf_results/haswell_2000.png')

% amd_10000 512,1024,2048
data_512=num(2:4, 14:15);
data_1024=num(6:8,14:15);
data_2048=num(10:12,14:15);
plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,'./perf_results/amd_10000.png')

ylabel_='Time (s)'
ylims={[0 250],[0 1000],[0 4000]}
data_512=num_orig(2:4,   17:18);
data_1024=num_orig(6:8,  17:18);
data_2048=num_orig(10:12,17:18);
plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,'./perf_results/amd_10000_sec.png')
s=1


function plot_3bars(data_512,data_1024,data_2048,ylabel_,ylims,filename) 
    fig = figure('units','normalized','position',[0,0,1,0.45]);
    subplot(1,3,1)
    b=bar(data_512)
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData,3));
    h=text(xtips1+0.1,ytips1/10*3,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold','Color','k')
    set(h,'Rotation',90,'Fontsize',18);
    xtips1 = b(2).XEndPoints;
    ytips1 = b(2).YEndPoints;
    labels1 = string(round(b(2).YData,3));
    h=text(xtips1+0.1,ytips1/10*3,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold','Color','k')
    set(h,'Rotation',90,'Fontsize',18);

    title('Grid size: 512*512*512')
    set(gca,'FontSize',20);grid on
    ylim(ylims{1})
    ylabel(ylabel_)
    xticklabels({'Stencil','Stencil+ABCs','Modeling'});
    xtickangle(0)
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',17,'FontWeight','bold')
    legend('SB','MWD')
    
    subplot(1,3,2)
    b=bar(data_1024)
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData,3));
    %% ytips1+1.6
    h=text(xtips1+0.1,ytips1/10*3,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold','Color','k')
    set(h,'Rotation',90,'Fontsize',18);
    xtips1 = b(2).XEndPoints;
    ytips1 = b(2).YEndPoints;
    labels1 = string(round(b(2).YData,3));
    h=text(xtips1+0.1,ytips1/10*3,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold','Color','k')
    set(h,'Rotation',90,'Fontsize',18);
    title('Grid size: 1024*1024*512')
    set(gca,'FontSize',20);grid on
    ylim(ylims{2})
    ylabel(ylabel_)
    xticklabels({'Stencil','Stencil+ABCs','Modeling'});
    xtickangle(0)
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',17,'FontWeight','bold')
    legend('SB','MWD')
    
    subplot(1,3,3)
    b=bar(data_2048)
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData,3));
    h=text(xtips1+0.1,ytips1/10*3,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold','Color','k')
    set(h,'Rotation',90,'Fontsize',18);
    xtips1 = b(2).XEndPoints;
    ytips1 = b(2).YEndPoints;
    labels1 = string(round(b(2).YData,3));
    h=text(xtips1+0.1,ytips1/10*3,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','FontWeight','bold','Color','k')
    set(h,'Rotation',90,'Fontsize',18);
    title('Grid size: 2048*2048*512')
    set(gca,'FontSize',20);grid on
    ylim(ylims{3})
    ylabel(ylabel_)
    xticklabels({'Stencil','Stencil+ABCs','Modeling'});
    xtickangle(0)
    a = get(gca,'XTickLabel');  
    set(gca,'XTickLabel',a,'fontsize',17,'FontWeight','bold')
    legend('SB','MWD')

    print(fig,filename,'-dpng','-r400');
    ss=1
end
function [h,a] = plotComLimiar(x,y,limiar)
    lim = y>limiar;
    a=area (x(lim),y(lim),'FaceColor',[0.5,0.5,0.5],'LineStyle',':');
    hold on
    h=plot (x, y,'LineWidth',1.5);
    hold off
    grid
end
function niceBox(fig)

ax = get(fig, 'currentaxes');

set(ax,'box','off','color','none','Units','pixels')
axes('Units', 'pixels', 'Position', get(ax,'Position'),'box','on','xtick',[],'ytick',[]);
axes(ax)

end


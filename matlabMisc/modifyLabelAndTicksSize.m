function ax = modifyLabelAndTicksSize(ax, fontsizeLab, fontsizeTicks)

set(get(ax, 'xlabel'), 'fontsize', fontsizeLab);
set(get(ax, 'ylabel'), 'fontsize', fontsizeLab);
set(ax, 'fontsize', fontsizeTicks);

return

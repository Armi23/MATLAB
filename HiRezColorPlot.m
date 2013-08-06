function HiRezColorPlot(colorplot,name,increaseFactor)

% Make data bigger and take moving averages to smooth edges
colorplot = SmoothMatrix(colorplot);
colorplot = Bigger(colorplot,increaseFactor);
colorplot(:,1:5) = colorplot(:,1:5) .* increaseFactor;
colorplot = SmoothMatrix(colorplot);
MakeColorPlot(colorplot,name);
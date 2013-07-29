function HiRezColorPlot(colorplot,name,increaseFactor)

colorplot = SmoothMatrix(colorplot);
% MakeColorPlot(colorplot,strcat('Smoothed-',name));
colorplot = Bigger(colorplot,increaseFactor);
colorplot(:,1:5) = colorplot(:,1:5) .* increaseFactor;
colorplot = SmoothMatrix(colorplot);
MakeColorPlot(colorplot,strcat('Bigger Smoothed-',name));
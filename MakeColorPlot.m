function MakeColorPlot(colorplot,outFile)

[m,~] = size(colorplot);
processed = zeros(m,80);
bregma = 0;
colorplot = SmoothMatrix(colorplot);
for i = 1:m 
   blank = colorplot(i,1); % first column has space to lesion
   lesion = colorplot(i,2); % second column has size of lesion
   cblank = colorplot(i,3);% distance to space to certain lesion
   clesion = colorplot(i,4); % size of certain lesion
   brain = colorplot(i,5);

   if (colorplot(i,6)) % Grab bregma's row
       bregma = i;
   end
   
   % first draw brain area
   processed(i,1:brain) = 1;

   % if all is 0, there is nothing here
   if ~(lesion == 0)
       processed(i,blank + 1:(lesion + blank)) = 2;
       processed(i,cblank + 1:(clesion + cblank)) = 3;
   end
end

% show generated image
figure, 
imagesc(processed);
title(outFile);
xlabel('Distance from Midline (100 microns)')
ylabel('Distance from Bregma (100 microns)')

if ~(bregma == 0) 
    above10 = mod(bregma,10);
    set(gca,'YTick',above10:5:m);
    set(gca,'YTickLabel',bregma - above10: -5: -100);
end


% % For testing, print out with slide numbering
% figure, 
% imagesc(processed);
% title(outFile);
% xlabel('Distance from Midline (100 microns)')
% ylabel('Distance from Bregma (100 microns)')
% set(gca,'YTick',0:5:m);
% set(gca,'YTickLabel',0:5:100)

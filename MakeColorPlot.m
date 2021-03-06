function MakeColorPlot(colorplot,outFile,testing)

if nargin < 3
    testing = 0; %testing is false by default
end

[m,~] = size(colorplot);
processed = zeros(m,80);
if (m > 100)
    processed = zeros(m,800);
end

bregma = 0;
for i = 1:m 
   blank = floor(colorplot(i,1)); % first column has space to lesion
   lesion = floor(colorplot(i,2)); % second column has size of lesion
   cblank = floor(colorplot(i,3));% distance to space to certain lesion
   clesion = floor(colorplot(i,4)); % size of certain lesion
   brain = floor(colorplot(i,5));
   
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
xlabel('Distance from Midline (1 mm)')
ylabel('Distance from Bregma (1 mm)')


if ~(bregma == 0) && (bregma < 100) 
    above10 = mod(bregma,10);
    set(gca,'YTick',above10:5:m);
    set(gca,'YTickLabel',(bregma - above10) / 10: -0.5: -100);
    set(gca,'XTick',0:10:80);
    set(gca,'XTickLabel',0:1:8);
    
elseif (bregma > 100) 
    % Find important values to modify graph
    modBregma = bregma;
    while (modBregma > 100)
        modBregma = modBregma / 10;
    end
    above10 = mod(modBregma,10);
    
    % Start y-axis on bregma and count ticks by 500 microns each way 
    set(gca,'YTick',above10:50:m);
    nearestTick = mod(above10,5);
    set(gca,'YTickLabel',(modBregma - nearestTick)/10: -0.5: -100);
    
    % Measure x-axis in 100 microns with 1mm ticks
    set(gca,'XTick',0:100:800);
    set(gca,'XTickLabel',0:1:8);
end

% For testing, print out with normal slide numbering
if (testing)
    figure, 
    imagesc(processed);
    title(outFile);
    xlabel('Distance from Midline (1 mm)')
    ylabel('Distance from Bregma (1 mm)')
    set(gca,'YTick',0:5:m);
    set(gca,'YTickLabel',0:5:100)
end
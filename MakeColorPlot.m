function MakeColorPlot(colorplot,outFile)

[m,~] = size(colorplot);
processed = zeros(m,80);
bregma = 0;
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

%      % first draw brain area
%    processed(i,1:floor(brain)) = 1;
%     
%     % if all is 0, there is nothing here
%    if ~(lesion == 0)
%        processed(i,floor(blank + 1):floor((lesion + blank))) = 2;
%        processed(i,floor(cblank + 1):floor(clesion + cblank)) = 3;
%    end

%    % first draw brain area
%    processed(i,1:round(brain)) = 1;
%     
%    % if all is 0, there is nothing here
%    if ~(lesion == 0)
%        processed(i,round(blank + 1):round((lesion + blank))) = 2;
%        processed(i,round(cblank + 1):round(clesion + cblank)) = 3;
%    end

end

% show generated image
figure, 
imagesc(processed);
title(outFile);
xlabel('Distance from Midline (100 microns)')
ylabel('Distance from Bregma (100 microns)')
if ~(bregma == 0) 
    above10 = mod(bregma,10);
    set(gca,'YTick',above10:10:m);
    set(gca,'YTickLabel',bregma - above10: -10: -100);
end
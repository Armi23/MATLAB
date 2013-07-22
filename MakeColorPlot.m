function MakeColorPlot(colorplot,outFile)

[m,~] = size(colorplot);
processed = zeros(m,80);
for i = 1:m 
   blank = colorplot(i,1); % first column has space to lesion
   lesion = colorplot(i,2); % second column has size of lesion
   cblank = colorplot(i,3);% distance to space to certain lesion
   clesion = colorplot(i,4); % size of certain lesion
   brain = colorplot(i,5);

   % first draw brain area
   processed(i,1:brain) = 1;

   % if all is 0, there is nothing here
   if ~(lesion == 0)
       processed(i,blank + 1:(lesion + blank)) = 2;
       processed(i,cblank + 1:(clesion + cblank)) = 3;
   end

end

save(strcat('ColorPlots\',outFile),'processed','colorplot');

% show generated image
figure, imagesc(processed);
    
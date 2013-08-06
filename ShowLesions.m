function ShowLesions(plots,name)

figure('Name',strcat(name,':','1-9')); %figures have their own name

[m,~] = size(plots); 
for row = 1:m
    
    % Retrieve information
    real = plots{row,1};
    boundset = plots{row,2};
    
    % Find the place to plot
    pos = mod(row,9);
    if (pos == 0); % Matlab is not 0-indexed
        pos = 9;
    end
    
    % Select subplot and add listener to call closer view
    h = subplot(3,3,pos);
    set(h,'ButtonDownFcn',{@ShowOnePlotOnNewFigure,real,boundset,row,name});
    
    % Plot if there is information. Hold in case of multiple objects
    hold all;
    set(gca,'YDir','reverse');
    title(strcat(name,{' '},num2str(row)));
    if ~isempty(real)       
        ShowOnePlot(real,boundset)
    end
    hold off;
    
    % Move to next figure after filling one
    if (pos == 9)
        figure('Name',strcat(name,':',num2str(row+1),'-',num2str(row+10)));
    end
end

% Creates figure with one big plot only
function ShowOnePlotOnNewFigure(~,~,real,boundset,row,name)
figure;
ShowOnePlot(real,boundset)
title(strcat(name,{' '},num2str(row)));
set(gca,'YDir','reverse');

% Plots objects on current figure
function ShowOnePlot(real,boundset)

[m1,~] = size(real);
[m2,~] = size(boundset);


% Iterate over all real boundaries
for i = 1:m1
    hold all;
    r = real{i};
    
    % If there are no edited sets, which would be the case if the boundset
    % at the index is empty of if the index is larger than the size of the
    % boundsets, then we only want the real in white
    if (i > m2) || isempty(boundset{i})
        fill(r(:,1),r(:,2),'w');
        
    % Otherwise, fill the real as black and then draw the edited on top in
    % white so the lesioned area is black    
    else            
        b = boundset{i};
        fill(r(:,1),r(:,2),'k',b(:,1),b(:,2),'w');
    
    end    
end

% 3D Plot - Experimental and more for fun
% for row = 1:m
%     real = file{row,1};
%     reals{row} = real;
%     boundset = file{row,2};
%     bounds{row} = boundset;
%     xlim([0 3500]);
%     ylim([0 2500]);
%     set(gca,'YDir','reverse');
%     hold all;
%     [m1,~] = size(real);
%     [m2,~] = size(boundset);
%     for i = 1:m1
%         r = real{i};
%         [f,~] = size(r)
%         z = ones([f,1]) .* row
%         if (i > m2)
%             fill3(r(:,1),r(:,2),z,'w', 'FaceAlpha', 0.4);
%         else            
%             if ~isempty(boundset{i})
%                 b = boundset{i};
%                 [f2,~] = size(b)
%                 z2 = ones([f2,1]) .* row
%                 fill3(r(:,1),r(:,2),z,'k',b(:,1),b(:,2),z2,'w', 'FaceAlpha', 0.4);   
%             else
%                 fill3(r(:,1),r(:,2),z,'w', 'FaceAlpha', 0.4);
%             end
%         end
%     end
% end
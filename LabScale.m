function LabScale

% Find files
[filenames,pathname,~] = uigetfile('*.tif','MultiSelect','on');

% If only one file, put it in an array to adhere to code
if ~iscellstr(filenames)
    filenames = {filenames};
end

% Go through all the files and bring them into our limit
[~,n] = size(filenames);
for i = 1:n
    location = strcat(pathname,filenames{i}) %#ok<NOPRT>
    I = imread(location);
    size(I)
    [m,~] = size(I);
    ratem = 2470/m;
    newI = imresize(I,ratem);
    size(newI)
    name = strrep(filenames{i},'.tif',''); % remove .tif to add new ending
    imwrite(newI,strcat(name,'_','down','_',num2str(ratem),'.tif'),'tiff');
end
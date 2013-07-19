function LabScale


[filenames,pathname,~] = uigetfile('*.tif','MultiSelect','on');

if ~iscellstr(filenames)
    filenames = {filenames};
end

[~,n] = size(filenames);
for i = 1:n
    location = strcat(pathname,filenames{i}) %#ok<NOPRT>
    I = imread(location);
    size(I)
    [m,~] = size(I);
    ratem = 2470/m;
    newI = imresize(I,ratem);
    size(newI)
    name = strrep(filenames{i},'.tif','');
    imwrite(newI,strcat(name,'_','down','_',num2str(ratem),'.tif'),'tiff');
end
function PickFile

[FileName,PathName,~] = uigetfile('*.mat');

fileinfo = load(strcat(PathName,FileName));
name = strtok(FileName,'.');

if isfield(fileinfo,'plots')
    file = fileinfo.plots; 
    ShowLesions(file,name);
elseif isfield(fileinfo,'colorplot')
    colorplot = fileinfo.colorplot;
    MakeColorPlot(colorplot,name);
    HiRezColorPlot(colorplot,name,10);

end
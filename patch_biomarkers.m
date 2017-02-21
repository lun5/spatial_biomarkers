%% script to convert from masks to individual nuclei
function patch_biomarkers(image_dir, patch_dir, spot_name, biomarker_name, tile_size)
%% INPUT: image in RGB, binary mask
%% OUTPUT: 1st order and 2nd order features 
%% morphology features
% use bwconncomp to get the connected components
coordinate_dir = fullfile(patch_dir,num2str(tile_size), 'coordinates');
tmp = load(fullfile(coordinate_dir, [spot_name '.mat']));
xy_numcells = tmp.data;
%sum(xy_numcells(:,3))
% read the image
%{
% This is for colon cancer
split_name = strsplit(spot_name, '_');
slide_name = strjoin(split_name(1:3),'_');
pos = sprintf('%03d', str2num(split_name{end}));

if strcmp(slide_name,'AGTA_264_3')
    singlecell_dir = '20140501_083300_738';
elseif strcmp(slide_name, 'AGA_260_3')
    singlecell_dir = '20140430_015736_305';
elseif strcmp(slide_name, 'AGTA_269_3')
    singlecell_dir = '20140501_083645_537';
else
    fprintf('There is no slide named %s\n', slide_name);
end
%}
     
%AGA_260_3\Results\20140430_015736_305\SingleCellAnalysis\001';
if strcmp(biomarker_name, 'CellBoundaryPics')%isempty(biomarker_name) % for colon
    %im = imread(fullfile(image_dir, slide_name,'Results',singlecell_dir,...
    %    'SingleCellAnalysis',pos, 'NUCLEI_SEG.TIF')); % for colon
    % for breast
    im = imread(fullfile(image_dir, biomarker_name, ['NucMembCyt_' spot_name '_Colocalization.tif']));
    output_dir = fullfile(patch_dir, num2str(tile_size), 'nuclei_images');
    %im = im >=1; %colon
else
    %im = imread(fullfile(image_dir, slide_name,'AFRemoved',...
    %    [biomarker_name,'_AFRemoved_', pos, '.tif']));
    imname = fullfile(image_dir, biomarker_name, [biomarker_name '_AFRemoved_' spot_name '.tif']);
    if ~exist(imname,'file')
        fprintf('%s does not exist \n', imname);
        return
    end
    im = imread(fullfile(image_dir, biomarker_name, [biomarker_name '_AFRemoved_' spot_name '.tif']));
    output_dir = fullfile(patch_dir,num2str(tile_size),'biomarker_images', biomarker_name);
end

if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

for i = 1:size(xy_numcells, 1)
    x = xy_numcells(i,1);
    y = xy_numcells(i,2);
    output_imname = fullfile(output_dir,[spot_name '_' num2str(x) '_' num2str(y) '.tif']);
    if exist(output_imname, 'file')
        continue;
    end
    cropped_im = imcrop(im, [x+1, y+1, tile_size-1, tile_size-1]);
    imwrite(cropped_im, output_imname);
end

end
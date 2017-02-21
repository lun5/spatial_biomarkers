%% script to extract patches from the nuclei masks
% and from different biomarkers
% 2/13/17

%{
% colon
image_dir = 'N:\ColonCancerStudy';
patch_dir = 'D:\Documents\multiplex\MIL\patches';
spot_name =  'AGA_260_3_2';
tile_sizes = [64,128];
bm_names = {'','BetaCatenin', 'E_cad', 'pck26'};
imlist = dir(fullfile(patch_dir,'64','coordinates','*.mat'));
imlist = {imlist.name}';

for i = 1:length(tile_sizes)
    tile_size = tile_sizes(i);
    T1 = tic;
    parfor j = 1:length(bm_names)
       T2 = tic;
       biomarker_name = bm_names{j};
       for k = 1:length(imlist)
          spot_name =  imlist{k}(1:end-4);
          patch_biomarkers(image_dir, patch_dir, spot_name, biomarker_name, tile_size);
       end
       fprintf('Done with biomarker %s, tile size %d, in %.2f seconds\n',...
           biomarker_name, tile_size, toc(T2));
    end
    fprintf('Done with tile size %d, in %.2f seconds\n',...
           tile_size, toc(T1));
end

%}
% 2/20/17
% script to tile up breast data
image_dir = 'D:\Documents\multiplex\breast\breast-data\images';
patch_dir = 'D:\Documents\multiplex\breast\patches';
tile_sizes = [128,256];
bm_names = {'CellBoundaryPics','ER', 'PR', 'HER2', 'EGFR','panCK','S6'};
imlist = dir(fullfile(patch_dir,'128','coordinates','*.mat'));
imlist = {imlist.name}';

for i = 1:length(tile_sizes)
    tile_size = tile_sizes(i);
    T1 = tic;
    for j = 1:length(bm_names)
       T2 = tic;
       biomarker_name = bm_names{j};
       parfor k = 1:length(imlist)
          spot_name =  imlist{k}(1:end-4);
          patch_biomarkers(image_dir, patch_dir, spot_name, biomarker_name, tile_size);
       end
       fprintf('Done with biomarker %s, tile size %d, in %.2f seconds\n',...
           biomarker_name, tile_size, toc(T2));
    end
    fprintf('Done with tile size %d, in %.2f seconds\n',...
           tile_size, toc(T1));
end

%{
imdir = 'N:\ColonCancerStudy\AGA_260_3\Results\20140430_015736_305\SingleCellAnalysis\002';
imlist = dir(fullfile(imdir, '*.tif'));
imlist = {imlist.name}';
imname = 'NUCLEI_SEG.TIF';
im = imread(fullfile(im_dir, imname));
slide_name = 'AGA_260_3';
pos = 2;
seg_im = imread(fullfile(imdir, imname));
figure; imshow(imadjust(seg_im));
bin_im = seg_im > 0;
figure; imshow(bin_im); hold on;
plot(xy(:,1), xy(:,2),'bx');
tmp = load(fullfile(data_dir,'AGA_260_3_2_xy.mat'));
xy = tmp.data;
%}


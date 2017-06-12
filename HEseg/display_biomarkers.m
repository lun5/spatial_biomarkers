% Display images with 3 different biomarkers
% BetaCatenin, COX2, and MAPK

%% script to segment the data
data_dir = 'D:\Documents\multiplex\coordinates';
image_dir = 'D:\Documents\multiplex\H&E_virtual';
src_dir = 'X:\ColonCancerStudy\';
mask_dir = 'D:\Documents\multiplex\seg_output_50';
se = strel('disk',7,4);

filelist = dir(fullfile(data_dir,'*.mat'));
filelist = {filelist.name}';


%bm_names = {'Memb.CD8','Cyt.COX2','Cell.NaKATPase', 'Cell.MSH2',...
%             'Cell.pERK','Cyt.BetaCatenin','Cyt.pMAPKAPK2','Nuc.p21',...
%             'Memb.BetaCatenin','Nuc.CD8','Nuc.CD3','Cell.PI3Kp110a'};
%coordinate_dir = 'D:\Documents\multiplex\coordinates';
bm_names = {'BetaCatenin','pMAPKAPK2','COX2'};
         
single_cell_dir = {'20140430_015736_305','20140501_083300_738','20140501_083645_537'};
output_dir = 'D:\Documents\multiplex\betacat_pMAPK_COX2';
if ~exist(output_dir,'dir')
    mkdir(output_dir);
end

parfor i = 1:length(filelist)
   spot_name = filelist{i}(1:end-4);
   split_name = strsplit(spot_name,'_');
   slide_name = strjoin(split_name(1:3),'_');
   pos = str2double(split_name{4});
   
   seg_dir = fullfile(src_dir,slide_name,'Results');
   switch slide_name
       case 'AGA_260_3'
           seg_dir = fullfile(seg_dir,single_cell_dir{1},'SingleCellAnalysis');
       case 'AGTA_264_3'
           seg_dir = fullfile(seg_dir,single_cell_dir{2},'SingleCellAnalysis');
       otherwise
           seg_dir = fullfile(seg_dir,single_cell_dir{3},'SingleCellAnalysis');
   end
   imname = 'NUCLEI_SEG.TIF';
   seg_im = imread(fullfile(seg_dir,sprintf('%03d',pos), imname));  
   output_im = uint16(zeros(size(seg_im,1), size(seg_im,2), 3));
   %figure; imshow(seg_im > 0);
   %hold on; plot(x,y,'x');
   
   % load the segmentation
   tmp = load(fullfile(mask_dir,[spot_name '.mat']));
   seg = tmp.data;
   bdry = seg2bdry(seg,'imageSize');
   bdry = imdilate(bdry,se);
   
   auto_dir = fullfile(src_dir,slide_name,'AFRemoved');
   for j = 1:length(bm_names)
       afr_im_name = sprintf('%s_AFRemoved_%03d.tif',bm_names{j},pos);
       output_im(:,:,j) = imread(fullfile(auto_dir, afr_im_name));
   end
   output_im = uint8(output_im);
   r = output_im(:,:,1); r(bdry>0) = 255;
   g = output_im(:,:,2); g(bdry>0) = 255;
   b = output_im(:,:,3); b(bdry>0) = 255;
   output_im = cat(3,r,g,b);
   imwrite(output_im, fullfile(output_dir,[spot_name '.png']));
end
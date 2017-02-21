%% script to create montage of images for the breast core tissue for viewing
% Luong Nguyen 2/17/2017

%IHC_dir = 'C:\Users\luong_nguyen\Box Sync\IHC data';
IHC_dir = 'D:\Documents\multiplex\breast\IHC';
er_dir = 'ER_spots_propername';
pr_dir = 'Multiplexed Biomarker in BrCa-PR (images)';
her2_dir = 'Multiplexed Biomarker in BrCa- HER2 (images)';
ki67_dir = 'Multiplexd Biomarker inBrCa- Ki67 (images)';
egfr_dir = 'Multiplexed Biomarker in BrCa- EGFR (images)';
new_er_dir = 'Multiplexd Biomarker inBrCa- ER (images)';
he_dir = 'TMA2 H&E Images';
% pr_dir only has 46 images
% fix the ER directory to match the name with the other folders


% ta

all_bm_dirs = {new_er_dir, pr_dir, her2_dir, ki67_dir, egfr_dir, he_dir};
all_im_lists = {};

for i = 1:length(all_bm_dirs)
    %if i == 1
    %    im_list = dir(fullfile(IHC_dir, all_bm_dirs{i}, '*.tif'));
    %else
    im_list = dir(fullfile(IHC_dir, all_bm_dirs{i}, '*.jpg'));
    %end
    im_list = {im_list.name}';
    for j = 1:length(im_list)
      imname = im_list{j};
      imname_new = strrep(imname,'position ','pos');
      imname_new = strrep(imname_new,')','');
      imname_new = strrep(imname_new,'(','');
      imname_new = strrep(imname_new,' ','_');
      if ~strcmp(imname, imname_new)
        movefile(fullfile(IHC_dir,all_bm_dirs{i},imname),...
              fullfile(IHC_dir,all_bm_dirs{i},imname_new));
      end
    end
    im_list = dir(fullfile(IHC_dir, all_bm_dirs{i}, '*.jpg'));
    im_list = {im_list.name}';
    all_im_lists{i} = im_list;
end

% clean up the position for H&E
for i = 1:length(all_im_lists{6})
   imname = all_im_lists{6}{i};
   imname_split = strsplit(imname,{'pos','.jpg'});
   if length(imname_split) > 2 && isempty(str2num(imname_split{2}(end)))
      imname_new = strrep(imname,imname_split{2},...
          [imname_split{2}(end),imname_split{2}(1:end-1)]);
      movefile(fullfile(IHC_dir,all_bm_dirs{6},imname),...
              fullfile(IHC_dir,all_bm_dirs{6},imname_new));
   end
end

i = 6;
im_list = dir(fullfile(IHC_dir, all_bm_dirs{i}, '*.jpg'));
im_list = {im_list.name}';
all_im_lists{i} = im_list;
non_match_images = {};
count = 0;
for i = 1:length(all_im_lists{1})
   imname = all_im_lists{1}{i};
   if ~ismember(imname, all_im_lists{6})
       count = count + 1;
       non_match_images{count} = imname;
   end
end

montage_dir = 'montage3';
if ~exist(fullfile(IHC_dir, montage_dir),'dir')
    mkdir(fullfile(IHC_dir, montage_dir));
end
im_sizes = [1200,1600];
bm_dirs = { he_dir, new_er_dir, her2_dir};
for i = 1:length(all_im_lists{1})
   imname = all_im_lists{1}{i};
   montage_im = zeros(im_sizes(1), im_sizes(2)*3, 3);
   for j = 1:length(bm_dirs)
       bm_im = imread(fullfile(IHC_dir, bm_dirs{j}, imname));
       if size(bm_im,1) > im_sizes(1)
          bm_im = imresize(bm_im, im_sizes); 
       end
       montage_im(:,(j-1)*im_sizes(2) + 1: j*im_sizes(2),:) = bm_im;
   end
   imwrite(uint8(montage_im), fullfile(IHC_dir, montage_dir, imname));
end

disp('done');

upload_dirs = {'he','er','her2'};
url = 'http://pitt.edu/~lun5/images/';
upload_list = {};

for i = 1:length(upload_dirs)
   if ~exist(fullfile(IHC_dir,upload_dirs{i}),'dir')
       mkdir(fullfile(IHC_dir,upload_dirs{i}))
   end
   
   for j = 1:72
       imname = all_im_lists{1}{j};
       im = imread(fullfile(IHC_dir,bm_dirs{i}, imname));
       if size(im,1) == 1200
           im_rz = imresize(im,0.7);
       else
           im_rz = imresize(im,0.7/3);
       end
       
       imwrite(im_rz,fullfile(IHC_dir,upload_dirs{i},imname),'BitDepth',8);
   end
end

for i = 1:72
    imname = all_im_lists{1}{i};
    for j = 1:length(upload_dirs)
       count = count+1;
       fprintf('%s\n',[url,upload_dirs{j},'/', imname]);
    end
end


%sftp://lun5@unixs.cssd.pitt.edu/afs/pitt.edu/home/l/u/lun5/public/html/images/he/AL13-1_1B_posA1.jpg
% flip ER by 180 degree
%{
for j = 1:length(all_im_lists{1})
   imname = all_im_lists{1}{j};
   bm_im =  imread(fullfile(IHC_dir, all_bm_dirs{1}, imname));
   bm_im = imrotate(bm_im, 180);
   imwrite(bm_im,fullfile(IHC_dir, all_bm_dirs{1}, imname));
end
%}
% clean up the h&e image

%{
montage_dir = 'montage4';
if ~exist(fullfile(IHC_dir, montage_dir),'dir')
    mkdir(fullfile(IHC_dir, montage_dir));
end
im_sizes = [1200,1600];
for i = 1:length(all_im_lists{1})
   imname = all_im_lists{1}{i};
   montage_im = zeros(im_sizes(1)*2, im_sizes(2)*2, 3);
   for j = 1:length(all_bm_dirs) - 1
       bm_im = imread(fullfile(IHC_dir, all_bm_dirs{j}, imname));
       if size(bm_im,1) > im_sizes(1)
          bm_im = imresize(bm_im, im_sizes); 
       end
       rc = de2bi(j-1, 2);
       montage_im(rc(1)*im_sizes(1) + 1: (rc(1)+1)*im_sizes(1), ...
           rc(2)*im_sizes(2) + 1: (rc(2)+1)*im_sizes(2),:) = bm_im;
   end
   imwrite(uint8(montage_im), fullfile(IHC_dir, montage_dir, imname));
end

montage_small_dir = 'montage4small';
if ~exist(fullfile(IHC_dir, montage_small_dir),'dir')
    mkdir(fullfile(IHC_dir, montage_small_dir));
end

for i = 1:length(all_im_lists{1})
   imname = all_im_lists{1}{i};
   montage_im = imread(fullfile(IHC_dir, montage_dir, imname));
   montage_im = imresize(montage_im,0.3);
   imwrite(uint8(montage_im), fullfile(IHC_dir, montage_small_dir, imname));
end


disp('done');

non_match_images = {};
count = 0;
for i = 1:length(her2_matches)
   imname = all_im_lists{6}{i};
   if ~ismember(imname, all_im_lists{2})
       count = count + 1;
       non_match_images{count} = imname;
   end
end



for i = 1:length(all_im_lists{2})
   imname = all_im_lists{2}{i};
   imname_split = strsplit(imname,{' ','.jpg','pos'});
   if length(imname_split) == 4
       er_imname = [imname_split{1}, ' ', imname_split{3}, '.tif'];
   else
       er_imname = [imname_split{1}, ' ', imname_split{2}, '.tif'];
   end
   if ismember(er_imname, all_im_lists{1})
       copyfile(fullfile(IHC_dir, er_dir, er_imname), ...
           fullfile(IHC_dir, new_er_dir, imname));
   else
      fprintf('Image %s does not exist in ER dir\n', imname); 
   end
end

non_er_matches = zeros(1, length(all_im_lists{2}));
for i = 1:length(non_er_matches)
   if strcmp(all_im_lists{2}{i},all_im_lists{3}{i}) && ...
        strcmp(all_im_lists{2}{i},all_im_lists{3}{i})
        non_er_matches(i) = 1;
   end
end

%her2_matches = zeros(1, length(all_im_lists{2}));
non_match_images = {};
count = 0;
for i = 1:length(her2_matches)
   imname = all_im_lists{3}{i};
   if ~ismember(imname, all_im_lists{2})
       count = count + 1;
       non_match_images{count} = imname;
   end
end


match_im_lists = {};
for i = 1:length(all_im_lists{2})
   imname = all_im_list{2}{i};
   imname_split = strsplit(imname,{' ','.jpg','pos'});
end


patient_rep_map = containers.Map();

pr_im_list = dir(fullfile(IHC_dir, pr_dir,'*.jpg'));
pr_im_list = {im_list.name}';

% create the mapping
for i = 1:length(im_list)
    imname = im_list{i};
    imname_split = strsplit(imname,{' ','.jpg','pos'});
    if ~ismember(imname_split{1}, keys(patient_rep_map))
        patient_rep_map(imname_split{1}) = imname_split(3);
    else
        patient_rep_map(imname_split{1}) =...
            [patient_rep_map(imname_split{1}), imname_split{3}];
    end
end

% loop through some of the images and create montage
all_keys = keys(patient_rep_map);
for i = 1:length(all_keys)
    key = all_keys{i};
    replicates = patient_rep_map(key);
    
end

%first_order_stats = cellfun(@(x) firstOrder(x, CC), all_channels, 'UniformOutput', false);
%}
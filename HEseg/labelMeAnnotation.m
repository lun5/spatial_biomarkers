%% extract annotations from Reet
% May 8, 2017

%% create groundTruth files

clearvars; close all;
LabelMe_dir = 'D:\Documents\GitHub\LabelMeToolbox';
addpath(genpath(LabelMe_dir));

source_dir = 'D:\Documents\multiplex\ReetAnnotation';
HOMEIMAGES = fullfile(source_dir,'Images'); % you can set here your default folder
HOMEANNOTATIONS = fullfile(source_dir,'Annotations'); % you can set here your default folder
%HOMELMSEGMENTS = fullfile(pwd,'LMSEGMENTS','downsample2');
%HOMELMSEGMENTS = 'Z:\HEproject\data\groundTruth_512_512';
HOMELMSEGMENTS = fullfile(source_dir,'segments');

if ~exist(HOMELMSEGMENTS,'dir')
    mkdir(HOMELMSEGMENTS);
end
%folderlist = { '10feb04_static_cars_techsquare_lot'};
%LMinstall(folderlist, HOMEIMAGES, HOMEANNOTATIONS);

folderlist = {'reet_meeting_may4','meeting_reet_may18','meeting_reet_may19'};
segnames = cell(length(folderlist),1);
for i =1 :length(folderlist)
    fnames = dir(fullfile(HOMELMSEGMENTS,'users','lun5',folderlist{i},'*.mat'));
    segnames{i} = {fnames.name}';
end

segnames = cat(1,segnames{:});

%{
for jj = 1:length(folderlist)
    database = LMdatabase(fullfile(HOMEANNOTATIONS,'users','lun5',folderlist{jj}));
    [D,j] = LMquery(database, 'folder',folderlist{jj});
    Nimages = length(D);
    
    %listfiles = {'dJUtEn6DHnfd','ErvrkyRmPB','Gg9wMyaHFpc0c','IITV3iXlbhih5q','q9VDQzxnxb'};
    
    for ndx = 1:Nimages
        %     if ~isfield(D(ndx).annotation, 'object')
        %         continue;
        %     end
        imname = D(ndx).annotation.filename(1:end-4);
        fileseg = fullfile(HOMELMSEGMENTS, D(ndx).annotation.folder, [imname '.mat']);
        if exist(fileseg,'file')
            continue;
        end
        %if ismember(imname,listfiles)
        fprintf('Start with %s\n',imname);
        fprintf('is this file there? %d\n', exist(fileseg,'file'));
        T = tic;
        [img, seg, names] = LM2segments(D(ndx), [], HOMEIMAGES, HOMELMSEGMENTS);
        mult = 1;
        seg = seg(1:mult:end,1:mult:end);
        %[annotation, img] = LMread(D,ndx,HOMEIMAGES);
        %LMplot(annotation, img)
        %pause(0.1)
        close all;
        %plot the annotations %figure;LMplot(annotation, img)
        bmap = logical(seg2bdry(seg,'imageSize'));  %figure; imshow(mat2gray(bmap));
        groundTruth = cell(1,1);
        groundTruth{1,1} = struct('Segmentation',seg,'Boundaries',bmap);
        groundTruth{1,1}.names = names;
        parsave_custom(fileseg,groundTruth,'groundTruth');
        elapsed_time = toc(T);
        fprintf('\n Done with %s in %.1f seconds\n',D(ndx).annotation.filename(1:end-4), elapsed_time);
        
    end
end
%}
%% not annotated images
image_dir = 'D:\Documents\multiplex\H&E_renamed';
folder_list = {'Stage_grade_1','Stage_grade_2','Stage_grade_3'};

for i = 1:length(folder_list)
    fnames = dir(fullfile(image_dir,folder_list{i},'*.jpg'));
    fnames = {fnames.name}';
    unannotated_folder = fullfile(image_dir,folder_list{i},'unannotated');
    if ~exist(unannotated_folder,'dir')
        mkdir(unannotated_folder);
    end
    
    for j = 1:length(fnames)
       if ~ismember(lower([fnames{j}(1:end-4) '.mat']), segnames)
          copyfile(fullfile(image_dir,folder_list{i},fnames{j}), unannotated_folder);
       end
    end
end


%% From ground truth --> 3 classes we care about
% i: immune, f: fibroblas rich, til: 
seg_dir = 'D:\Documents\multiplex\ReetAnnotation\segments\users\lun5\reet_meeting_may4';
data_dir = 'D:\Documents\multiplex';
%coordinate_dir = 'D:\Documents\multiplex\coordinates_all_bm';
coordinate_dir = 'D:\Documents\multiplex\NN_data\test_biomarker_intensity\nuc_morph';
train_dir = fullfile(coordinate_dir,'train');
if ~exist(train_dir,'dir')
    mkdir(train_dir);
end

filelist = dir(fullfile(seg_dir,'*.mat'));
filelist = {filelist.name}';

labels = {'i','f','til'}; 

counts = [0,0,0];
%fibroblast = {};
%immune = {};
%til = {};

% find the index of lineage biomarkers
tmp = load(fullfile(coordinate_dir,filelist{1}));
bm_names = cellstr(tmp.bm_names);
immune_bm_names = {'Median.Nuc.CD20','Median.Nuc.CD8','Median.Nuc.CD79',...
    'Median.Nuc.CD3','Median.Nuc.EPCAM', 'Median.Nuc.SMA'};
indx = find(ismember(bm_names,immune_bm_names));
se = strel('disk',4,8);
for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    if ~exist(fullfile(coordinate_dir, [upper(imname) '.mat']),'file')
        fprintf('Do not have this file: %s, i = %d\n',imname, i);
        continue;
    end
    
    %if exist(fullfile(output_dir,[imname '.png']),'file')
    %    continue;
    %end
    
    tmp = load(fullfile(seg_dir, [imname '.mat']));
    seg = tmp.groundTruth{1}.Segmentation;
    numSegs = max(seg(:));
    names = tmp.groundTruth{1}.names;
    
    % read in the coordinates
    tmp = load(fullfile(coordinate_dir, [upper(imname) '.mat']));
    x = max(1,tmp.x);
    y = max(1,tmp.y);
    epi_stroma = tmp.epithelial == 1; % 1 is ep, 0 is stroma
    bm_data = tmp.bm_data(:,indx);
    indx_cells = sub2ind(size(seg),y,x);
    fibroblast = zeros(1,length(x));
    immune = zeros(1,length(x));
    til = zeros(1,length(x));
    for n = 1:length(names)
        name = names{n};
        if ~ismember(name,labels)
            continue;
        end
        mask = seg == n;
        dilated_mask = imdilate(mask,se);
        cells_in_mask = dilated_mask(indx_cells);
        
        if strcmp(name,'f')
            fibroblast(cells_in_mask) = 1;
        elseif strcmp(name,'i')
            immune(cells_in_mask) = 1;
        elseif strcmp(name,'til')
            til(cells_in_mask) = 1;
        else
            error('Labels have to be f, i, or til');
        end
    end
    % regrouping with the stroma/epithelial 
    fibroblast = (~epi_stroma) & logical(fibroblast);
    immune = (~epi_stroma) & logical(immune);
    til = (~epi_stroma) & logical(til);
    tmp.fibroblast = fibroblast;
    tmp.immune = immune;
    tmp.til = til;
    save(fullfile(train_dir, [upper(imname) '.mat']),'tmp');
end

%{
for i = 1:length(filelist)
    imname = filelist{i}(1:end-4);
    if ~exist(fullfile(coordinate_dir, [upper(imname) '.mat']),'file')
        fprintf('Do not have this file: %s, i = %d\n',imname, i);
        continue;
    end
    
    %if exist(fullfile(output_dir,[imname '.png']),'file')
    %    continue;
    %end
    
    tmp = load(fullfile(seg_dir, [imname '.mat']));
    seg = tmp.groundTruth{1}.Segmentation;
    numSegs = max(seg(:));
    names = tmp.groundTruth{1}.names;
    
    % read in the coordinates
    tmp = load(fullfile(coordinate_dir, [upper(imname) '.mat']));
    x = max(1,tmp.x);
    y = max(1,tmp.y);
    epi_stroma = tmp.epithelial == 1; % 1 is ep, 0 is stroma
    bm_data = tmp.bm_data(:,indx);
    indx_cells = sub2ind(size(seg),y,x);
    for n = 1:length(names)
        name = names{n};
        if ~ismember(name,labels)
            continue;
        end
        mask = seg == n;
        dilated_mask = imdilate(mask,se);
        cells_in_mask = dilated_mask(indx_cells);
        if strcmp(name,'f')
            counts(1) = counts(1) + 1;
            fibroblast{counts(1),1} = bm_data(cells_in_mask & epi_stroma,:); % epi in
            fibroblast{counts(1),2} = bm_data(cells_in_mask & (~epi_stroma),:); % stroma in
            fibroblast{counts(1),3} = bm_data(~cells_in_mask & (epi_stroma),:); % epi out
            fibroblast{counts(1),4} = bm_data(~cells_in_mask & (~epi_stroma),:); % stroma out
        elseif strcmp(name,'i')
            counts(2) = counts(2) + 1;
            immune{counts(2),1} = bm_data(cells_in_mask & epi_stroma,:); % epi in
            immune{counts(2),2} = bm_data(cells_in_mask & (~epi_stroma),:); % stroma in
            immune{counts(2),3} = bm_data(~cells_in_mask & (epi_stroma),:); % epi out
            immune{counts(2),4} = bm_data(~cells_in_mask & (~epi_stroma),:); % stroma out
        elseif strcmp(name,'til')
            counts(3) = counts(3) + 1;
            til{counts(3),1} = bm_data(cells_in_mask & epi_stroma,:); % epi in
            til{counts(3),2} = bm_data(cells_in_mask & (~epi_stroma),:); % stroma in
            til{counts(3),3} = bm_data(~cells_in_mask & (epi_stroma),:); % epi out
            til{counts(3),4} = bm_data(~cells_in_mask & (~epi_stroma),:); % stroma out

        else
            error('Labels have to be f, i, or til');
        end
        
    end
end
%}

%% plot the results
%{
bm_all_labels = {fibroblast, immune, til};
all_titles = {'epi in', 'stroma in', 'epi out','stroma out'};
all_suptitles = {'fibroblast','immune','tils'};
for i = 1:length(bm_all_labels)
   curr_bm = bm_all_labels{i};
   num_cells = 0;
   
   bm_data = cell(1,4);
   
   for j = 1:4
       bm_data{j} = cat(1,curr_bm{:,j});
       num_cells = num_cells + size(bm_data{j},1);    
   end
   
   figure('Position', [100, 100, 1000, 1000]);
   %[counts1, binCenters1] = histcounts(bm_data{1}(:,1),16);
   
   % 1: epi in; 2: stroma in; 3: epi out; 4: stroma out
   for j = 1:4
       [counts1, ~] = histcounts(bm_data{j}(:,1),0:1:16);
       [counts2, ~] = histcounts(bm_data{j}(:,2), 0:1:16);
       [counts3, ~] = histcounts(bm_data{j}(:,3), 0:1:16);
       [counts4, ~] = histcounts(bm_data{j}(:,4), 0:1:16);
       [counts5, ~] = histcounts(bm_data{j}(:,5), 0:1:16);
       [counts6, ~] = histcounts(bm_data{j}(:,6), 0:1:16);
      
       subplot(2,2,j);
%        plot(binCenters1(2:end), counts1/num_cells, 'r-','LineWidth',2);
%        hold on;
%        plot(binCenters1(2:end), counts2/num_cells, 'b-','LineWidth',2);
%        plot(binCenters1(2:end), counts3/num_cells, 'r-v','LineWidth',2);
%        plot(binCenters1(2:end), counts4/num_cells, 'b-o','LineWidth',2);
%        plot(binCenters1(2:end), counts5/num_cells, 'g-x','LineWidth',2);
%        plot(binCenters1(2:end), counts6/num_cells, 'y-v','LineWidth',2);
       
       plot(binCenters1(2:end), counts1/sum(counts1), 'r-','LineWidth',2);
       hold on;
       plot(binCenters1(2:end), counts2/sum(counts1), 'b-','LineWidth',2);
       plot(binCenters1(2:end), counts3/sum(counts1), 'r-v','LineWidth',2);
       plot(binCenters1(2:end), counts4/sum(counts1), 'b-o','LineWidth',2);
       plot(binCenters1(2:end), counts5/sum(counts1), 'g-x','LineWidth',2);
       plot(binCenters1(2:end), counts6/sum(counts1), 'y-v','LineWidth',2);

       
       legend({'CD20','CD8','CD79','CD3','EpCAM','SMA'})
       hold off;
       ylim([0 1]); ylabel('Probability')
       xlim([0 10]); xlabel('Biomarker Exp. (log scale)')
       set(gca,'FontSize',16);axis square;
       title(all_titles{j});
       
   end
   suptitle(all_suptitles{i});

end
%}


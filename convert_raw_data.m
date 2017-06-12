% process the raw data from Dan
data_dir = 'C:\Users\luong_nguyen\Box Sync\Colon Cancer Projects\data';
%tmp = load(fullfile(data_dir,'slide_struct_total56_recur5.mat'));
im_dir = 'D:\Documents\multiplex\H&E_renamed';
raw_dir = fullfile(data_dir,'raw_data_Dan');
if ~exist(raw_dir,'dir')
    mkdir(raw_dir)
end
slide = tmp.slide;

output_raw_dir = fullfile(data_dir,'raw_data_filtered');
if ~exist(output_raw_dir,'dir')
    mkdir(output_raw_dir)
end

output_csv_dir = fullfile(data_dir,'csv_data_filtered');
if ~exist(output_csv_dir,'dir')
    mkdir(output_csv_dir)
end

output_image_dir = fullfile(data_dir,'image_data_filtered');
if ~exist(output_image_dir,'dir')
    mkdir(output_image_dir)
end
%% convert to mat files
%{
for i = 1:length(tmp.slide)
    slide_name = slide(i).name;
    all_spots = slide(i).spot;
    bm_names = slide(i).biomarkers;
    t = tic;
    for j = 1:length(all_spots)
        spot_name = [slide_name '_' num2str(str2double(all_spots(j).name))];
        if exist(fullfile(output_dir,[spot_name '.mat']),'file')
            continue;
        end
        if isempty(all_spots(j).epi) 
            fprintf('Spot %s does not have epithelial cells\n',spot_name);
            continue
        end
        bm_data = cat(2,all_spots(j).epi.feats,all_spots(j).stroma.feats)';
        x = cat(2,all_spots(j).epi.xy(1,:),all_spots(j).stroma.xy(1,:));
        y = cat(2,all_spots(j).epi.xy(2,:),all_spots(j).stroma.xy(2,:));
        area = cat(2,all_spots(j).epi.area,all_spots(j).stroma.area);
        epithelial = cat(2, ones(1,length(all_spots(j).epi.area)),...
            2*ones(1,length(all_spots(j).stroma.area)));
        save(fullfile(output_dir,[spot_name '.mat']),'area','bm_data',...
            'bm_names','epithelial','x','y');
        %fprintf('Done with image %s\n',spot_name);
    end
    fprintf('Done with slide %s in %.2f seconds\n',slide_name,toc(t));
end
%}

%% compare csv and raw data
csv_dir = fullfile(data_dir,'coordinates_all_bm');
T = readtable(fullfile(data_dir,'clinical_data_all_spots.csv'),'Delimiter',',');

filelist = T.spot_name;
% match the coordinates
%{
for i = 1:length(filelist)
   fname = filelist{i};
   if ~exist(fullfile(raw_dir,[fname '.mat']),'file')
       fprintf('File %s is not in the raw data\n',fname);
       continue;
   elseif ~exist(fullfile(csv_dir,[fname '.mat']),'file')
       fprintf('File %s is not in the csv data\n',fname);
       continue;
   elseif ~exist(fullfile(im_dir,[fname '.jpg']),'file')
       fprintf('Image %s does not exist\n',fname);
       continue;
   end
   
   im = imread(fullfile(im_dir,[fname '.jpg']));
   tmp1 = load(fullfile(raw_dir,[fname '.mat']));   
   tmp2 = load(fullfile(csv_dir,[fname '.mat']));
  
   % display the position on the image
   figure; imshow(im); hold on;
   plot(tmp1.x,tmp1.y,'o');
   plot(tmp2.x,tmp2.y,'x'); hold off;
   print(fullfile(output_image_dir,[fname '_before']) ,'-dpng');

   xy = {[tmp1.x; tmp1.y]', [tmp2.x; tmp2.y]'};
   area = {tmp1.area, tmp2.area};
   indx = cell(2,1);
   fprintf('Image %s: raw has %d more cells than csv\n',fname,...
       length(area{1}) - length(area{2}));
   % matching between the two sets
   match_id = zeros(size(area{1}));
   distances = pdist2(xy{1},xy{2});
   [min_val, min_indx] = min(distances,[],2);
   min_indx(min_val > 1e-6) = 0;
   indx{1} = find(min_indx > 0);
   indx{2} = min_indx(min_indx > 0);
   fprintf('Image %s: csv has %d more cells than after filtered\n',fname,...
       length(area{2}) - length(indx{1}));
   
   % only keep the matched cells between raw and csv
   figure; imshow(im); hold on;
   plot(tmp1.x(indx{1}),tmp1.y(indx{1}),'o');
   plot(tmp2.x(indx{2}),tmp2.y(indx{2}),'x'); hold off;
   print(fullfile(output_image_dir,[fname '_after']) ,'-dpng');
   pause(.1)
   close all;
   
   % save the new set into the raw output
   area = tmp1.area(indx{1});
   x = tmp1.x(indx{1}); y = tmp1.y(indx{1});
   bm_data = tmp1.bm_data(indx{1},:);
   bm_names = tmp1.bm_names;
   epithelial = tmp1.epithelial(indx{1});
   
   save(fullfile(output_raw_dir,[fname '.mat']),'area','bm_data',...
       'bm_names','epithelial','x','y');
   
   % save the new set into the csv output
   area = tmp2.area(indx{2});
   x = tmp2.x(indx{2}); y = tmp2.y(indx{2});
   bm_data = tmp2.bm_data(indx{2},:);
   bm_names = tmp2.bm_names;
   epithelial = tmp2.epithelial(indx{2});
   
   save(fullfile(output_csv_dir,[fname '.mat']),'area','bm_data',...
       'bm_names','epithelial','x','y');
end
%}

% examine biomarker expressions of matched cells
bm_compare_dir = fullfile(data_dir,'compare_bm_expression');
if ~exist(bm_compare_dir,'dir')
    mkdir(bm_compare_dir);
end

for i = 1:length(filelist)
   fname = filelist{i};
   if ~exist(fullfile(output_raw_dir,[fname '.mat']),'file')
       fprintf('File %s is not in the raw data\n',fname);
       continue;
   elseif ~exist(fullfile(output_csv_dir,[fname '.mat']),'file')
       fprintf('File %s is not in the csv data\n',fname);
       continue;
   elseif ~exist(fullfile(im_dir,[fname '.jpg']),'file')
       fprintf('Image %s does not exist\n',fname);
       continue;
   end
   
   im = imread(fullfile(im_dir,[fname '.jpg']));
   tmp1 = load(fullfile(output_raw_dir,[fname '.mat']));   
   tmp2 = load(fullfile(output_csv_dir,[fname '.mat']));
  
   bm_names = tmp1.bm_names;
   bm_names{12} = 'ColIV';
   long_bm_names = cellfun(@(x) ['Median.Cell.' x], bm_names,'UniformOutput',false);
   indx_bm = cellfun(@(x) find(ismember(cellstr(tmp2.bm_names)',x)),...
       long_bm_names,'UniformOutput',false);
   indx_bm = cat(1,indx_bm{:});
   
   bm_data = {tmp1.bm_data,tmp2.bm_data(:,indx_bm)};
   
   % scatter plot for all 56
   figure('position',[100 100 1200 1000]); 
   ha = tight_subplot(7,8,[.03 .03],[.08 .01],[.07 .01]);
   
   for j = 1:length(bm_names)
       axes(ha(j));
       plot(bm_data{1}(:,j), bm_data{2}(:,j),'.');
       xlim([-1 15]); ylim([-1 15]);
       text(3,10,bm_names{j});grid on;
       set(gca, 'LooseInset', get(gca,'TightInset'))
   end
   print(fullfile(bm_compare_dir,fname ) ,'-dpng');
   close all;
end
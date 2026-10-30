clc, clear all, close all
%% nx 201

% Input parameters
file_start = 2;
file_end = 2158;
path = '/blue/bala/nzgheib/COVID/room_27/par';

nx = 201;
ny = 65;
nz = 201;

%% Generate grid
edgex = linspace(0,10,nx);
edgey = linspace(0,3.2,ny);
edgez = linspace(0,10,nz);
midx = (edgex(1:end-1) + edgex(2:end)) / 2;
midy = (edgey(1:end-1) + edgey(2:end)) / 2;
midz = (edgez(1:end-1) + edgez(2:end)) / 2;
[midY,midX,midZ] = meshgrid(midy,midx,midz);

%%
flag = 0;
for n=file_start:file_end
    % Read hdf5 files
    tic
    infile = strcat(path,sprintf('/par%05d.h5', n));
    disp(['Read .h5 file ',infile])
    time  = h5read(infile,'/g1/time');
    data_id  = h5read(infile,'/g1/tag');
    data_dp  = h5read(infile,'/g1/dp');
    data_loc = h5read(infile,'/g1/loc');
    data = [double(data_id) double(data_dp) double(data_loc)];
    clear data_id data_dp data_loc
    
    % filter per diameter   
    if (flag ==0)
        flag = 1;
        diam = unique(data(:,2));
        count_fltr_all_avg = zeros(size(midX,1),...
                               size(midX,2),...
                               size(midX,3),...
                               size(diam,1)+1);
    end
    count_fltr_all = zeros(size(midX,1),...
                           size(midX,2),...
                           size(midX,3),...
                           size(diam,1)+1);
    
    for i=1:size(diam,1)
        idx_diam = (data(:,2) == diam(i,1));
        data_diam = data(idx_diam,:);
        
        % Count total
        [count edges mid loc] = histcn(data_diam(:,3:5),edgex,edgey,edgez);
      
%         % Filter count
%         count_fltr = imgaussfilt3(count,sigma);
         count_fltr_all(:,:,:,i) = count;%_fltr;
    end
    
    % Count total
    [count edges mid loc] = histcn(data(:,3:5),edgex,edgey,edgez);
      
    % Filter count
%     count_fltr = imgaussfilt3(count,sigma);
    count_fltr_all(:,:,:,end) = count;%_fltr;
    
    % Add to total count
    count_fltr_all_avg = count_fltr_all_avg + count_fltr_all;
    
    toc
    
    clear count count_fltr data count_fltr_all
end
outfile = sprintf('par_count_all_nx%03d_ny%03d_nz%03d.mat', nx, ny, nz);
save(outfile,'count_fltr_all_avg','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all
%% nx 401

% Input parameters
file_start = 2;
file_end = 2158;
path = '/blue/bala/nzgheib/COVID/room_27/par';

nx = 401;
ny = 129;
nz = 401;

%% Generate grid
edgex = linspace(0,10,nx);
edgey = linspace(0,3.2,ny);
edgez = linspace(0,10,nz);
midx = (edgex(1:end-1) + edgex(2:end)) / 2;
midy = (edgey(1:end-1) + edgey(2:end)) / 2;
midz = (edgez(1:end-1) + edgez(2:end)) / 2;
[midY,midX,midZ] = meshgrid(midy,midx,midz);

%%
flag = 0;
for n=file_start:file_end
    % Read hdf5 files
    tic
    infile = strcat(path,sprintf('/par%05d.h5', n));
    disp(['Read .h5 file ',infile])
    time  = h5read(infile,'/g1/time');
    data_id  = h5read(infile,'/g1/tag');
    data_dp  = h5read(infile,'/g1/dp');
    data_loc = h5read(infile,'/g1/loc');
    data = [double(data_id) double(data_dp) double(data_loc)];
    clear data_id data_dp data_loc
    
    % filter per diameter   
    if (flag ==0)
        flag = 1;
        diam = unique(data(:,2));
        count_fltr_all_avg = zeros(size(midX,1),...
                               size(midX,2),...
                               size(midX,3),...
                               size(diam,1)+1);
    end
    count_fltr_all = zeros(size(midX,1),...
                           size(midX,2),...
                           size(midX,3),...
                           size(diam,1)+1);
    
    for i=1:size(diam,1)
        idx_diam = (data(:,2) == diam(i,1));
        data_diam = data(idx_diam,:);
        
        % Count total
        [count edges mid loc] = histcn(data_diam(:,3:5),edgex,edgey,edgez);
      
%         % Filter count
%         count_fltr = imgaussfilt3(count,sigma);
         count_fltr_all(:,:,:,i) = count;%_fltr;
    end
    
    % Count total
    [count edges mid loc] = histcn(data(:,3:5),edgex,edgey,edgez);
      
    % Filter count
%     count_fltr = imgaussfilt3(count,sigma);
    count_fltr_all(:,:,:,end) = count;%_fltr;
    
    % Add to total count
    count_fltr_all_avg = count_fltr_all_avg + count_fltr_all;
    
    toc
    
    clear count count_fltr data count_fltr_all
end
outfile = sprintf('par_count_all_nx%03d_ny%03d_nz%03d.mat', nx, ny, nz);
save(outfile,'count_fltr_all_avg','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc, clear all, close all
%% nx 801

% Input parameters
file_start = 2;
file_end = 2158;
path = '/blue/bala/nzgheib/COVID/room_27/par';

nx = 801;
ny = 257;
nz = 801;

%% Generate grid
edgex = linspace(0,10,nx);
edgey = linspace(0,3.2,ny);
edgez = linspace(0,10,nz);
midx = (edgex(1:end-1) + edgex(2:end)) / 2;
midy = (edgey(1:end-1) + edgey(2:end)) / 2;
midz = (edgez(1:end-1) + edgez(2:end)) / 2;
[midY,midX,midZ] = meshgrid(midy,midx,midz);

%%
flag = 0;
for n=file_start:file_end
    % Read hdf5 files
    tic
    infile = strcat(path,sprintf('/par%05d.h5', n));
    disp(['Read .h5 file ',infile])
    time  = h5read(infile,'/g1/time');
    data_id  = h5read(infile,'/g1/tag');
    data_dp  = h5read(infile,'/g1/dp');
    data_loc = h5read(infile,'/g1/loc');
    data = [double(data_id) double(data_dp) double(data_loc)];
    clear data_id data_dp data_loc
    
    % filter per diameter   
    if (flag ==0)
        flag = 1;
        diam = unique(data(:,2));
        count_fltr_all_avg = zeros(size(midX,1),...
                               size(midX,2),...
                               size(midX,3),...
                               size(diam,1)+1);
    end
    count_fltr_all = zeros(size(midX,1),...
                           size(midX,2),...
                           size(midX,3),...
                           size(diam,1)+1);
    
    for i=1:size(diam,1)
        idx_diam = (data(:,2) == diam(i,1));
        data_diam = data(idx_diam,:);
        
        % Count total
        [count edges mid loc] = histcn(data_diam(:,3:5),edgex,edgey,edgez);
      
%         % Filter count
%         count_fltr = imgaussfilt3(count,sigma);
         count_fltr_all(:,:,:,i) = count;%_fltr;
    end
    
    % Count total
    [count edges mid loc] = histcn(data(:,3:5),edgex,edgey,edgez);
      
    % Filter count
%     count_fltr = imgaussfilt3(count,sigma);
    count_fltr_all(:,:,:,end) = count;%_fltr;
    
    % Add to total count
    count_fltr_all_avg = count_fltr_all_avg + count_fltr_all;
    
    toc
    
    clear count count_fltr data count_fltr_all
end
outfile = sprintf('par_count_all_nx%03d_ny%03d_nz%03d.mat', nx, ny, nz);
save(outfile,'count_fltr_all_avg','-v7.3')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Filter time average
% tic
% sigma = 0.5;
% count_fltr_all_avg1 = zeros(size(count_fltr_all_avg));
% for i=1:size(count_fltr_all_avg,4)
%     count_fltr_all_avg1(:,:,:,i) = ...
%     imgaussfilt3(count_fltr_all_avg(:,:,:,i),sigma);
% end
% toc
% 
% %% Output time average
% filename = sprintf('par_s05_acum_%05d_%05d.vtk', file_start, file_end);
% vtkwrite(filename,'structured_grid',midX,midY,midZ,...
%                   'scalars','nd1',count_fltr_all_avg1(:,:,:,1),   ...
%                   'scalars','nd2',count_fltr_all_avg1(:,:,:,2),   ...
%                   'scalars','nd3',count_fltr_all_avg1(:,:,:,3),   ...
%                   'scalars','nd4',count_fltr_all_avg1(:,:,:,4),   ...
%                   'scalars','nd5',count_fltr_all_avg1(:,:,:,5),   ...
%                   'scalars','nd6',count_fltr_all_avg1(:,:,:,6),   ...
%                   'scalars','nd7',count_fltr_all_avg1(:,:,:,7),   ...
%                   'scalars','nd8',count_fltr_all_avg1(:,:,:,8),   ...
%                   'scalars','nd9',count_fltr_all_avg1(:,:,:,9),   ...
%                   'scalars','nd10',count_fltr_all_avg1(:,:,:,10),   ...                      
%                   'scalars','ndtot',count_fltr_all_avg1(:,:,:,11),   ...
%                   'BINARY')  

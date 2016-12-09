close all
clear all
clc

files = {'balls','football','baseball','simple','football2'};
filePath = '/home/pramit/Documents/MATLAB/Modified CV/segmentationimages/';
resizing = 200;
iterationsO = 600;
iterationsD = 600;
muO = 0.2;
muD = 0.2;

resize = [256, 512];

adjust = true;
num_channels = 3;
channels = {'r','g','b'};

% for number=1:size(files,2)
for number=5:5

    fprintf('File Number - %d\n',number);
    file = files{number};
    original = imread(strcat(filePath,file,'/',file,'.jpg'));

    % I = imresize(original,resize);


    tic;
    [seg,phiO] = chenvese_org(original,'whole',iterationsO,muO,'multiphase',filePath,file,resizing);
    timeTakenOriginal = toc;
    fprintf('Original CV timeTakenOriginal %d\n',timeTakenOriginal);    

    params.sV1 = [255,0,0];
    params.sV2 = [255,255,255];
    params.sV3 = [0,255,0];
    params.sV4 = [0,0,255];

    lambdasVector = [4,1];
    lambdasMultiphase = [3,3,2,1];

    tic;
    [segm,phiD] = chenvese_distance(original,'whole',iterationsD,muD,'multiphase',params,filePath,file,lambdasMultiphase,resizing);
    timeTakenDistance  = toc;
    fprintf('Modified CV timeTakenDistance %d\n',timeTakenDistance);


end
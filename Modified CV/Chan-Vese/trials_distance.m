close all
clear all
clc

files = {'simple-1','simple-2','simple-3','simple-1-2-cells','simple-2-2-cells','joined-1','joined-2','joined-3','complex-1','complex-2'};
filePath = '../segmentation-cell-images/cropped/';

% files = {'DSLC0001_2_normalized_rot','DSLC0002_2_normalized_rot','DSLC0002_normalized_rot','DSLC0004_normalized_rot','DSLC0006_normalized_rot','DSLC0009_normalized_rot','DSLC0010_normalized_rot','DSLC0011_normalized_rot','DSLC0012_normalized_rot','DSLC0013_2_normalized_rot','DSLC0014_normalized_rot'};
% filePath = '../segmentation-cell-images/Normalized/';

% files = {'1','2','3','4','5','6','7','8','9'};
% filePath = '../segmentation-cell-images/results-compilation/';

trim = 40000;
resizing = 200;
iterationsO = 600;
iterationsD = 600;
muO = 0.2;
muD = 0.2;

resize = [256, 512];

adjust = true;
num_channels = 3;
channels = {'r','g','b'};

[nV,nVa,cV,cVa,bV1,bV1a,bV2,bV2a] = getStainVectors();

for number=1:size(files,2)

    close all;

    fprintf('File Number - %d\n',number);
    file = files{number};
    original = imread(strcat(filePath,file,'/',file,'.bmp'));
    I = imread(strcat(filePath,file,'/',file,'.bmp'));
    [height,width,dim] = size(I);
    Ihsv = rgb2hsv(I);
    Ilab = rgb2lab(I);

    if adjust
        % I = imadjust(I,stretchlim(I),[]);
        I(:,:,2) = imadjust(I(:,:,2),stretchlim(I(:,:,2)),[]);
        I(:,:,3) = imadjust(I(:,:,3),stretchlim(I(:,:,3)),[]);
    end

    % I = Ihsv;

    % I = imresize(I,resize);


    tic;
    [seg,phiO] = chenvese_org(original,'whole',iterationsO,muO,'multiphase',filePath,file,resizing);
    timeTakenOriginal = toc;
    fprintf('Original CV timeTakenOriginal %d\n',timeTakenOriginal);    

    params.sV1 = cV;
    params.sV2 = nV;
    params.sV3 = bV1;
    params.sV4 = bV2;

    lambdasVector = [4,1];
    lambdasMultiphase = [3,3,2,1];

    tic;
    [segm,phiD] = chenvese_distance(I,'whole',iterationsD,muD,'multiphase',params,filePath,file,lambdasMultiphase,resizing);
    timeTakenDistance  = toc;
    fprintf('Modified CV timeTakenDistance %d\n',timeTakenDistance);


    s = resizing./min(size(I,1),size(I,2));
    resizedOriginal = imresize(original,s);
    % fig = figure('units','normalized','outerposition',[0 0 1 1]);
    fig = figure;
    subplot(1,2,1);
    showphi(resizedOriginal, phiO, number, filePath, file);
    title('Original');
    subplot(1,2,2);
    showphi(resizedOriginal, phiD, number, filePath, file);
    title('Modified');

    pause(5);

    % saveas(fig,strcat(filePath,file,'/',file,'-comparison-r-ge-be.bmp'));

end
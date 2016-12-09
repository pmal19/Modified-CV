function [nV,nVa,cV,cVa,bV1,bV1a,bV2,bV2a] = getStainVectors()

  % files = {'simple-1','simple-2','simple-3','simple-1-2-cells','simple-2-2-cells','joined-1','joined-2','joined-3','complex-1','complex-2'};
  % filePath = '/home/pramit/Documents/MATLAB/Modified CV/segmentation-cell-images/cropped/';

  filePath = '/home/pramit/Documents/Git/Modified CV/segmentation-cell-images/final/distribution/';
  files = {'DSLC0002-crop-1','DSLC0002-crop-2','DSLC0002-crop-3'};

  adjust = true;

  num_channels = 3;
  channels = {'r','g','b'};

  nuc.r = []; nuc.g = []; nuc.b = []; nuc.t = 0;
  cyt.r = []; cyt.g = []; cyt.b = []; cyt.t = 0;
  b1.r = []; b1.g = []; b1.b = []; b1.t = 0;
  b2.r = []; b2.g = []; b2.b = []; b2.t = 0;

  for number=1:size(files,2)
    
    file = files{number};
    original = imread(strcat(filePath,file,'/',file,'.png'));
    I = imread(strcat(filePath,file,'/',file,'.png'));
    [height,width,dim] = size(I);

    Ihsv = rgb2hsv(I);
    Ilab = rgb2lab(I);
    
    if adjust
      % I = imadjust(I,stretchlim(I),[]);
      I(:,:,2) = imadjust(I(:,:,2),stretchlim(I(:,:,2)),[]);
      I(:,:,3) = imadjust(I(:,:,3),stretchlim(I(:,:,3)),[]);
    end

    % I = Ihsv;

    Iintermediate = imread(strcat(filePath,file,'/',file,'-nucleus-mask-pdf.png'));
    if adjust
      % Iintermediate = imadjust(Iintermediate,stretchlim(Iintermediate),[]);
      Iintermediate(:,:,2) = imadjust(Iintermediate(:,:,2),stretchlim(Iintermediate(:,:,2)),[]);
      Iintermediate(:,:,3) = imadjust(Iintermediate(:,:,3),stretchlim(Iintermediate(:,:,3)),[]);
    end
    nucleusMask = zeros(size(I));
    nucleusMaskp = imresize(Iintermediate,1);
    flatNucleusMask = bitand(bitand(nucleusMaskp(:,:,1),nucleusMaskp(:,:,2)),nucleusMaskp(:,:,3));
    flatNucleusMask = uint8(bitand(uint16(flatNucleusMask)+1,256) / 256);
    for i=1:num_channels
        nucleusMask(:,:,i) = flatNucleusMask;
    end
    nuc.t = nuc.t + sum(sum(flatNucleusMask));

    for k = 1:num_channels
        temp = double(nucleusMask(:,:,k));
        num = sum(sum(temp));
        c = double(I(:,:,k));
        c = c.*temp;
        d = c(:);
        nuc.(channels{k}) = cat(1,nuc.(channels{k}),d);
    end

    Iintermediate = imread(strcat(filePath,file,'/',file,'-cyto-mask-pdf.png'));
    if adjust
      % Iintermediate = imadjust(Iintermediate,stretchlim(Iintermediate),[]);
      Iintermediate(:,:,2) = imadjust(Iintermediate(:,:,2),stretchlim(Iintermediate(:,:,2)),[]);
      Iintermediate(:,:,3) = imadjust(Iintermediate(:,:,3),stretchlim(Iintermediate(:,:,3)),[]);
    end
    cytoMask = zeros(size(I));
    cytoMaskp = imresize(Iintermediate,1);
    flatCytoMask = bitand(bitand(cytoMaskp(:,:,1),cytoMaskp(:,:,2)),cytoMaskp(:,:,3));
    flatCytoMask = uint8(bitand(uint16(flatCytoMask)+1,256) / 256);
    for i=1:num_channels
        cytoMask(:,:,i) = flatCytoMask;
    end
    cyt.t = cyt.t + sum(sum(flatCytoMask));

    for k = 1:num_channels
        temp = double(cytoMask(:,:,k));
        num = sum(sum(temp));
        c = double(I(:,:,k));
        c = c.*temp;
        d = c(:);
        cyt.(channels{k}) = cat(1,cyt.(channels{k}),d);
    end

    Iintermediate = imread(strcat(filePath,file,'/',file,'-background-1-mask-pdf.png'));
    if adjust
      % Iintermediate = imadjust(Iintermediate,stretchlim(Iintermediate),[]);
      Iintermediate(:,:,2) = imadjust(Iintermediate(:,:,2),stretchlim(Iintermediate(:,:,2)),[]);
      Iintermediate(:,:,3) = imadjust(Iintermediate(:,:,3),stretchlim(Iintermediate(:,:,3)),[]);
    end
    backgroundMask1 = zeros(size(I));
    backgroundMask1p = imresize(Iintermediate,1);
    flatBackgroundMask1 = bitand(bitand(backgroundMask1p(:,:,1),backgroundMask1p(:,:,2)),backgroundMask1p(:,:,3));
    flatBackgroundMask1 = uint8(bitand(uint16(flatBackgroundMask1)+1,256) / 256);
    for i=1:num_channels
        backgroundMask1(:,:,i) = flatBackgroundMask1;
    end
    b1.t = b1.t + sum(sum(flatBackgroundMask1));

    bV1 = zeros(1,num_channels);
    bV1a = zeros(1,num_channels);
    for k = 1:num_channels
        temp = double(backgroundMask1(:,:,k));
        num = sum(sum(temp));
        c = double(I(:,:,k));
        c = c.*temp;
        d = c(:);
        b1.(channels{k}) = cat(1,b1.(channels{k}),d);
    end

    Iintermediate = imread(strcat(filePath,file,'/',file,'-background-2-mask-pdf.png'));
    if adjust
      % Iintermediate = imadjust(Iintermediate,stretchlim(Iintermediate),[]);
      Iintermediate(:,:,2) = imadjust(Iintermediate(:,:,2),stretchlim(Iintermediate(:,:,2)),[]);
      Iintermediate(:,:,3) = imadjust(Iintermediate(:,:,3),stretchlim(Iintermediate(:,:,3)),[]);
    end
    backgroundMask2 = zeros(size(I));
    backgroundMask2p = imresize(Iintermediate,1);
    flatBackgroundMask2 = bitand(bitand(backgroundMask2p(:,:,1),backgroundMask2p(:,:,2)),backgroundMask2p(:,:,3));
    flatBackgroundMask2 = uint8(bitand(uint16(flatBackgroundMask2)+1,256) / 256);
    for i=1:num_channels
        backgroundMask2(:,:,i) = flatBackgroundMask2;
    end
    b2.t = b2.t + sum(sum(flatBackgroundMask2));

    bV2 = zeros(1,num_channels);
    bV2a = zeros(1,num_channels);
    for k = 1:num_channels
        temp = double(backgroundMask2(:,:,k));
        num = sum(sum(temp));
        c = double(I(:,:,k));
        c = c.*temp;
        d = c(:);
        b2.(channels{k}) = cat(1,b2.(channels{k}),d);
    end


    nV = zeros(1,num_channels);
    nVa = zeros(1,num_channels);
    cV = zeros(1,num_channels);
    cVa = zeros(1,num_channels);
    bV1 = zeros(1,num_channels);
    bV1a = zeros(1,num_channels);
    bV2 = zeros(1,num_channels);
    bV2a = zeros(1,num_channels);

    for ch = 1:num_channels
      nV(ch) = sum(nuc.(channels{ch}))/nuc.t;
      modeN = nuc.(channels{ch});
      nVa(ch) = mode(modeN(modeN>0));
      cV(ch) = sum(cyt.(channels{ch}))/cyt.t;
      modeC = cyt.(channels{ch});
      cVa(ch) = mode(modeC(modeC>0));      
      bV1(ch) = sum(b1.(channels{ch}))/b1.t;
      modeB1 = b1.(channels{ch});
      bV1a(ch) = mode(modeB1(modeB1>0));
      bV2(ch) = sum(b2.(channels{ch}))/b2.t;
      modeB2 = b2.(channels{ch});
      bV2a(ch) = mode(modeB2(modeB2>0));
    end

  end

return
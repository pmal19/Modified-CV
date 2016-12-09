function showphi(I, phi, i, filePath, file)
% show curve evolution of phi

% Copyright (c) 2009, 
% Yue Wu @ ECE Department, Tufts University
% All Rights Reserved  

for j = 1:size(phi,3)
    phi_{j} = phi(:,:,j);
end

original_image = imread(strcat(filePath,file,'/',file,'.bmp'));
% original_image = imread(strcat(filePath,file,'/',file,'.jpg'));

resize = 200./min(size(original_image,1),size(original_image,2));
% resize = [256,512];

imshow((imresize(original_image,resize)),'initialmagnification','fit','displayrange',[0 255]);
% imshow(I,'initialmagnification','fit','displayrange',[0 255]);
  hold on;

  if size(phi,3) == 1
      contour(phi_{1}, [0 0], 'r','LineWidth',4);
      contour(phi_{1}, [0 0], 'g','LineWidth',1.3);
  else
      contour(phi_{1}, [0 0], 'r','LineWidth',4);
      contour(phi_{1}, [0 0], 'x','LineWidth',1.3);
      contour(phi_{2}, [0 0], 'g','LineWidth',4);
      contour(phi_{2}, [0 0], 'x','LineWidth',1.3);
  end
  hold off; 
  title([num2str(i) ' Iterations']); 
  drawnow;

return  
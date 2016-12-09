%==========================================================================
%
%   Active contour with Chen-Vese Method 
%   for image segementation
%
%   Implemented by Yue Wu (yue.wu@tufts.edu)
%   Tufts University
%   Feb 2009
%   http://sites.google.com/site/rexstribeofimageprocessing/
% 
%   all rights reserved 
%   Last update 02/26/2009
%--------------------------------------------------------------------------
%   Usage of varibles:
%   input: 
%       I           = any gray/double/RGB input image
%       mask        = initial mask, either customerlized or built-in
%       num_iter    = total number of iterations
%       mu          = weight of length term
%       method      = submethods pick from ('chen','vector','multiphase')
%
%   Types of built-in mask functions
%       'small'     = create a small circular mask
%       'medium'    = create a medium circular mask
%       'large'     = create a large circular mask
%       'whole'     = create a mask with holes around
%       'whole+small' = create a two layer mask with one layer small
%                       circular mask and the other layer with holes around
%                       (only work for method 'multiphase')
%   Types of methods
%       'chen'      = general CV method
%       'vector'    = CV method for vector image
%       'multiphase'= CV method for multiphase (2 phases applied here)
%
%   output: 
%       phi0        = updated level set function 
%
%--------------------------------------------------------------------------
%
% Description: This code implements the paper: "Active Contours Without
% Edges" by Chan and Vese for method 'chen', the paper:"Active Contours Without
% Edges for vector image" by Chan and Vese for method 'vector', and the paper
% "A Multiphase Level Set Framework for Image Segmentation Using the 
% Mumford and Shah Model" by Chan and Vese. 
%
%--------------------------------------------------------------------------
% Deomo: Please see HELP file for details
%==========================================================================

function [seg,phi] = chenvese_distance(I,mask,num_iter,mu,method,params,filePath,file,lambdas,resizing)

  if(~exist('mu','var')) 
    mu=0.2; 
  end
  
  if(~exist('method','var')) 
    method = 'chan'; 
  end

   s = resizing./min(size(I,1),size(I,2));
   if s<1
       I = imresize(I,s);
   end
  
  if ischar(mask)
      switch lower (mask)
          case 'small'
              mask = maskcircle2(I,'small');
          case 'medium'
              mask = maskcircle2(I,'medium');
          case 'large'
              mask = maskcircle2(I,'large');              
          case 'whole'
              mask = maskcircle2(I,'whole'); 
          case 'whole+small'
              m1 = maskcircle2(I,'whole');
              m2 = maskcircle2(I,'small');
              mask = zeros(size(I,1),size(I,2),2);
              mask(:,:,1) = m1(:,:,1);
              mask(:,:,2) = m2(:,:,2);
          otherwise
              error('unrecognized mask shape name (MASK).');
      end
  else
      if s<1
          mask = imresize(mask,s);
      end
      if size(mask,1)>size(I,1) || size(mask,2)>size(I,2)
          error('dimensions of mask unmathch those of the image.')
      end
      switch lower(method)
          case 'multiphase'
              if  (size(mask,3) == 1)  
                  error('multiphase requires two masks but only gets one.')
              end
      end

  end       

  
switch lower(method)
    case 'chan'
        if size(I,3)== 3
            P = rgb2gray((I));
            P = double(P);
        elseif size(I,3) == 2
            P = 0.5.*(double(I(:,:,1))+double(I(:,:,2)));
        else
            P = double(I);
        end
        layer = 1;
        
    case 'vector'
        s = resizing./min(size(I,1),size(I,2));
        I = imresize(I,s);
        mask = imresize(mask,s);
        layer = size(I,3);
        if layer == 1
            display('only one image component for vector image')
        end
        P = double(I);
            
    case 'multiphase'
        layer = size(I,3);
        if size(I,1)*size(I,2)>resizing^2
            s = resizing./min(size(I,1),size(I,2));
            I = imresize(I,s);
            mask = imresize(mask,s);
        end
            
        P = double(I);
    otherwise
        error('!invalid method')
end

switch lower(method)
    case {'chan','vector'}
        
        mask = mask(:,:,1);
        phi0 = bwdist(mask)-bwdist(1-mask)+im2double(mask)-.5;
        force = eps; 
        figure();
        subplot(2,2,1); imshow((I)); title('Input Image');
        subplot(2,2,2); contour(flipud(phi0), [0 0], 'r','LineWidth',1); title('initial contour');
        subplot(2,2,3); title('Segmentation');

          for n=1:num_iter
              inidx = find(phi0>=0);
              outidx = find(phi0<0);
              force_image = 0;
              for i=1:layer
                  L = im2double(P(:,:,i));
                  c1 = sum(sum(L.*Heaviside(phi0)))/(length(inidx)+eps);
                  c1_new = (c1 + lambdas(1)*params.sV1(i))/(1+lambdas(1));
                  c2 = sum(sum(L.*(1-Heaviside(phi0))))/(length(outidx)+eps);
                  % c2_new = (c2 + lambdas(2)*params.sV2(i))/(1+lambdas(2));

                  force_image=-(L-c1_new).^2+(L-c2).^2-lambdas(1)*(c1_new-params.sV1(i))^2+force_image; 
                  % force_image=-(L-c1_new).^2+(L-c2_new).^2-lambdas(1)*(c1_new-params.sV1(i))^2+lambdas(2)*(c2_new-params.sV2(i))^2 +force_image; 
              end

              force = mu*kappa(phi0)./max(max(abs(kappa(phi0))))+1/layer.*force_image;

              force = force./(max(max(abs(force))));

              dt=0.5;
              
              old = phi0;
              phi0 = phi0+dt.*force;
              new = phi0;
              indicator = checkstop(old,new,dt);

              if(mod(n,20) == 0) 
                 showphi(I,phi0,n,filePath,file);  
              end;
              if indicator
                  showphi(I,phi0,n,filePath,file);

                  seg = phi0<=0;

                  subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');

                  return;
              end
          end;
          showphi(I,phi0,n,filePath,file);

          seg = phi0<=0;
          phi = phi0;

          subplot(2,2,4); imshow(seg); title('Global Region-Based Segmentation');
    case 'multiphase'

        mask1 = mask(:,:,1);
        mask2 = mask(:,:,2);
        phi1=bwdist(mask1)-bwdist(1-mask1)+im2double(mask1)-.5;
        phi2=bwdist(mask2)-bwdist(1-mask2)+im2double(mask2)-.5;
        
        fig = figure(28);
        subplot(2,2,1); 
        if layer ~= 1
            imshow((I(:,:,1:3))); title('Input Image');
        else
            imagesc(P); axis image; colormap(gray);title('Input Image');
        end
        subplot(2,2,2);
        hold on
        contour(flipud(mask1),[0,0],'r','LineWidth',2.5); 
        contour(flipud(mask1),[0,0],'x','LineWidth',1);
        contour(flipud(mask2),[0,0],'g','LineWidth',2.5);
        contour(flipud(mask2),[0,0],'x','LineWidth',1);
        title('initial contour');
        hold off
        subplot(2,2,3); title('Segmentation');
        
        
        for n=1:num_iter
              nb1 = find(phi1<1.2 & phi1>=-1.2);
              inidx1 = find(phi1>=0);
              outidx1 = find(phi1<0);

              nb2 = find(phi2<1.2 & phi2>=-1.2);
              inidx2 = find(phi2>=0);
              outidx2 = find(phi2<0);

              cc11 = intersect(inidx1,inidx2);
              cc12 = intersect(inidx1,outidx2);
              cc21 = intersect(outidx1,inidx2);
              cc22 = intersect(outidx1,outidx2);
              
              f_image11 = 0;
              f_image12 = 0;
              f_image21 = 0;
              f_image22 = 0;
              
              for i=1:layer
                  L = im2double(P(:,:,i));
          
              if isempty(cc11)
                  c11 = eps;
              else
                  %c11 = mean(L(cc11));
                  c11 = (mean(L(cc11)) + lambdas(1)*params.sV1(i))/(1+lambdas(1));
              end
              
              if isempty(cc12)
                  c12 = eps;
              else
                  %c12 = mean(L(cc12));
                  c12 = (mean(L(cc12)) + lambdas(2)*params.sV2(i))/(1+lambdas(2));
              end
              
              if isempty(cc21)
                  c21 = eps;
              else
                  %c21 = mean(L(cc21));
                  c21 = (mean(L(cc21)) + lambdas(3)*params.sV3(i))/(1+lambdas(3));
              end
              
              if isempty(cc22)
                  c22 = eps;
              else
                  %c22 = mean(L(cc22));
                  c22 = (mean(L(cc22)) + lambdas(4)*params.sV4(i))/(1+lambdas(4));
              end
              
              
              f_image11=(((L-c11).^2 + lambdas(1)*(c11-params.sV1(i))^2).*Heaviside(phi1).*Heaviside(phi2)) +f_image11;
              f_image12=(((L-c12).^2 + lambdas(2)*(c12-params.sV2(i))^2).*Heaviside(phi1).*(1-Heaviside(phi2))) +f_image12;
              f_image21=(((L-c21).^2 + lambdas(3)*(c21-params.sV3(i))^2).*Heaviside(phi2).*(1-Heaviside(phi1))) +f_image21;
              f_image22=(((L-c22).^2 + lambdas(4)*(c22-params.sV4(i))^2).*(1-Heaviside(phi2)).*(1-Heaviside(phi1))) +f_image22;

              end
              
              curvature1 = mu*kappa(phi1);
              curvature1 = curvature1(nb1);

              fim1 = (-f_image11(nb1)+f_image21(nb1)-f_image12(nb1)+f_image22(nb1));
              fim1 = fim1./max(abs(fim1)+eps);

              curvature2 = mu*kappa(phi2);
              curvature2 = curvature2(nb2);

              fim2 = (-f_image11(nb2)+f_image12(nb2)-f_image21(nb2)+f_image22(nb2));
              fim2 = fim2./max(abs(fim2)+eps);

              force1 = curvature1+fim1;
              force2 = curvature2+fim2;

              dt = 1.5;
              
              old(:,:,1) = phi1;
              old(:,:,2) = phi2;

              phi1(nb1) = phi1(nb1)+dt.*force1;
              phi2(nb2) = phi2(nb2)+dt.*force2;
              
              new(:,:,1) = phi1;
              new(:,:,2) = phi2;
              
              indicator = checkstop(old,new,dt);

              if indicator 
                 showphi(I, new, n, filePath, file);
                 seg11 = (phi1>0 & phi2>0);
                 seg12 = (phi1>0 & phi2<0);
                 seg21 = (phi1<0 & phi2>0);
                 seg22 = (phi1<0 & phi2<0);

                 se = strel('disk',1);
                 aa1 = imerode(seg11,se);
                 aa2 = imerode(seg12,se);
                 aa3 = imerode(seg21,se);
                 aa4 = imerode(seg22,se);
                 seg = aa1+2*aa2+3*aa3+4*aa4;
                
                 subplot(2,2,4); imagesc(seg);axis image;title('Global Region-Based Segmentation');
%                  saveas(fig,strcat(filePath,file,'\',file,'-multiphase-segmented-g-b-channel-enhanced-024-0.1.bmp'));

                  return
              end

              phi1 = reinitialization(phi1, 0.6);
              phi2 = reinitialization(phi2, 0.6);

              if(mod(n,20) == 0) 
                 phi(:,:,1) = phi1;
                 phi(:,:,2) = phi2;
                 showphi(I, phi, n, filePath, file);
                 [n,c11,c12,c21,c22];
              end;
          end;
          phi(:,:,1) = phi1;
          phi(:,:,2) = phi2;
          showphi(I, phi, n, filePath, file);

        seg11 = (phi1>0 & phi2>0);
        seg12 = (phi1>0 & phi2<0);
        seg21 = (phi1<0 & phi2>0);
        seg22 = (phi1<0 & phi2<0);

        se = strel('disk',1);
        aa1 = imerode(seg11,se);
        aa2 = imerode(seg12,se);
        aa3 = imerode(seg21,se);
        aa4 = imerode(seg22,se);
        seg = aa1+2*aa2+3*aa3+4*aa4;

        subplot(2,2,4); imagesc(seg);axis image;title('Global Region-Based Segmentation');
        
%         saveas(fig,strcat(filePath,file,'\',file,'-multiphase-segmented-g-b-channel-enhanced-024-0.1.bmp'));
        
end
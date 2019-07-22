clear all; close all; clc;

doReadImage = true;

%% Make image
if(doReadImage)
    filename = 'mon1thr1.gif';
    IM = logical(imread(filename,'GIF'));
    nI = size(IM,1);
    nJ = size(IM,2);
    IMeroded = IM;
else
    nI = 100;
    nJ = 100;

    % Square
    border = 10;
    IM = false(nI,nJ);
    IM(1+border:end-border,1+border:end-border) = true;
    IMeroded = IM;
end

figure
subplot(1,3,1)
imshow(IM,'InitialMagnification','fit')

%% Erode image with a 3x3 rectangle
se = true(11,11);
seSize = size(se);

if(any(mod(seSize,2)~=1))
    fprintf(2,'Please pick a structure element with an odd number of cells per dim\n');
end
seWidth = (seSize-1)/2;
% [I,J] = meshgrid(1:seSize(2),1:seSize(1));
% Rsqr = (J-seWidth(1)-1).^2 + (I-seWidth(2)-1).^2;
% se = Rsqr <= seWidth(1)^2;

subplot(1,3,2)
imshow(se,'InitialMagnification','fit')


% loop over image pixels
for i=1+seWidth(1):nI-seWidth(1)
    for j=1+seWidth(2):nJ-seWidth(2)

        % loop over neighboring pixels inside of the structure element
        allPixelsWhite = true;
        
%         % C++ way
%         for ii=1:seSize(1)
%             for jj=1:seSize(2)
%                
%                 iPixel = i-seWidth(1)+ii;
%                 jPixel = j-seWidth(2)+jj;
%                 
%                 if(iPixel < 1 || iPixel > nI)
%                     continue;
%                 end
%                 if(jPixel < 1 || jPixel > nJ)
%                     continue;
%                 end
%                 
%                 % Skip this cell in the structure element if it doesn't
%                 % contribute
%                 if(~se(ii,jj))
%                     continue;
%                 end
%                 
%             
%                 if(~IM(iPixel,jPixel))
%                     allPixelsWhite = false;
%                 end
% 
%                 
%             end
%         end

        % Matlab way
        
        iMin = max(1,i-seWidth(1));
        iMax = min(nI,i+seWidth(1));
        jMin = max(1,j-seWidth(2));
        jMax = min(nJ,j+seWidth(2));
        allPixelsWhite = ~any(any(~IM(iMin:iMax, jMin:jMax) & se));
        
        
        % Set all pixels in the structure element to false if any
        % are found to be a background color
        if(~allPixelsWhite)
%             % C++ way
%             for ii=1:seSize(1)
%                 for jj=1:seSize(2)
% 
%                     iPixel = i-seWidth(1)+ii;
%                     jPixel = j-seWidth(2)+jj;
%                     if(iPixel < 1 || iPixel > nI)
%                         continue;
%                     end
%                     if(jPixel < 1 || jPixel > nJ)
%                         continue;
%                     end
% 
%                     IMeroded(iPixel,jPixel) = false;
%                 end
%             end
            
            % MATLAB way
            iMin = max(1,i-seWidth(1));
            iMax = min(nI,i+seWidth(1));
            jMin = max(1,j-seWidth(2));
            jMax = min(nJ,j+seWidth(2));
            IMeroded(iMin:iMax, jMin:jMax) = false;

%             subplot(1,3,3)
%             imshow(IMeroded,'InitialMagnification','fit')
%             drawnow;
        end
        
    end
end

%%

subplot(1,3,3)
imshow(IMeroded,'InitialMagnification','fit')

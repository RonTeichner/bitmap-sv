function VerilogBitmapArray(inputImageFileName,outputVerilogFileName,sProcessing)
% function VerilogBitmapArray(inputImageFileName,outputVerilogFileName,sProcessing)
%
% Function Description:
% --------------------------------
% This function receives an image and outputs a systemVerilog bitmap by the
% name outputVerilogFileName.
% It supports the image formats jpg, tif and png. Other formats by be
% supported as well, but only jpg, tif and png were tested.
% The function allows cropping, resizing and, if the image has
% transparency, performing binary transparency.
%
% Inputs:
% --------------------------------
% inputImageFileName - a string, image input file name. image should be
%   located in the working directory. for example "corn.tif".
% outputVerilogFileName - a string, output file name. for example
%   "cornSvBitmap.sv". if empty ([]) no system verilog file will be
%   written.
% sProcessing.
%   sCrop.
%       enable - boolean, if true image will be cropped to user specified
%           values: the center of the cropped image will be at
%           sProcessing.sCrop.xyCenter and the cropped picture will contain
%           sProcessing.sCrop.xyPortions of the original picture.
%       xyPortions - [portionOf_x , portionOf_y] units are [%] (i.e values from 0 to 100)
%       xyCenter - [center_x , center_y] -  units are [%] (i.e values from 0 to 100)
%   sResize.
%       enable - boolean, if true the cropped image is resized to user
%           specified values in sProcessing.sResize.new_xy.
%       new_xy - [#columns , #rows] = no. of rows and columns in the
%           resized picture
%   quantize_nBits - no. of bits in the quantized image:
%       case (quantize_nBits == 8): red: 3bits, green: 3bits, blue: 3bits.
%       case (quantize_nBits == 4): red: 2bit, green: 1bit1, blue: 1bit1
%       case (quantize_nBits == 1): black&white image
%   binaryTransparencyTh - for png images there sometimes is a transparency
%       value for every pixel. In the output verilog file we support a
%       boolean transparency. binaryTransparencyTh is a value in the range
%       [0,100] that is used as a boolean transparency threshold.
%
% Outputs:
% --------------------------------
% A system-verilog file with the name specified in outputVerilogFileName is
%   written to the working directory.
%
% Ron Teichner, 25.03.2019

%% read from txt if no inputs to function

if nargin < 3
    filename = 'VerilogBitmapArrayTxtNumericInpus.txt';
    if ~isfile(filename)
        disp('VerilogBitmapArrayTxtNumericInpus.txt was not found in working directory');
        return;
    end
    VerilogBitmapArrayTxtNumericInpus = VerilogBitmapArrayImportTxtNumericInpus(filename, 1, 7);
    
    sProcessing.binaryTransparencyTh = VerilogBitmapArrayTxtNumericInpus(1,1);
    sProcessing.sCrop.enable = VerilogBitmapArrayTxtNumericInpus(2,1);
    sProcessing.sCrop.xyPortions = VerilogBitmapArrayTxtNumericInpus(3,:);
    sProcessing.sCrop.xyCenter = VerilogBitmapArrayTxtNumericInpus(4,:);
    sProcessing.sResize.enable = VerilogBitmapArrayTxtNumericInpus(5,1);
    sProcessing.sResize.new_xy = VerilogBitmapArrayTxtNumericInpus(6,:);
    sProcessing.quantize_nBits = VerilogBitmapArrayTxtNumericInpus(7,1);
    
    filename = 'VerilogBitmapArrayTxtFileNameInpus.txt';
    if ~isfile(filename)
        disp('VerilogBitmapArrayTxtFileNameInpus.txt was not found in working directory');
        return;
    end
    VerilogBitmapArrayTxtFileNameInpus = VerilogBitmapArrayTxtFileNameInpusFunc(filename, 1, 2);
    
    inputImageFileName = VerilogBitmapArrayTxtFileNameInpus{1};
    outputVerilogFileName = VerilogBitmapArrayTxtFileNameInpus{2};
end

%% read image
if ~isstring(inputImageFileName)
    disp('inputImageFileName should be a string, for example: inputImageFileName="ee.png"')
end

if ~isfile(inputImageFileName)
    errStr = strcat(inputImageFileName,' was not found in the current working directory');
    disp(errStr);
    disp('exiting');
    return;
end

[A,map,transparency] = imread(inputImageFileName);

% convert to rgb:
if ~isempty(map)
    A = ind2rgb(A,map);
end
A = double(A);
A = A./max(A(:)); % scaling A values to capture full [0,1] range

if numel(size(A)) == 1
    disp('image has only one dimension')
    return;
end

binaryTransp = false(size(A,1) , size(A,2));
if ~isempty(transparency)
    disp('image has transparency; scale binary threshold using sProcessing.binaryTransparencyTh');
    transpA = A;
    aboveThIdx = transparency < sProcessing.binaryTransparencyTh/100*max(transparency(:));
    binaryTransp(aboveThIdx) = true;
    for colorIdx=1:size(transpA,3)
        tmp = transpA(:,:,colorIdx);
        tmp(aboveThIdx) = 1;
        transpA(:,:,colorIdx) = tmp;
    end
else
    disp('image has no transparency');
    transpA = A;
end

%% crop image
[croppedA,croppedBinaryTransp] = BitmapCropping(transpA,binaryTransp,sProcessing);

%% resize image
[croppedResizedA,croppedResizedBinaryTransp] = BitmapResize(croppedA,croppedBinaryTransp,sProcessing);

%% Quantize image
% scale values to [0,255]:

switch sProcessing.quantize_nBits
    case 8
        nBitsRed = 3; nBitsGreen = 3; nBitsBlue = 2;
        [croppedResizedQuantizedFixedPointA,croppedResizedQuantizedDoubleA] = ...
            BitmapQuant(nBitsRed,nBitsGreen,nBitsBlue,croppedResizedA);
    case 4
        nBitsRed = 2; nBitsGreen = 1; nBitsBlue = 1;
        [croppedResizedQuantizedFixedPointA,croppedResizedQuantizedDoubleA] = ...
            BitmapQuant(nBitsRed,nBitsGreen,nBitsBlue,croppedResizedA);
    case 1
        nBitsRed = 0; nBitsGreen = 0; nBitsBlue = 0;
        [Ridx,Gidx,Bidx] = deal(1,2,3);
        whiteTh = 0.5;
        greyScaleImage = 0.2989 * croppedResizedA(:,:,Ridx) + 0.5870 * croppedResizedA(:,:,Gidx) + 0.1140 * croppedResizedA(:,:,Bidx);
        croppedResizedQuantizedFixedPointA = greyScaleImage > whiteTh;
        croppedResizedQuantizedDoubleA = croppedResizedQuantizedFixedPointA;
    otherwise
        disp('Quantization level @ sProcessing.quantize_nBits not supported');
end

%% write systemVerilog file and mif file
if ~isstring(outputVerilogFileName)
    disp('outputVerilogFileName is not a string - no systemVerilog file was created; example: outputVerilogFileName = "bitmapArray.sv"');
else
    BitmapWriteSvFile(outputVerilogFileName,croppedResizedQuantizedFixedPointA,croppedResizedBinaryTransp,sProcessing.quantize_nBits,nBitsRed,nBitsGreen,nBitsBlue);
    ramSize = numel(croppedResizedQuantizedFixedPointA);
    wordSize = sProcessing.quantize_nBits;
    BitmapMifGen(outputVerilogFileName,croppedResizedQuantizedFixedPointA(:),wordSize,ramSize);
end
%% figures
[Ridx,Gidx,Bidx] = deal(1,2,3);
f1 = figure;
movegui(f1,'west');
subplot(2,3,1); im=imshow(A); title('original image');
if ~isempty(transparency)
    im.AlphaData = transparency;
end
subplot(2,3,2); imshow(transpA); title('binary transparency');
subplot(2,3,3); imshow(croppedA); title('cropped image');
subplot(2,3,4); imshow(croppedResizedA); title('resized image');
subplot(2,3,5); imshow(croppedResizedQuantizedDoubleA); title('Quantized image');

f2 = figure;
movegui(f2,'east');
rValues = croppedResizedA(:,:,Ridx); rValues = rValues(:);
gValues = croppedResizedA(:,:,Gidx); gValues = gValues(:);
bValues = croppedResizedA(:,:,Bidx); bValues = bValues(:);



subplot(2,3,1); hist(rValues,50); title('hist red');
subplot(2,3,2); hist(gValues,50); title('hist green');
subplot(2,3,3); hist(bValues,50); title('hist blue');
switch sProcessing.quantize_nBits
    case {8,4}
        rValuesQ = croppedResizedQuantizedDoubleA(:,:,Ridx); rValuesQ = rValuesQ(:);
        gValuesQ = croppedResizedQuantizedDoubleA(:,:,Gidx); gValuesQ = gValuesQ(:);
        bValuesQ = croppedResizedQuantizedDoubleA(:,:,Bidx); bValuesQ = bValuesQ(:);
        subplot(2,3,4); hist(rValuesQ,50);title('hist red quantized');
        subplot(2,3,5); hist(gValuesQ,50);title('hist green quantized');
        subplot(2,3,6); hist(bValuesQ,50);title('hist blue quantized');
end

end

function [croppedResizedQuantizedFixedPointA,croppedResizedQuantizedDoubleA] = BitmapQuant(nBitsRed,nBitsGreen,nBitsBlue,croppedResizedA)
[Ridx, Gidx, Bidx] = deal(1,2,3);
nColors = size(croppedResizedA,3);

if Ridx <= nColors
    nLevelsRed = 2^nBitsRed;
    croppedResizedQuantizedFixedPointA(:,:,Ridx) = round(croppedResizedA(:,:,Ridx) * (nLevelsRed - 1));
else
    croppedResizedQuantizedFixedPointA(:,:,Ridx) = zeros(size(croppedResizedA(:,:,1)));
end
croppedResizedQuantizedDoubleA(:,:,Ridx) = croppedResizedQuantizedFixedPointA(:,:,Ridx) ./ (nLevelsRed - 1);

if Gidx <= nColors
    nLevelsGreen = 2^nBitsGreen;
    croppedResizedQuantizedFixedPointA(:,:,Gidx) = round(croppedResizedA(:,:,Gidx) * (nLevelsGreen - 1));
else
    croppedResizedQuantizedFixedPointA(:,:,Gidx) = zeros(size(croppedResizedA(:,:,1)));
end
croppedResizedQuantizedDoubleA(:,:,Gidx) = croppedResizedQuantizedFixedPointA(:,:,Gidx) ./ (nLevelsGreen - 1);

if Bidx <= nColors
    nLevelsBlue = 2^nBitsBlue;
    croppedResizedQuantizedFixedPointA(:,:,Bidx) = round(croppedResizedA(:,:,Bidx) * (nLevelsBlue - 1));
else
    croppedResizedQuantizedFixedPointA(:,:,Bidx) = zeros(size(croppedResizedA(:,:,1)));
end
croppedResizedQuantizedDoubleA(:,:,Bidx) = croppedResizedQuantizedFixedPointA(:,:,Bidx) ./ (nLevelsBlue - 1);

end

function [croppedA,croppedbinaryTransp] = BitmapCropping(A,binaryTransp,sProcessing)
% function [croppedA,croppedbinaryTransp] = BitmapCropping(A,binaryTransp,sProcessing)
%
% Function Description:
% --------------------------------
% This function receives an image, a trasparency matrix and cropping values
% and outputs a cropped image.
%
% Inputs:
% --------------------------------
% A - an RGB image
% binaryTransp - a transparency matrix with size(A)
% sProcessing.
%   sCrop.
%       enable - boolean, if true image will be cropped to user specified
%           values: the center of the cropped image will be at
%           sProcessing.sCrop.xyCenter and the cropped picture will contain
%           sProcessing.sCrop.xyPortions of the original picture.
%       xyPortions - [portionOf_x , portionOf_y] units are [%] (i.e values from 0 to 100)
%       xyCenter - [center_x , center_y] -  units are [%] (i.e values from 0 to 100)
%
% Outputs:
% --------------------------------
% croppedA - cropped image A
% croppedbinaryTransp - cropped transparency matrix
%
% Ron Teichner, 25.03.2019

croppedA = A; croppedbinaryTransp = binaryTransp;
if isfield(sProcessing , 'sCrop')
    if isfield(sProcessing.sCrop , 'enable')
        if sProcessing.sCrop.enable
            if ~isfield(sProcessing.sCrop , 'xyPortions')
                disp('field sProcessing.sCrop.xyPortions is missing - image was not cropped');
            elseif ~isfield(sProcessing.sCrop , 'xyCenter')
                disp('field sProcessing.sCrop.xyCenter is missing - image was not cropped');
            else
                if ( (numel(sProcessing.sCrop.xyPortions) == 2) && (numel(sProcessing.sCrop.xyCenter) == 2) )
                    r = max(min(sProcessing.sCrop.xyPortions(2),100),0);
                    c = max(min(sProcessing.sCrop.xyPortions(1),100),0);
                    
                    rC = max(min(sProcessing.sCrop.xyCenter(2),100),0);
                    cC = max(min(sProcessing.sCrop.xyCenter(1),100),0);
                    
                    halfRowsKeep = floor(size(A,1) * (r/100) / 2);
                    halfColsKeep = floor(size(A,2) * (c/100) / 2);
                    rCenter = round(size(A,1) * (rC/100));
                    cCenter = round(size(A,2) * (cC/100));
                    
                    minRow = max(rCenter - halfRowsKeep , 1);
                    maxRow = min(rCenter + halfRowsKeep , size(A,1));
                    minCol = max(cCenter - halfColsKeep , 1);
                    maxCol = min(cCenter + halfColsKeep , size(A,2));
                    
                    croppedA = A(minRow:maxRow , minCol:maxCol , :);
                    croppedbinaryTransp = binaryTransp(minRow:maxRow , minCol:maxCol , :);
                else
                    disp('sProcessing.sCrop.xyPortions or sProcessing.sCrop.xyCenter are not a [2x1] vector - image was not cropped');
                end
            end
        else
            disp('sProcessing.sCrop.enable is false - image was not cropped');
        end
    else
        disp('no field enable in sProcessing.sCrop - image was not cropped');
    end
else
    disp('no field sCrop in sProcessing - image was not cropped');
end

aboveThIdx = (croppedbinaryTransp == true);
for colorIdx=1:size(croppedA,3)
    tmp = croppedA(:,:,colorIdx);
    tmp(aboveThIdx) = 1;
    croppedA(:,:,colorIdx) = tmp;
end

end

function [croppedResizedA,croppedResizedBinaryTransp] = BitmapResize(croppedA,croppedbinaryTransp,sProcessing)
% function [croppedResizedA,croppedResizedBinaryTransp] = BitmapResize(croppedA,croppedbinaryTransp,sProcessing)
%
% Function Description:
% --------------------------------
% This function receives an image, a trasparency matrix and resizing values
% and outputs a resized image.
%
% Inputs:
% --------------------------------
% croppedA - an RGB image
% croppedbinaryTransp - a transparency matrix with size(A)
% sProcessing.
%   sResize.
%       enable - boolean, if true the cropped image is resized to user
%           specified values in sProcessing.sResize.new_xy.
%       new_xy - [#columns , #rows] = no. of rows and columns in the
%           resized picture
% Outputs:
% --------------------------------
% croppedResizedA - resized image A
% croppedResizedBinaryTransp - resized transparency matrix
%
% Ron Teichner, 25.03.2019

croppedResizedA = croppedA; croppedResizedBinaryTransp = croppedbinaryTransp;
if isfield(sProcessing , 'sResize')
    if isfield(sProcessing.sResize , 'enable')
        if sProcessing.sResize.enable
            if ~isfield(sProcessing.sResize , 'new_xy')
                disp('field sProcessing.sResize.new_xy is missing - image was not resized');
            else
                if numel(sProcessing.sResize.new_xy) == 2
                    croppedResizedA = imresize(croppedA,[sProcessing.sResize.new_xy(2) sProcessing.sResize.new_xy(1)]);
                    croppedResizedA(croppedResizedA < 0) = 0; 
                    croppedResizedA(croppedResizedA > 1) = 1; 
                    croppedResizedBinaryTransp = imresize(croppedbinaryTransp,[sProcessing.sResize.new_xy(2) sProcessing.sResize.new_xy(1)]);
                else
                    disp('sProcessing.sResize.new_xy is not a [2x1] vector - image was not resized');
                end
            end
        else
            disp('sProcessing.sResize.enable is false - image was not resized');
        end
    else
        disp('no field enable in sProcessing.sResize - image was not resized');
    end
else
    disp('no field sResize in sProcessing - image was not resized');
end

aboveThIdx = (croppedResizedBinaryTransp == true);
for colorIdx=1:size(croppedResizedA,3)
    tmp = croppedResizedA(:,:,colorIdx);
    tmp(aboveThIdx) = 1;
    croppedResizedA(:,:,colorIdx) = tmp;
end

end

function BitmapWriteSvFile(outputVerilogFileName,croppedResizedQuantizedFixedPointA,...
    croppedResizedBinaryTransp,quantize_nBits,nBitsRed,nBitsGreen,nBitsBlue)
% function BitmapWriteSvFile(outputVerilogFileName,croppedResizedQuantizedFixedPointA,croppedResizedBinaryTransp,quantize_nBits,nBitsRed,nBitsGreen,nBitsBlue)
% This function writes a bitmap in system verilog to a file named outputVerilogFileName
%
% Inputs
% ---------------------------------
% outputVerilogFileName - a string
% croppedResizedQuantizedFixedPointA - RGB image
% transparency matrix with size(croppedResizedQuantizedFixedPointA)
% quantize_nBits - numberof bits in quantized image. 8,4,1 are supported
% nBitsRed,nBitsGreen,nBitsBlue - number of bits for every color
%
% Outputs
% -------------------------------
% A system verilog file is written to the working directory
%
% Ron Teichner, 25.03.2019

enableFullSvFile = false;

nRows = size(croppedResizedQuantizedFixedPointA,1);
nCols = size(croppedResizedQuantizedFixedPointA,2);
switch quantize_nBits
    case {8,4}
        totalNumBitsPerPixel = nBitsRed + nBitsGreen + nBitsBlue;
    case 1
        totalNumBitsPerPixel = 1;
end

[Ridx, Gidx, Bidx] = deal(1,2,3);

dotIdx = strfind(outputVerilogFileName,'.');
if numel(dotIdx) > 0
    moduleName = outputVerilogFileName{1}(1:dotIdx-1);
else
    moduleName = outputVerilogFileName{1};
end

fileID = fopen(outputVerilogFileName,'w');
if enableFullSvFile
    fprintf(fileID,'module %s (\n',moduleName);
    
    fprintf(fileID,'    input logic clk,\n');
    fprintf(fileID,'    input logic resetN,\n');
    fprintf(fileID,'    input logic [10:0] offsetX, // offset from top left  position \n');
    fprintf(fileID,'    input logic [10:0] offsetY,\n');
    fprintf(fileID,'    input logic InsideRectangle, //input that the pixel is within a bracket \n');
    fprintf(fileID,'    output logic drawingRequest, //output that the pixel should be dispalyed \n');
    fprintf(fileID,'    output logic [23:0] RGBout //rgb value form the bitmap \n');
    fprintf(fileID,');\n');
    fprintf(fileID,'\n');
    
    fprintf(fileID,'localparam logic [%d:0] TRANSPARENT_ENCODING = %d''hFF ;// RGB value in the bitmap representing a transparent pixel ',[(totalNumBitsPerPixel-1) ,totalNumBitsPerPixel] );
    fprintf(fileID,'\n');
end
fprintf(fileID,'localparam  int OBJECT_HEIGHT_Y = %d;\n',nRows);
fprintf(fileID,'localparam  int OBJECT_WIDTH_X = %d;\n',nCols);
fprintf(fileID,'\n');

switch quantize_nBits
    case {8,4}
        fprintf(fileID,'logic [0:OBJECT_WIDTH_X-1] [0:OBJECT_HEIGHT_Y-1] [%d-1:0] object_colors = {\n',totalNumBitsPerPixel);
    case 1
        fprintf(fileID,'logic [0:OBJECT_WIDTH_X-1] [%d*%d-1:0] object_colors = {\n',[totalNumBitsPerPixel,nCols]);
end


for r=1:nRows
    fprintf(fileID,'{');
    for c=1:nCols
        % get pixel val:
        switch quantize_nBits
            case {8,4}
                val = croppedResizedQuantizedFixedPointA(r,c,Bidx) + croppedResizedQuantizedFixedPointA(r,c,Gidx)*2^nBitsBlue + croppedResizedQuantizedFixedPointA(r,c,Ridx)*2^(nBitsBlue+nBitsGreen);
            case 1
                val = croppedResizedQuantizedFixedPointA(r,c);
        end
        if croppedResizedBinaryTransp(r,c)
            val = 2^totalNumBitsPerPixel - 1; % value for transparency
        end
        
        
        switch quantize_nBits
            case {8}
                if c < nCols
                    fprintf(fileID,'%d''h%02X, ',[totalNumBitsPerPixel,val]);
                else
                    fprintf(fileID,'%d''h%02X ',[totalNumBitsPerPixel,val]);
                end
            case 4
                if c==1
                    fprintf(fileID,'%d''h%01X',[totalNumBitsPerPixel*nCols,val]);
                else
                    fprintf(fileID,'%01X',val);
                end
            case 1
                if c==1
                    fprintf(fileID,'%d''b%01X',[totalNumBitsPerPixel*nCols,val]);
                else
                    fprintf(fileID,'%01X',val);
                end
        end
    end
    if r < nRows
        fprintf(fileID,'},\n');
    else
        fprintf(fileID,'}\n');
    end
end

fprintf(fileID,'};');
fprintf(fileID,'\n');
fprintf(fileID,'\n');

fprintf(fileID,'wire [7:0] red_sig, green_sig, blue_sig;\n');

switch quantize_nBits
    case 8
        fprintf(fileID,'assign red_sig     = {object_colors[offsetY][offsetX][7:5] , 5''d0};\n');
        fprintf(fileID,'assign green_sig   = {object_colors[offsetY][offsetX][4:2] , 5''d0};\n');
        fprintf(fileID,'assign blue_sig    = {object_colors[offsetY][offsetX][1:0] , 6''d0};\n');
    case 4
        fprintf(fileID,'assign red_sig     = {object_colors[offsetY][offsetX][3:2] , 6''d0};\n');
        fprintf(fileID,'assign green_sig   = {object_colors[offsetY][offsetX][1:1] , 7''d0};\n');
        fprintf(fileID,'assign blue_sig    = {object_colors[offsetY][offsetX][0:0] , 7''d0};\n');
        %         fprintf(fileID,'case(offsetX)\n');
        %         for c=1:nCols
        %             fprintf(fileID,'    11''d%d :  red_sig     = {object_colors[offsetY][%d:%d] , 6''d0};\n',[c,c*4+3,c*4+2]);
        %         end
        %         fprintf(fileID,'    default :  red_sig = 8''h00;\n');
        %         fprintf(fileID,'endcase\n');
        %
        %         fprintf(fileID,'case(offsetX)\n');
        %         for c=1:nCols
        %             fprintf(fileID,'    11''d%d :  green_sig     = {object_colors[offsetY][%d:%d] , 7''d0};\n',[c,c*4+1,c*4+1]);
        %         end
        %         fprintf(fileID,'    default :  green_sig = 8''h00;\n');
        %         fprintf(fileID,'endcase\n');
        %
        %         fprintf(fileID,'case(offsetX)\n');
        %         for c=1:nCols
        %             fprintf(fileID,'    11''d%d :  blue_sig     = {object_colors[offsetY][%d:%d] , 7''d0};\n',[c,c*4+0,c*4+0]);
        %         end
        %         fprintf(fileID,'    default :  blue_sig = 8''h00;\n');
        %         fprintf(fileID,'endcase\n');
        %
        %         fprintf(fileID,'assign red_sig     = {object_colors[offsetY][offsetX*%d+3:offsetX*%d+2] , 6''d0};\n',[totalNumBitsPerPixel,totalNumBitsPerPixel]);
        %         fprintf(fileID,'assign green_sig   = {object_colors[offsetY][offsetX*%d+1:offsetX*%d+1] , 7''d0};\n',[totalNumBitsPerPixel,totalNumBitsPerPixel]);
        %         fprintf(fileID,'assign blue_sig    = {object_colors[offsetY][offsetX*%d+0:offsetX*%d+0] , 7''d0};\n',[totalNumBitsPerPixel,totalNumBitsPerPixel]);
    case 1
        fprintf(fileID,'assign red_sig     = {object_colors[offsetY][offsetX] , 6''d0};\n');
        fprintf(fileID,'assign green_sig   = {object_colors[offsetY][offsetX] , 7''d0};\n');
        fprintf(fileID,'assign blue_sig    = {object_colors[offsetY][offsetX] , 7''d0};\n');
end

fprintf(fileID,'\n');
fprintf(fileID,'\n');
if enableFullSvFile
    fprintf(fileID,'always_ff@(posedge clk)\n');
    fprintf(fileID,'begin\n');
    
    fprintf(fileID,'       RGBout      <= {red_sig,green_sig,blue_sig};\n');
    fprintf(fileID,'       if (InsideRectangle == 1''b1 ) begin // inside an external bracket \n');
    
    switch quantize_nBits
        case {8,4}
            fprintf(fileID,'            if (object_colors[offsetY][offsetX] != TRANSPARENT_ENCODING)\n');
            %     case 4
            %         fprintf(fileID,'            if (object_colors[offsetY][offsetX*(%d+1):offsetX*%d] != TRANSPARENT_ENCODING)\n',[totalNumBitsPerPixel,totalNumBitsPerPixel]);
        case 1
            fprintf(fileID,'            if (object_colors[offsetY][offsetX] == 1''b1)\n');
    end
    fprintf(fileID,'                drawingRequest <= 1''b1;\n');
    fprintf(fileID,'            else\n');
    fprintf(fileID,'                drawingRequest <= 1''b0;\n');
    fprintf(fileID,'       end\n');
    fprintf(fileID,'       else\n');
    fprintf(fileID,'            drawingRequest <= 1''b0;\n');
    
    
    fprintf(fileID,'end\n');
    
    fprintf(fileID,'\n');
    
    fprintf(fileID,'endmodule');
    
end
fclose(fileID);
disp(strcat(outputVerilogFileName, ' written to disk'));
end

function BitmapMifGen(outputVerilogFileName,values,word_size,ram_size)
% ram_size [words]
% word_size [bits]
% values - column vector
%-------------------------------------------------------------------
% Open mif file for write
%

dotIdx = strfind(outputVerilogFileName,'.');
if numel(dotIdx) > 0
    mifFileName = [outputVerilogFileName{1}(1:dotIdx-1), '.mif'];
else
    mifFileName = [outputVerilogFileName{1}, '.mif'];
end

fileID = fopen(mifFileName, 'w');

fprintf( fileID, '%s%d%s\n', 'WIDTH=', word_size, ';');
fprintf( fileID, '%s%d%s\n\n', 'DEPTH=', ram_size, ';');
fprintf( fileID, '%s\n', 'ADDRESS_RADIX=HEX;');
fprintf( fileID, '%s\n\n', 'DATA_RADIX=HEX;');
fprintf( fileID, '%s\n', 'CONTENT BEGIN');

% data format
fdata = int2str(word_size/4);
faddr = int2str(log2(ram_size)/4);
format_str = strcat('    %0', faddr, 'X  :  %0', fdata, 'X;\n');

% write values to file
idx = 0;
for ii=1:ram_size
    fprintf( fileID, format_str, idx, values(idx+1));
    idx = idx+1;
end
fprintf( fileID, '%s\n', 'END;');

fclose( fileID);
disp([mifFileName,' written to disk']);
end

function VerilogBitmapArrayTxtNumericInpus = VerilogBitmapArrayImportTxtNumericInpus(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   VERILOGBITMAPARRAYTXTNUMERICINPUS = IMPORTFILE(FILENAME) Reads data
%   from text file FILENAME for the default selection.
%
%   VERILOGBITMAPARRAYTXTNUMERICINPUS = IMPORTFILE(FILENAME, STARTROW,
%   ENDROW) Reads data from rows STARTROW through ENDROW of text file
%   FILENAME.
%
% Example:
%   VerilogBitmapArrayTxtNumericInpus = importfile('VerilogBitmapArrayTxtNumericInpus.txt', 1, 7);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/03/28 09:58:41

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
VerilogBitmapArrayTxtNumericInpus = cell2mat(raw);


end

function VerilogBitmapArrayTxtFileNameInpus = VerilogBitmapArrayTxtFileNameInpusFunc(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   VERILOGBITMAPARRAYTXTFILENAMEINPUS = IMPORTFILE(FILENAME) Reads data
%   from text file FILENAME for the default selection.
%
%   VERILOGBITMAPARRAYTXTFILENAMEINPUS = IMPORTFILE(FILENAME, STARTROW,
%   ENDROW) Reads data from rows STARTROW through ENDROW of text file
%   FILENAME.
%
% Example:
%   VerilogBitmapArrayTxtFileNameInpus = importfile('VerilogBitmapArrayTxtFileNameInpus.txt', 1, 2);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/03/28 09:54:53

%% Initialize variables.
delimiter = {''};
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));


%% Split data into numeric and string columns.
rawNumericColumns = {};
rawStringColumns = string(raw(:, 1));


%% Create output variable
VerilogBitmapArrayTxtFileNameInpus = raw;
end
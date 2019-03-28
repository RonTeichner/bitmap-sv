clear; close all; clc; 

inputImageFileName      = "corn.tif";
%inputImageFileName      = "ngc6543a.jpg";
%inputImageFileName      = "EEEL LOGO TRANSPERANT.png";
%inputImageFileName      = "peppers.png";
%inputImageFileName      = "sad-468923_960_720.jpg";

outputVerilogFileName   = "smileyBitMap.sv";

sProcessing.binaryTransparencyTh = 1; % [%]

sProcessing.sCrop.enable        = true;
sProcessing.sCrop.xyPortions    = [25,25]; % [%]
sProcessing.sCrop.xyCenter      = [50,45]; % [%] % x values are left to right; y values are up to down

sProcessing.sResize.enable  = true;
sProcessing.sResize.new_xy  = [5,5];

sProcessing.quantize_nBits = 4; % {8 - 3Red, 3Green, 2Blue ; 4 - 2Red, 1Green, 1Blue ; 1 - black&white image}

VerilogBitmapArray(inputImageFileName,outputVerilogFileName,sProcessing)
%VerilogBitmapArray

% create the exe:
% mcc -m VerilogBitmapArray.m
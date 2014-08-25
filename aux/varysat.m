function colorRGB = varysat(origColor,satValues)
% varies the hue of a specified color
% origColor is n x 3
% satValues is 1 x m, specifying range of saturation
% output is is n x m x 3 
% 

colorHSV = zeros([size(origColor),length(satValues)]);
for icond = 1:size(origColor,1)
    temp = rgb2hsv(origColor(icond,:));
    for nn = 1:length(satValues)
        temp2 = temp; temp2(2) = temp(2)*satValues(nn);
        colorHSV(icond,:,nn) = temp2;
    end
end
colorHSV(colorHSV>1) = 1;
colorHSV = permute(colorHSV,[1,3,2]);
% convert back to RGB
colorRGB = hsv2rgb(colorHSV);
function colors = dfdGetColors(nrColors)

%% Function to get colors corresponding to conditions used in manuscript figures
%
% dfdGetColors(nrColors)
% 
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This function will return a matrix with either 3 RGB color codes (blue, green, red) 
% or 4 RGB colors codes (blue, green, red, and grey). Corresponding
% to conditions Full,  Left and right (Blank).

switch nrColors
    case 3
        colors = [63, 121, 204;  116,183,74; 228, 65, 69;]/255;
    case 4
        colors = [63, 121, 204; 116,183,74; 228, 65, 69; 127,127,127]/255;
        
end

return
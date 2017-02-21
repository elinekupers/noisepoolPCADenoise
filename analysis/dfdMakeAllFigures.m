%% Script to make all figures from the paper:
% 
% AUTHORS. TITLE. JOURNAL. YEAR.

% DESCRIPTION HERE


if isempty(which('ft_prepare_layout')), dfdAddFieldtripPath, end

% Make all figures
dfdMakeFigure2;
dfdMakeFigure3;

dfdMakeFigure5;
dfdMakeFigure6;

dfdMakeFigure7;
dfdMakeFigure8;

dfdMakeFigure9;
dfdMakeFigure10;
dfdMakeFigure11; 
dfdMakeFigure12;
dfdMakeFigure13;
dfdMakeFigure14;
dfdMakeFigureVaryEpochLength;
dfdMakeFigureGridNPCvsNoisepool;


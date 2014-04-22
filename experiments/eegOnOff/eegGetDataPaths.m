function [sessionDir, conditionNames, conditionNumbers] = eegGetDataPaths(sessionum, conditionNumbers)
% [sessionDir, conditions, condisionNames] = eegGetDataPaths(sessionum)

sessionDirsAll = {...
    'Farzin_20121219_1657' ...  1
    'Farzin_20130108_1518' ...  2
    'Winawer_20121129_1224' ... 3
    'Cell_20130109_1135' ...    4
    'Gomez_20130116_1101' ...   5
    'Winawer_20130123_1632' ... 6
    'Cooper_20130228_1456' ...  7
    'Winawer_20130411_1403' ... 8
    'Norcia_20130412_1340' ...  9
    'Farzin_20130505_1813' ... 10
    'Winawer_20130530_1705' ...11
    'Hermes_20130705_1306' ... 12
    'Voyles_20130713_1105' ... 13
    'Hermes_20130719_1600' ... 14
    'Hermes_20130726_1827' ... 15
    'Hermes_20130801_1316' ... 16
    'Liu_20130807_1700' ...    17 with notch filter
    'Liu_20130807_1633' ...    18 without notch filter    
    'Simulation_20130712'  ... 19   
    };

% CONDITIONS:
conditionNamesAll = { ...
    {'on' 'off'}                        ...  1
    {'on' 'off'}                        ...  2
    {'' '' '' '' 'on' 'off'}            ...  3
    {'' '' '' '' '' '' 'on' 'off'}      ...  4
    {'on' 'off'}                        ...  5
    {'on' 'off'}                        ...  6
    {'on' 'off'}                        ...  7
    {'' '' 'on' 'off'}                  ...  8
    {'' '' 'on' 'off'}                  ...  9
    {'' '' 'on' 'off'}                  ... 10
    {'' '' 'on' 'off'}                  ... 11
    {'' '' 'on' 'off'}                  ... 12
    {'' '' 'on' 'off'}                  ... 13
    {'' '' 'on' 'off' 'left' 'right'}   ... 14
    {'on' 'off' 'left' 'right'}         ... 15
    {'on' 'off' 'upper' 'lower'}        ... 16
    {'on' 'off' 'upper' 'lower'}        ... 17
    {'on' 'off' 'upper' 'lower'}        ... 18
    {'' '' 'on' 'off'}                  ... 19
    };




if ~exist('conditionNumbers', 'var') || isempty(conditionNumbers),     
    conditionNumbers(1) = find(strcmp(conditionNamesAll{sessionum}, 'on'));
    conditionNumbers(2) = find(strcmp(conditionNamesAll{sessionum}, 'off'));
end

sessionDir     = sessionDirsAll{sessionum};
conditionNames  = conditionNamesAll{sessionum}(conditionNumbers);

disp(sessionDir); disp(conditionNames); pause(1)

end
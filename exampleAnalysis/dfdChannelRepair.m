function sensorData = dfdChannelRepair(sensorDataIn, outliers, 'nearest')

% get data sizes
nEpochs = size(sensorDataIn,2);
nChannels = size(sensorDataIn,3);

% load the xyz positions of the sensor
  net=load('meg160xyz.mat');
  
  connectivityMatrix = eye
            
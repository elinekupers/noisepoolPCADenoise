function [rho,theta,XY] = msGetThetaRho(posXY,sac,p1_n)

if ~exist('p1_n','var')
    p1_n = 1;
end

nsac = size(sac,1);
rho  = zeros(nsac,1);
theta= zeros(nsac,1);
XY   = zeros(nsac,2);


for j = 1:nsac
    % onset and offset
    a = sac(j,1);
    b = sac(j,2);
    
    % maximum displacement
        
    [~,ind] = max( sqrt ( (posXY(a:b,1) - posXY(a,1)).^2 + ...
                          (posXY(a:b,2) - posXY(a,2)).^2 ));
                      
    dx = posXY(a+ind-1,1) - nanmean(posXY(a:a+p1_n-1,1));
    dy = posXY(a+ind-1,2) - nanmean(posXY(a:a+p1_n-1,2));
    %dy = dy*3;
    
    % compute the angle at maximum displacement 
    theta(j) = mod(atan2(dy,dx),2*pi);
    rho(j)  = sqrt(dx.^2+dy.^2);
    
    XY(j,:) = [dx,dy];
end

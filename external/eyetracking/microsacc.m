function [sac, radius] = microsacc(x,vel,VFAC,MINDUR)
%-------------------------------------------------------------------
%
%  FUNCTION microsacc.m
%  Detection of monocular candidates for microsaccades;
%  Please cite: Engbert, R., & Mergenthaler, K. (2006) Microsaccades 
%  are triggered by low retinal image slip. Proceedings of the National 
%  Academy of Sciences of the United States of America, 103: 7192-7197.
%
%  (Version 2.1, 03 OCT 05)
%
%-------------------------------------------------------------------
%
%  INPUT:
%
%  x(:,1:2)         position vector
%  vel(:,1:2)       velocity vector
%  VFAC             relative velocity threshold
%  MINDUR           minimal saccade duration
%
%  OUTPUT:
%
%  sac(1:num,1)   onset of saccade
%  sac(1:num,2)   end of saccade
%  sac(1:num,3)   peak velocity of saccade (vpeak)
%  sac(1:num,4)   horizontal component     (dx)
%  sac(1:num,5)   vertical component       (dy)
%  sac(1:num,6)   horizontal amplitude     (dX)
%  sac(1:num,7)   vertical amplitude       (dY)
%
%---------------------------------------------------------------------

% compute threshold
msdx = sqrt( nanmedian(vel(:,1).^2) - (nanmedian(vel(:,1)))^2 );
msdy = sqrt( nanmedian(vel(:,2).^2) - (nanmedian(vel(:,2)))^2 );
if msdx<realmin
    msdx = sqrt( nanmean(vel(:,1).^2) - (nanmean(vel(:,1)))^2 );
    if msdx<realmin
        error('msdx<realmin in microsacc.m');
    end
end
if msdy<realmin
    msdy = sqrt( nanmean(vel(:,2).^2) - (nanmean(vel(:,2)))^2 );
    if msdy<realmin
        error('msdy<realmin in microsacc.m');
    end
end
radiusx = VFAC*msdx;
radiusy = VFAC*msdy;
radius = [radiusx radiusy];

% compute test criterion: ellipse equation
test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2;
indx = find(test>1);

% determine saccades
N = length(indx); 
sac = [];
nsac = 0;
dur = 1;
a = 1;
k = 1;
while k<N
    if indx(k+1)-indx(k)==1
        dur = dur + 1;
    else
        if dur>=MINDUR
            nsac = nsac + 1;
            b = k;
            sac(nsac,:) = [indx(a) indx(b)];
        end
        a = k+1;
        dur = 1;
    end
    k = k + 1;
end

% check for minimum duration
if dur>=MINDUR
    nsac = nsac + 1;
    b = k;
    sac(nsac,:) = [indx(a) indx(b)];
end

% compute peak velocity, horizonal and vertical components
for s=1:nsac
    % onset and offset
    a = sac(s,1); 
    b = sac(s,2); 
    % saccade peak velocity (vpeak)
    vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
    sac(s,3) = vpeak;
    % saccade vector (dx,dy)
    dx = x(b,1)-x(a,1); 
    dy = x(b,2)-x(a,2); 
    sac(s,4) = dx;
    sac(s,5) = dy;
    % saccade amplitude (dX,dY)
    i = sac(s,1):sac(s,2);
    [minx, ix1] = min(x(i,1));
    [maxx, ix2] = max(x(i,1));
    [miny, iy1] = min(x(i,2));
    [maxy, iy2] = max(x(i,2));
    dX = sign(ix2-ix1)*(maxx-minx);
    dY = sign(iy2-iy1)*(maxy-miny);
    sac(s,6:7) = [dX dY];
    
    %[th,r] = cart2pol(x(i,1),x(i,2));
    %sac(s,8) = max(r)-min(r);
end


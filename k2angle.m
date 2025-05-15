function angle=k2angle(y,varargin)
    if nargin == 1
        dt=y.end-y.start;
        y=dt(2);
        x=dt(1);
    else
        x=varargin{1};
    end
    angle=atan2d(y,x);
    if angle<0
        angle=angle+360;
    end
end
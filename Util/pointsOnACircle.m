function out = pointsOnACircle(p1,p2,p3,nump)
%POINTSONACIRCLE Summary of this function goes here
%   Detailed explanation goes here
    
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    x3 = p3(1); y3 = p3(2);

    t1 = x3^2-x2^2+y3^2-y2^2;
    t2 = x1^2-x3^2+y1^2-y3^2;
    t3 = x2^2-x1^2+y2^2-y1^2;
    d = x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3;
    xc = 1/2*(t1*y1+t2*y2+t3*y3)/d;
    yc = -1/2*(t1*x1+t2*x2+t3*x3)/d;
    radius = sqrt((xc-p1(1))^2 + (yc-p1(2))^2);
    
    a1 = atan2(yc-y1,x1-xc);
    if a1<0,a1=a1+(2*pi);end
    a3 = atan2(yc-y3,x3-xc);
    if a3<0,a3=a3+(2*pi);end
    
    points_r = a1:(a3-a1)/(nump-1):a3;
    if sum(points_r>pi)>0, points_r(points_r>pi) = points_r(points_r>pi/2)-(2*pi); end
    out = [radius*cos(points_r)+xc; radius-(radius*sin(points_r))+(yc-radius)];
end


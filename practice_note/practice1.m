x = [-3/4*pi -pi/2 0 pi/4 5/4*pi 3/2*pi];
standard = pi/4;
x = -10:0.1:10;
y = wrapToPi(x - standard) + standard;
plot(x,y)
% Grating compressor dispersion calculator

c=0.3;

D=1.48; % lines per um
lambda=0.773; %central wavelength
omega=(2*pi*c/lambda); %central angular frequency
theta1=50.5; %incident angle on first grating (degrees) original 47.8

% calculate 1st order diffraction angle:
theta2=asin(sin(theta1*pi/180)-lambda*D)*180/pi;

GVD=-2*(1/cos(theta2*pi/180)^2)*(4*pi^2*c*D^2)/omega^3*1000; % compressor GVD in fs^2/mm
TOD=2*(12*pi^2*c*D^2/omega^4)*(cos(theta2*pi/180)^-2+(sin(theta2*pi/180)*2*pi*c*D/(omega*cos(theta2*pi/180)^4)))*1000; %compressor third order dispersion in fs^3/mm

disp(['GVD: ',num2str(GVD),' fs^2/mm'])
disp(['TOD: ',num2str(TOD),' fs^3/mm'])


% angle dependence of dispersion:
thetaArray=(theta1-10):0.01:(theta1+10);
theta2Array=asin(sin(thetaArray*pi/180)-lambda*D)*180/pi;
GVDarray=-2*(1./cos(theta2Array*pi/180).^2)*(4*pi^2*c*D^2)/omega^3*1000;
TODarray=2*(12*pi^2*c*D^2/omega^4)*(cos(theta2Array*pi/180).^-2+(sin(theta2Array*pi/180)*2*pi*c*D./(omega*cos(theta2Array*pi/180).^4)))*1000;

func1=1+(GVD-GVDarray)./GVDarray;
func2=TOD./TODarray;

figure('Name','Second and third order compressor dispersions vs. incident angle')
plot(thetaArray,GVDarray,'r-',thetaArray,TODarray,'b-')


% % I don't trust that my calculations in the code below are correct, so
% % use it at your own peril
% figure('Name','red curve < blue curve in order to lower TOD of pulse')
% plot(thetaArray,func1,'r-',thetaArray,func2,'b-')
% 
% 
% % determine how to adjust grating angle and distance in order to reduce TOD
% % of pulse while also maintaining GDD = 0
% L=475; % distance between gratings (mm)
% phiTripPrime=1.08*10^5; %third order pulse dispersion that you need to get rid of (fs^3)
% DeltaL=(L*(TOD-TODarray)-phiTripPrime)./TODarray;
% DeltaL2=L*((GVD-GVDarray)./GVDarray);
% 
% figure('Name','Intersection gives angle and change in compressor distance to get rid of TOD of pulse')
% plot(thetaArray,DeltaL,thetaArray,DeltaL2)
% 
% % calculate where the two curves intersect:
% [~, index]=min(abs(DeltaL-DeltaL2));
% NewTheta=thetaArray(index);
% changeInL=mean([DeltaL(index),DeltaL2(index)]);
% 
% disp(['to eliminate ',num2str(phiTripPrime),' fs^3 of pulse dispersion,'])
% disp(['adjust grating incident angle from ',num2str(theta1),' to ',num2str(NewTheta),' degrees, and'])
% disp(['change distance between gratings by ',num2str(changeInL),' mm.'])
% disp(['GVD at new position: ',num2str(GVDarray(index)),' fs^2/mm'])
% disp(['TOD at new position: ',num2str(TODarray(index)),' fs^3/mm'])

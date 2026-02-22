clear all
close all

% Process D-scan measurement into a common D-scan trace format to input
% into retrieval algorithm codes. 
% 
% Version 2 is designed to take d-scan traces output by the automated
% LabVIEW code. 
%
% This version calibrates input spectra with Bright and Dark calibration 
% spectra from the LabVIEW spectrometer program.
%
% The expected input is a single matrix with the following format:
% 
% 0  l1  l2  l3  l4 ... (one zero followed by the wavelength array)
% 0  f1  f2  f3  f4 ... (one zero followed by the fundamental spectrum)
% p1 s11 s12 s13 s14 ... (initial stage pos; first row of raw d-scan)
% p2 s21 s22 s23 s24 ... (second stage pos; second row of d-scan)
% .  .
% .    .
% .      .
%
% This version of the code CALIBRATES ALL INPUT SPECTRA useing calibration
% spectra.
%
% Each collected spectrum needs to be 
% 1) background-subtracted
% 2) plotted against frequency (instead of wavelength) and resampled for
% uniform sampling in frequency space
% 3) truncated to include only the SH, centering the frequency on that
% which corresponds to the central frequency of the fundamental spectrum
%
% The final output will be a .txt file with a matrix for the body,
% representing the D-scan trace (mm of dispersive material/compressor
% grating distance vs. SH frequency).  The first row of the matrix will
% correspond to the fundamental spectrum.
% The header will specify the central SH frequency (f0SH), frequency spacing (df), 
% starting position of dispersive material/stage position (Pi), spacing of 
% positions (dP):
%
% output file format:
% header:      f0SH, df, Pi, dP,
% fund spec:    s1, s2, s3, ... ,   sN
% Dscan trace: S11, S12, S13, ... , S1N
%              S21, S22, S23, ... , S2N
%               .    .    .     .    .
%               .    .    .     .    .
%              SM1, SM2, SM3, ... , SMN
%
% NxM matrix with N spectral samples and M stage position samples.
%
% The output file is TAB-DELIMITED

%D-scan data and background spectrum:
path = '/Users/conradkuz/Downloads/dscanMatlab/';
fn='2021-07-20_25_30.txt';
bgFN='2021-07-20_25_30_background.txt';

dat=dlmread([path,fn],'\t');

bgDat=dlmread([path,bgFN],'\t');
bgRaw=bgDat(2,:);

% File path for calibration spectra:
calibPath='/Users/conradkuz/Downloads/dscanMatlab/';
DarkFN='background22021-06-25.txt';
BrightFN='ls2-1-2021-06-28.txt';

% extract stage position, wavelength, fundamental spectrum, and dscan
% arrays
pos=dat(3:end,1);
lambda=dat(1,2:end);
fundspecRaw=dat(2,2:end);
dscanRaw=dat(3:end,2:end);

% Compute corresponding frequency array:
rawFreqs=300./lambda; % PHz

% calculate average frequency spacing to be used as the sampling rate for
% the uniformly resampled frequency axis
df=mean(rawFreqs(1:end-1)-rawFreqs(2:end));
uniFreqs=min(rawFreqs):df:max(rawFreqs);

% Calibration spectra: Bright (BlackBody Lamp) and Dark (electronic noise)
% spectra, and theoretical blackbody spectrum:
brightDat=dlmread([calibPath,BrightFN],'\t');
bright=brightDat(2,:);
darkDat=dlmread([calibPath,DarkFN],'\t');
dark=darkDat(2,:);

h=6.626*10^-34; %Planck's Constant
T=3100; %Lamp temp
k = 1.38*10^-23; %Boltzmann constant
c=300; %speed of light (nm/fs)
BBspec=2*h*(c*10^6)^2./((lambda*10^-9).^5.*(exp(h*(c*10^6)./((lambda*10^-9)*k*T)) - 1));

% Calibrate Fundamental Spectrum and d-scan trace:
bg=BBspec.*((bgRaw-dark)./(bright-dark))/10^10;
%fundspec=BBspec.*((fundspecRaw-dark)./(bright-dark))/10^10;
fundspec=fundspecRaw;
fundspec=fundspec-bg;
%dscan=BBspec.*((dscanRaw-dark)./(bright-dark))/10^10;
dscan=dscanRaw;
dscan=dscan-bg;

% Resample over uniform frequency grid:
RStrace1=zeros(length(pos),length(lambda));
for i=1:length(pos)
    F=griddedInterpolant(flip(rawFreqs),flip(dscan(i,:)));
    RStrace1(i,:)=300./(uniFreqs.^2).*F(uniFreqs);
end
F=griddedInterpolant(flip(rawFreqs),flip(fundspec));
RSfundSpec=300./(uniFreqs.^2).*F(uniFreqs);

% pad resampled spectra on the left with zeros all the way to the
% frequency origin:
extraFreqs=(0:df:uniFreqs(1)-df);
uniFreqs=[extraFreqs,uniFreqs];
RStrace=zeros(length(pos),length(uniFreqs));
for i=1:length(pos)
    RStrace(i,:)=[zeros(1,length(extraFreqs)),RStrace1(i,:)];
end
RSfundSpec=[zeros(1,length(extraFreqs)),RSfundSpec];


% suppress high noise at long wavelengths due to calibration
cutoffF=0.33;
width=0.01;
filter=0.5*(1+erf((uniFreqs-cutoffF)/width));
%RSfundSpec=filter.*RSfundSpec;
for i=1:length(pos)
    RStrace(i,:)=filter.*RStrace(i,:);
end
RSfundSpec=sum(RStrace);
RSfundSpec=RSfundSpec/max(RSfundSpec);

figure(1)
plot(uniFreqs(711:1300),RSfundSpec(711:1300))
%511:871 moved to 781:1300
% Find weighted average fundamental frequency:
fundField=sign(RSfundSpec).*sqrt(abs(RSfundSpec));
freq0fund=sum(uniFreqs(711:1300).*fundField(711:1300))/sum(fundField(711:1300));
[~, indf0]=min(abs(uniFreqs-freq0fund));

% Use twice the central frequency of the fundamental spectrum as the
% central frequency of the SH spectrum
[val, ind] = min(abs(uniFreqs-freq0fund*2));
f0SH=uniFreqs(ind);

% choose number of samples to truncate about the central SH frequency:
Nsamps=floor(2^8);
finalFreqArray=uniFreqs(ind-Nsamps/2+1:ind+Nsamps/2)-f0SH;
finalTrace=RStrace(:,ind-Nsamps/2+1:ind+Nsamps/2);
finalFundSpec=RSfundSpec(indf0-Nsamps/2+1:indf0+Nsamps/2);

% Plot d-scan trace:
figure('Position',[250,200,500,400])
imagesc(finalFreqArray,pos,finalTrace);
xlabel('SH frequency (PHz)'); ylabel('Compressor Stage Position (mm)'); zlabel('SH signal');
colormap('jet');
shading interp
set(gca,'YDir','normal')

% Plot Fundamental Spectum:
figure('Position',[780,200,500,400])
plot(finalFreqArray,finalFundSpec)


% Save data to output file (SEE TOP OF CODE FOR FILE FORMAT DETAILS)
initialStagePos=pos(1);
ScanIncrement=pos(2)-pos(1);

str=input('Save data to disc? y/n: ','s');
if str == 'n'
    disp('rerun program with new settings');
elseif str == 'y'
    header=[f0SH, df, initialStagePos, ScanIncrement];
    fid=fopen([path,fn(1:end-4),'_Processed.txt'],'w');
    fprintf(fid, '%f\t',header);
    tline=fgets(fid);
    fclose(fid);
    
    body=cat(1,finalFundSpec,finalTrace);
    dlmwrite([path,fn(1:end-4),'_Processed.txt'],body,'-append','delimiter','\t','roffset',1);
    
else
    disp('incorrect input.  Rerun the program.');
end



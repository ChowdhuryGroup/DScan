% THIS VERSION OF THE CODE CAN BE USED TO RETREIVE PULSES FROM MEASURED
% D-SCANS from the lab.
%
% D-scan retreival algorithm based on numerical minimization of the fitness
% function G (rms difference between measured and simulated D-scan traces).
% The algorithm is based on differential evolution, outlined in Das 2011.
% This version of the algorithm is based more closely on that described in
% Escoto 2017:
% - the DE part of the algorithm is identical to that used in DScan_DE_2.m
% (i.e. DE\BoR\1\bin, and with F = rand(0,1) for each donor vector)
% - the main difference is that the resolution of the spectral phase begins
% low and is gradually increased until the sampling limit is reached. For
% each resolution, a spline interpolation between the DE phase vector
% elements is used to calculate the spectral phase for each vector with the
% full resolution of the original spectrum in order to calculate a D-scan
% trace and evaluate the G value.
%
% Version 5 Updates:
% - Takes data from LabVIEW measurement
% - uses Sellmeier equation for BK7 (from refractiveindex.info) in order to
% generated d-scan traces (zeroth and first order dispersions subtracted)
% - implements frequency-dependent mu-factor to take into account
% nonuniformities in the spectral response of the measurement

close all;
clear all;

% Open d-scan data file:
filepath='/Users/conradkuz/Downloads/dscanMatlab/';
filename='2021-07-20_25_30_Processed.txt';

allDat=dlmread([filepath,filename]);

% Extract central SH frequency, calculate central fund freq, extract frequency spacing, initial stage pos, stage
% increment, and number of frequency and position samples
f0SH=allDat(1,1);
fcenter=f0SH/2;
df=allDat(1,2);
initialStagePos=allDat(1,3);
dPos=allDat(1,4);

Nsamps=length(allDat(2,:)); % number of frequency samples
NzSamps=length(allDat(:,1)')-2; % number of z pos samples

% define frequency and position arrays.  These arrays are really "offset"
% arrays, where the middle element of each array is zero.
freqs=([1:1:Nsamps]-Nsamps/2)*df;
z=([1:1:NzSamps]-NzSamps/2)*dPos; % inches of stage travel

% time array:
t=([1:1:Nsamps]-Nsamps/2)*(1/(2*max(freqs)));

% Fundamental spectrum:
fundSpec=allDat(2,:);
Espec=sign(fundSpec).*sqrt(abs(fundSpec)); % spectrum of the field (assuming flat phase)
Espec=Espec.';

% NOTE: in the algorithm, each time the d-scan trace is computed for a
% target or trial vector, the retreived field spectrum is set equal to the
% absolute value of the measured field spectrum (i.e. Ertrvd = abs(Espec)).
% However, because there is noise in the measured spectrum, the absolute
% value signs SHOULD BE REMOVED so that there is not a nonzero background
% level on average.  The Espec array in this version of the code is always
% associated with the measured spectrum of the fundamental and carries no
% phase information.

% Measured D-scan trace:
dScanMeasd=allDat(3:end,:); % (z on vertical, freq on horizontal)
dScanMeasd=dScanMeasd/max(dScanMeasd,[],'all');

% recenter the position axis so that z=0 mm corresponds to the position
% where the position-marginal of the d-scan trace is a maximum:
[~, zMargMaxInd]=max(sum(dScanMeasd,2));
zMargMax=z(zMargMaxInd);

% frequency and z meshgrids.  NOTE that these meshgrids correspond to a
% d-scan matrix that has frequency on the vertical and z on the horizontal.
% This is because the fft() function performs the fft along columns (freq.
% Also, the Newport BK7 wedges have an 8 degree wedge angle, which is
% factored into zmat to calculate mm of glass
freqmat=repmat(freqs.',1,length(z));
zmat=repmat(z-zMargMax,length(freqs),1); % recentered as explained above.
zmat=zmat*tan(8*pi/180)*25.4; % mm of glass, taking into account wedge angle

zGlass=zmat(1,:); % array of mm of glass, to be used later for plotting against glass insertion instead of stage position

[freqGrid,zGrid]=meshgrid(freqs,z);
figure('Name','Input D-scan Trace','Position',[100 50 400 325])
imagesc(freqs,z-zMargMax,dScanMeasd)
colormap('jet');
shading interp
set(gca,'YDir','normal')
xlim([min(freqs) max(freqs)])

dScanMeasd=circshift(dScanMeasd,Nsamps/2-1,2);


% Sellmeier coefficients for BK7 (from refractiveIndex.info)
A1=1.03961212;
B1=-0.00600069867;
A2=0.231792344;
B2=-0.0200179144;
A3=1.01046945;
B3=-103.560653;

% BK7 index of refraction MINUS zeroth and first order phase terms to
% subtract the group delay and phase delay at the central frequency caused
% by propagation through the glass.
nBK7full=sqrt(1+(A1./(1+B1.*((fcenter+freqs)/0.3).^2))+(A2./(1+B2.*((fcenter+freqs)/0.3).^2))+(A3./(1+B3.*((fcenter+freqs)/0.3).^2)));
nBK7_0= sqrt(1+(A1./(1+B1.*(fcenter/0.3).^2))+(A2./(1+B2.*(fcenter/0.3).^2))+(A3./(1+B3.*(fcenter/0.3).^2)));
dndf_0=(1/nBK7_0)*(-A1*B1*(fcenter/0.3^2)/(1+B1*(fcenter/0.3)^2)^2-A2*B2*(fcenter/0.3^2)/(1+B2*(fcenter/0.3)^2)^2-A3*B3*(fcenter/0.3^2)/(1+B3*(fcenter/0.3)^2)^2);

% Newport BK7 wedges dispersion phase matrix
phaseMat=exp(1i*2*pi*zmat.*((fcenter+freqmat).*(nBK7full')-(fcenter*nBK7_0+(nBK7_0+fcenter*dndf_0)*(freqmat)))/(3*10^-4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin d-scan retrieval algorithm: Differential Evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set coarsest sampling rate for first iteration of spectral resolution.
% At every stage of increase of resolution, the number of samples will
% increase by a factor of 2 until the original sampling rate is reached.
N1=2^7; % (i.e. take every N1-th phase) orignal - 2^4
indArray1=1:N1:Nsamps;

NresIncreases=log2(N1);

% Define search window.  Since spectral phases are defined from -pi to pi,
% we can just make that the range of the initial targets phases.
phimin=-pi;
phimax=pi;

% Population size:
NP=400; %original - 40
% initialize generation zero, represented by a NPxN1 matrix, where each row vector is a target vector:
X=phimin+rand(NP,Nsamps/N1)*(phimax-phimin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUTATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a donor vector associated with each target vector:
% Note that the mutation scaling factor, F, is a randomly generated number
% between 0 and 1 for every donor vector (see last line in first nested for
% loop below)

% initialize donor vector matrix:
V=zeros(NP,Nsamps/N1);

% interpolate phases with full resolution of original spectrum:
Xrs=interp1(freqs(indArray1),X.',freqs,'spline',0); % this is an NsampsxNP matrix now
    
NitersPerRes = 300; %This is the number of iterations to make for every spectral resolution used. Original - 300
Gvals=ones(NitersPerRes*NresIncreases,1);
for q=0:1:NresIncreases
    indArray1=1:N1/(2^q):Nsamps;
    X=Xrs(indArray1,:).';
    V=zeros(NP,Nsamps/(N1/(2^q)));
    for j=1:1:NitersPerRes

        for i=1:1:NP
            % 1) for each row of X (i.e. each target vector), randomly choose three
            % other target vectors
            mutInd=randperm(NP);
            mutInd(mutInd==i)=[];
            mutInd=mutInd(1:3);

            % 2) calculate G-value for the chosen target vectors

            G3rand=ones(3,1);
            for m=1:1:3
               Ertrvd=Espec; % column vector
               Ertrvd=Ertrvd.*exp(1i*Xrs(:,mutInd(m)));
               Ematrtrvd=repmat(Ertrvd,1,length(z));
               Ematrtrvd=Ematrtrvd.*phaseMat;
               ft=fft(Ematrtrvd);
               shg=ft.^2;
               shgSpec=ifft(shg);
               dScanRtrvd=(abs(shgSpec).^2).';

               mu=sum(dScanMeasd.*dScanRtrvd)./sum(dScanRtrvd.^2); % This is an array of mu values associated with each frequency
               G3rand(m)=sqrt((1/(length(freqs)*length(z)))*sum((dScanMeasd-mu.*dScanRtrvd).^2,'all')); 
            end

            % 3) sort G-values of the three target vectors in ascending order
            % of fitness (i.e. if v=x1 + rand(0,1)*(x2-x3), then
            % G(x1)<G(x2)<G(x3))
            [G3randSorted, G3randSortInd]=sort(G3rand);

            % 4) for each target vector, associate a donor vector v constructed
            % from the three other randomly chosen, sorted targets.
            V(i,:)=X(mutInd(G3randSortInd(1)),:)+rand.*(X(mutInd(G3randSortInd(2)),:)-X(mutInd(G3randSortInd(3)),:));
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % CROSSOVER:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generate trial vectors by having there be a chance that components of a
        % given target vector are swapped with those of the associated donor
        % vector.
        % initialize trial vector matrix to be equal to the target vector matrix:
        U=X;
        % set crossover rate:
        Cr = 0.5;
        % initialize crossover event outcomes (rand[0,1] to compare to Cr to decide
        % whether crossover for a given vector component should occur
        crRand=rand(NP,Nsamps/(N1/(2^q)));
        % logical array, where each element indicates whether crossover should or
        % shouldn't be done for the corresponding vector element (swYes = "switch?
        % Yes.")
        swYes=crRand<=Cr;

        % apply crossover for elements that passed the crossover test (for swYes =
        % 1)
        U(swYes)=V(swYes);
        % For elements that do not crossover, the trial vector element should be set equal to
        % the target vector element, and it already is because we initialized the
        % trial vector matrix to be equal to the target vector matrix.

        % the crossover process should also be done for a randomly selected element
        % in each vector (so that if the crossover test fails for all elements of a
        % given vector, some amount of crossover still occurs). So we randomly
        % choose an element to automatically crossover for each vector:
        jrand=randi(Nsamps/(N1/(2^q)),NP,1);
        % use linear indexing to pick out each element that is to be automatically
        % crossed over:
        jrandLinIndx=NP.*jrand-NP+(1:NP)';
        U(jrandLinIndx)=V(jrandLinIndx);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SELECTION:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Replace target vectors with the corresponding trial vectors if the
        % function value of the trial vector is less than that of the target
        % vector:

        % evaluate the function at each target and trial:

        % 1) initialize lists GX and GU of the G-values for the target vectors
        % X and trial vectors U:
        GX=ones(NP,1);
        GU=ones(NP,1);

        % 2) Calculate G value for each target vector:
        for m=1:1:NP
           Ertrvd=Espec; % column vector
           Ertrvd=Ertrvd.*exp(1i*Xrs(:,m));
           Ematrtrvd=repmat(Ertrvd,1,length(z));
           Ematrtrvd=Ematrtrvd.*phaseMat;
           ft=fft(Ematrtrvd);
           shg=ft.^2;
           shgSpec=ifft(shg);
           dScanRtrvd=(abs(shgSpec).^2).';

           mu=sum(dScanMeasd.*dScanRtrvd)./sum(dScanRtrvd.^2);
           GX(m)=sqrt((1/(length(freqs)*length(z)))*sum((dScanMeasd-mu.*dScanRtrvd).^2,'all')); 
        end

        % 3) Calculate G value for each trial vector:
        % interpolate phases first:
        Urs=interp1(freqs(indArray1),U.',freqs,'spline',0); % this is an NsampsxNP matrix now
        for m=1:1:NP
           Ertrvd=Espec; % column vector
           Ertrvd=Ertrvd.*exp(1i*Urs(:,m));
           Ematrtrvd=repmat(Ertrvd,1,length(z));
           Ematrtrvd=Ematrtrvd.*phaseMat;
           ft=fft(Ematrtrvd);
           shg=ft.^2;
           shgSpec=ifft(shg);
           dScanRtrvd=(abs(shgSpec).^2).';

           mu=sum(dScanMeasd.*dScanRtrvd)./sum(dScanRtrvd.^2);
           GU(m)=sqrt((1/(length(freqs)*length(z)))*sum((dScanMeasd-mu.*dScanRtrvd).^2,'all')); 
        end

        [Gvals(j+q*NitersPerRes), index1]=min([GX;GU]);
        index=mod(index1(1)-1,NP)+1;

        % logical array indicating whether trial vector has lower function value
        % than target:
        fUlsthnfX=GU<=GX;

        % replace target vectors with trial vectors that have lower function value:
        X(fUlsthnfX,:)=U(fUlsthnfX,:);
        Xrs(:,fUlsthnfX)=Urs(:,fUlsthnfX);
    end
    
end


% plot G values:
figure('Name','G values')
plot([1:1:NitersPerRes*(NresIncreases+1)],log10(Gvals),'.b')

% Display lowest G value:
disp(['G = ',num2str(Gvals(end))])

% Calculate D-scan trace of best guess:
Ertrvd=Espec; % column vector
Ertrvd=Ertrvd.*exp(1i*Xrs(:,index));
Ematrtrvd=repmat(Ertrvd,1,length(z));
Ematrtrvd=Ematrtrvd.*phaseMat;
ft=fft(Ematrtrvd);
shg=ft.^2;
shgSpec=ifft(shg);
dScanRtrvd=(abs(shgSpec).^2).';
dScanRtrvd=dScanRtrvd/max(dScanRtrvd,[],'all');
dScanRtrvd=circshift(dScanRtrvd,Nsamps/2,2);

% Calculate spectral response function based on frequency-dependent mu:
dScanMeasd=circshift(dScanMeasd,Nsamps/2,2);
mu=sum(dScanMeasd.*dScanRtrvd)./sum(dScanRtrvd.^2);
mu=mu/max(mu);
figure('Name','Spectral response function')
plot(freqs,mu,'r-')

% Plot Best guess full D-scan Trace (i.e. removing spectral response nonuniformity):
figure('Name','Retrieved D-scan trace (full)')
imagesc(freqs,zGlass,dScanRtrvd)
colormap('jet');
shading interp
set(gca,'YDir','normal')
xlim([min(freqs) max(freqs)])

% Plot Best guess D-scan Trace with retreived nonuniform spectral response applied:
figure('Name','Retrieved D-scan trace (w\ spectral response)')
imagesc(freqs,zGlass,mu.*dScanRtrvd)
colormap('jet');
shading interp
set(gca,'YDir','normal')
xlim([min(freqs) max(freqs)])

% calculate temporal intensity distributions for each stage position:
Imatrtrvd=abs(ft).^2;
Imatrtrvd=Imatrtrvd.';
Imatrtrvd=circshift(Imatrtrvd,Nsamps/2+1,2);

% Find time associated with largest peak intensity, and recenter intensity
% distribution so that the peak intensity lies at t=0 (i.e. get rid of
% linear component of spectral phase:
[maxI, Iind]=max(Imatrtrvd,[],'all','linear');
t0Ind=floor((Iind+1)/length(z));
t0=t(t0Ind+1);
zeroOriginPhaseMat=exp(1i*2*pi*freqmat*(-t0));
Ematrtrvd=Ematrtrvd.*zeroOriginPhaseMat;
Ematrtrvd=Ematrtrvd*(max(Espec)/max(abs(Ematrtrvd),[],'all'));
ft=fft(Ematrtrvd);
Imatrtrvd=abs(ft).^2;
Imatrtrvd=Imatrtrvd.';
Imatrtrvd=circshift(Imatrtrvd,Nsamps/2,2);

% % Plot temporal distribution (vertical axis z, horizontal axis t)
% figure('Name','Retrieved intensity distributions','Position',[900 50 400 325])
% imagesc(t,z,Imatrtrvd)
% colormap('jet');
% shading interp
% set(gca,'YDir','normal')
% xlim([-25 25])

% Calculate BW-limited pulse intensity envelope:
ftBWL=fft(Espec);
ftBWL=circshift(ftBWL,Nsamps/2);
BWLintensity=abs(ftBWL).^2;

% Display z position corresponding to highest peak intensity compression:
Iind=mod(Iind-1,length(z))+1;
disp(['Stage position for best compression: ',num2str(z(Iind)+(NzSamps/2-1)*dPos+initialStagePos),' inches'])
disp(['(z = ',num2str(zGlass(Iind)),' mm on plot)'])

% Determine at what time the shortest retrieved pulse has it's peak
% intensity
ftlr=fft(Ematrtrvd);
ftlr=circshift(ftlr,Nsamps/2,1);
[~,maxIndlowRes]=max(abs(ftlr(:,Iind)).^2);
tlr=t(maxIndlowRes);

% Compute the spectral phase associated with the stage position that gives
% the largest peak intensity (i.e. most compressed pulse):
bestPhase=unwrap(Xrs(:,index)); 
bestPhase=bestPhase.'+ 2*pi*freqs*(-tlr)+2*pi*zGlass(Iind)*((fcenter+freqs).*(nBK7full)-(fcenter*nBK7_0+(nBK7_0+fcenter*dndf_0)*(freqs)))/(3*10^-4);
bestPhase=bestPhase-bestPhase(Nsamps/2)+pi;

% Plot flattest spectral phase for retrieved d-scan:
figure('Name','Flattest Retrieved Phase')
plot(freqs,bestPhase/(2*pi),'-r',freqs,fundSpec/max(fundSpec),'-b')
xlim([-0.12 0.12])
ylim([-0.3 1.3])

% Plot original and retrieved pulses in time-domain, but first pad the
% spectrum with zeros so that time-resolution is increased
RSNsamps=Nsamps*2^5;
RSfreqs=([1:1:RSNsamps]-RSNsamps/2)*df;
Nnewpoints=length(RSfreqs)-length(freqs);

EmatrtrvdRS=cat(1,Ematrtrvd,zeros(Nnewpoints/2,length(z)));
EmatrtrvdRS=cat(1,zeros(Nnewpoints/2,length(z)),EmatrtrvdRS);

EOrig=cat(1,Espec,zeros(Nnewpoints/2,1));
EOrig=cat(1,zeros(Nnewpoints/2,1),EOrig);

time=([1:1:RSNsamps]-RSNsamps/2)*(1/(2*max(RSfreqs)));

ft2=fft(EOrig);
ft2=circshift(ft2,RSNsamps/2);
BWLintRS=abs(ft2).^2;
[BWLmax, maxInd0]=max(abs(ft2).^2);
t0=time(maxInd0);
BWLintRS=BWLintRS/BWLmax;

ft3=fft(EmatrtrvdRS);
ft3=circshift(ft3,RSNsamps/2,1);
[~, maxInd]=max(abs(ft3(:,Iind)).^2);
tr=time(maxInd);

EmatrtrvdRS=EmatrtrvdRS.*(exp(1i*2*pi*RSfreqs'.*(-tr+t0)));
ft4=fft(EmatrtrvdRS);
ft4=circshift(ft4,RSNsamps/2,1);

ImatrtrvdRS=abs(ft4).^2;
ImatrtrvdRS=(ImatrtrvdRS.')/BWLmax;

% plot higher res intensity envelope vs time for both BW-limited and
% highest peak intensity retrieved pulses:
figure('Name','Intensity Envelope vs Time (higher resolution)')
plot(time,BWLintRS,'-b',time,ImatrtrvdRS(Iind,:),'-r')
xlim([-100 100])

% Plot high temporal resolution intensity vs time and z distribution (vertical axis z, horizontal axis t)
figure('Name','Retrieved intensity distributions','Position',[900 50 400 325])
imagesc(time,zGlass,ImatrtrvdRS)
colormap('jet');
shading interp
set(gca,'YDir','normal')
xlim([-45 45])


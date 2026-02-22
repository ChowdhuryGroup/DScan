% THIS VERSION OF THE CODE CAN BE USED TO RETREIVE PULSES FROM MEASURED
% D-SCANS from the lab.
%
% This version uses calculated values for the compressor second and third
% order dispersion (i.e. not "self-calibrating").
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
%n
% NOTE: Be sure to make note of what direction corresponds to increasing
% grating spacing with respect to the stage micrometer.  Currently,
% INCREASING z corresponds to DECREASING grating separation, meaning that
% the GVD of the compressor with respect to z is POSITIVE.  This means that
% the initial range of guess GVD values should bound a set of positive
% numbers.
%
% Updates in version 2:
% uses a frequency-dependent mu factor to take into account non-uniformity
% in spectral response.

close all;
clear all;

% Open d-scan data file:
filepath='/Users/conradkuz/Documents/dscanMatlab/';
filename='11.085mm dscan_Processed_NoShift.txt';

% GRATING COMPRESSOR DISPERSION (FROM
% GratingCompressorDispersionCalculator.m):
GVD=-4154.4; %fs^2/mm original -4270
TOD=2582.8; %fs^3/mm original 2360

allDat=dlmread([filepath,filename]);

% Extract central frequency, frequency spacing, initial stage pos, stage
% increment, and number of frequency and position samples
f0SH=allDat(1,1);
df=allDat(1,2);
initialStagePos=allDat(1,3);
dPos=allDat(1,4);

Nsamps=length(allDat(2,:)); % number of frequency samples
NzSamps=length(allDat(:,1)')-2; % number of z pos samples

% define frequency and position arrays.  These arrays are really "offset"
% arrays, where the middle element of each array is zero.
freqs=([1:1:Nsamps]-Nsamps/2)*df;
z=([1:1:NzSamps]-NzSamps/2)*dPos;

% time array:
t=([1:1:Nsamps]-Nsamps/2)*(1/(2*max(freqs)));

% Fundamental spectrum:
fundSpec=allDat(2,:);
Espec=sign(fundSpec).*sqrt(abs(fundSpec)); % spectrum of the field (assuming flat phase)
Espec=Espec.';

% Calculate bandwidth of fund spec (to compute range over which phase fit
% is to be done)
Df=sqrt(sum(freqs.^2.*fundSpec)/sum(fundSpec)-sum(freqs.*fundSpec)^2/sum(fundSpec)^2);

% NOTE: in the algorithm, each time the d-scan trace is computed for a
% target or trial vector, the retreived field spectrum is set equal to the
% absolute value of the measured field spectrum (i.e. Ertrvd = abs(Espec)).
% However, because there is noise in the measured spectrum, the absolute
% value signs SHOULD BE REMOVED so that there is not a nonzero background
% level on average.  The Espec array in this version of the code is always
% associated with the measured spectrum of the fundamental and carries no
% phase information.

% frequency and z meshgrids.  NOTE that these meshgrids correspond to a
% d-scan matrix that has frequency on the vertical and z on the horizontal.
% This is because the fft() function performs the fft along columns (freq)
freqmat=repmat(freqs.',1,length(z));
zmat=repmat(z,length(freqs),1);


% Measured D-scan trace:
dScanMeasd=allDat(3:end,:); % (z on vertical, freq on horizontal)
dScanMeasd=dScanMeasd/max(dScanMeasd,[],'all');

[freqGrid,zGrid]=meshgrid(freqs,z);
figure('Name','Input D-scan Trace','Position',[100 50 400 325])
imagesc(freqs,z,dScanMeasd)
colormap('jet');
shading interp
set(gca,'YDir','normal')
xlim([min(freqs) max(freqs)])

dScanMeasd=circshift(dScanMeasd,Nsamps/2-1,2);

% grating compressor dispersion phase matrix:
phaseMat=exp(1i*zmat.*(0.5*GVD*((2*pi)*freqmat).^2+(1/6)*TOD*((2*pi)*freqmat).^3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Begin d-scan retrieval algorithm: Differential Evolution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set coarsest sampling rate for first iteration of spectral resolution.
% At every stage of increase of resolution, the number of samples will
% increase by a factor of 2 until the original sampling rate is reached.
N1=2^4; % (i.e. take every N1-th phase) original 2^4
indArray1=1:N1:Nsamps;

NresIncreases=log2(N1);

% Define search windows.  Since spectral phases are defined from -pi to pi,
% we can just make that the range of the initial targets phases.
phimin=-pi;
phimax=pi;

% Population size:
NP=80; %orignal 40
% initialize generation zero, represented by a NPxN1 matrix, where each row vector is a target vector:
X=phimin+rand(NP,Nsamps/N1)*(phimax-phimin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MUTATION:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a donor vector associated with each target vector:
% Note that the mutation scaling factor, F, is a randomly generated number
% between 0 and 1 for every donor vector (see last line in first nested for
% loop below)
% Donor vector initialized for every spectral resolution iteration.

% interpolate phases with full resolution of original spectrum:
Xrs=interp1(freqs(indArray1),X.',freqs,'spline',0); % this is an NsampsxNP matrix now
    
NitersPerRes = [50, 100, 200, 200, 250]; %This is the number of iterations to make for every spectral resolution used.
Gvals=ones(sum(NitersPerRes),1);
NitersSoFar=0;
Progress=0;
disp('Progress:');
disp([num2str(Progress),'%']);
progIncrement=floor(0.05*sum(NitersPerRes));
for q=0:1:NresIncreases
    indArray1=1:N1/(2^q):Nsamps;
    X=Xrs(indArray1,:).';
    V=zeros(NP,Nsamps/(N1/(2^q)));
    for j=1:1:NitersPerRes(q+1)

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

               mu=sum(dScanMeasd.*dScanRtrvd)./sum(dScanRtrvd.^2);
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
        Cr = 0.65;
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

        [Gvals(1+NitersSoFar), index1]=min([GX;GU]);
        index=mod(index1(1)-1,NP)+1;

        % logical array indicating whether trial vector has lower function value
        % than target:
        fUlsthnfX=GU<=GX;

        % replace target vectors with trial vectors that have lower function value:
        X(fUlsthnfX,:)=U(fUlsthnfX,:);
        Xrs(:,fUlsthnfX)=Urs(:,fUlsthnfX);
        
        
        % Progress counter:
        NitersSoFar=NitersSoFar+1;
        if mod(NitersSoFar,progIncrement)==0
            Progress=Progress+5;
            disp([num2str(Progress),'%'])
        else
            %do nothing
        end
    end
    
end


% plot G values:
figure('Name','G values')
plot([1:1:sum(NitersPerRes)],log10(Gvals),'.b')

%display final G value:
display(['Final G value: G = ',num2str(Gvals(end))])

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
mu=sum(dScanMeasd.*dScanRtrvd)./sum(dScanRtrvd.^2);
mu=circshift(mu,Nsamps/2+1,2);
dScanRtrvd=circshift(dScanRtrvd,Nsamps/2+1,2);

% spectral response:
figure('Name','Spectral response mu vs freq.')
plot(freqs,mu,'r-')

% Plot best guess D-scan trace:
figure('Name','Retrieved D-scan trace','Position',[500 50 400 325])
imagesc(freqs,z,mu.*dScanRtrvd)
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
Imatrtrvd=circshift(Imatrtrvd,Nsamps/2+1,2);

% Plot temporal distribution (vertical axis z, horizontal axis t)
figure('Name','Retrieved intensity distributions','Position',[900 50 400 325])
imagesc(t,z,Imatrtrvd)
colormap('jet');
shading interp
set(gca,'YDir','normal')

% Calculate BW-limited pulse intensity envelope:
ftBWL=fft(Espec);
ftBWL=circshift(ftBWL,Nsamps/2);
BWLintensity=abs(ftBWL).^2;

% Plot the intensity vs time corresponding to the largest peak intensity:
Iind=mod(Iind-1,length(z))+1;
figure('Name','Highest Intensity Retrieved pulse')
plot(t,Imatrtrvd(Iind,:),'-r',t,BWLintensity,'-b')
xlim([-250 250])

% Display z position corresponding to highest peak intensity compression:
zForBestCompression=z(Iind)+(NzSamps/2-1)*dPos+initialStagePos;
disp(['Stage position for best compression: ',num2str(zForBestCompression),' mm'])
disp(['(z = ',num2str(z(Iind)),' mm on plot)'])

% Compute the spectral phase associated with the stage position that gives
% the largest peak intensity (i.e. most compressed pulse):
bestPhase=unwrap(Xrs(:,index)); 
bestPhase=bestPhase.'+ 2*pi*freqs*(-t0)+z(Iind)*(0.5*GVD*((2*pi)*freqs).^2+(1/6)*TOD*((2*pi)*freqs).^3);
bestPhase=bestPhase-bestPhase(Nsamps/2)+pi;

% polynomial fit of spectral phase 
inRange=(abs(freqs)<=2.0*Df);
fit=polyfit(2*pi*freqs(inRange),bestPhase(inRange),6);
fitPhase=polyval(fit,2*pi*freqs(inRange));
disp(['Polynomial fit to phase: ','phi(w)= (1/4!)',num2str((factorial(4))*fit(3)),'w^4 + (1/6)',num2str(6*fit(4)),'w^3 + (1/2)',num2str(2*fit(5)),'w^2 + ',num2str(fit(6)),'w + ',num2str(fit(7))]);

% Plot flattest spectral phase for retrieved d-scan:
figure('Name','Flattest Retrieved Phase')
plot(freqs,bestPhase/(2*pi),'-r',freqs,fundSpec/max(fundSpec),'-b',freqs(inRange),fitPhase/(2*pi),'.k')
xlim([-0.02 0.02])
ylim([-0.2 1.2])

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
xlim([-250 250])

% Plot high temporal resolution intensity vs time and z distribution (vertical axis z, horizontal axis t)
figure('Name','Retrieved intensity distributions','Position',[900 50 400 325])
imagesc(time,z,ImatrtrvdRS)
colormap('jet');
shading interp
set(gca,'YDir','normal')


% Save to disk?
str=input('Save data to disc? y/n: ','s');
if str == 'n'
    disp('rerun program if desired');
elseif str == 'y'
    mkdir(filepath,[filename(1:end-4),'_Retrievals']);
    str=datestr(now);
    dateAndTimeStr=[date,'_',str(end-7:end-6),str(end-4:end-3)];    
    mkdir([filepath,filename(1:end-4),'_Retrievals'],dateAndTimeStr);
    savepath=[filepath,filename(1:end-4),'_Retrievals\',dateAndTimeStr,'\'];
    % save measured dscan trace:
    writematrix(circshift(dScanMeasd,Nsamps/2,2),[savepath,'Retrieved_',dateAndTimeStr,'_MeasuredDscanTrace.txt']);
    % save retrieved dscan trace:
    writematrix(dScanRtrvd,[savepath,'Retrieved_',dateAndTimeStr,'_RetrievedDscanTrace.txt']);
    % save matrix with rows: [central SH frequency, zeros(1,lenght(freqs)-1)],
    % freqs(offset), fundamental spectrum, retrieved phase, response
    % function mu
    writematrix([[f0SH,zeros(1,length(freqs)-1)];freqs;fundSpec;bestPhase;mu],[savepath,'Retrieved_',dateAndTimeStr,'_f0SH_f_spec_phi_mu.txt']);
    % save Resampled intensity vs time vs z matrix:
    writematrix(ImatrtrvdRS,[savepath,'Retrieved_',dateAndTimeStr,'_RetrievedResampledIvsZvsT.txt']);
    % save matrix with rows: resampled time, BW limited intensity, highest
    % peak intensity from retrieval
    writematrix([time;BWLintRS';ImatrtrvdRS(Iind,:)],[savepath,'Retrieved_',dateAndTimeStr,'_Resampled_time_BWLint_BestInt.txt']);
    % save G values:
    writematrix(Gvals,[savepath,'Retrieved_',dateAndTimeStr,'_Gvals.txt']);
    % save matrix with rows: z(offsets), [GVD, TOD,initialStagePos, stage position for highest peak intensity,zeros(1,length(z)-4)]
    writematrix([z;[GVD,TOD,initialStagePos,zForBestCompression,zeros(1,length(z)-4)]] ,[savepath,'Retrieved_',dateAndTimeStr,'_z_GVDandTODandPosiandBestZ.txt']);
    
    
else
    disp('incorrect input.  Rerun the program.');
end


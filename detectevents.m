% detectevents.m: 
% use definition-based approach to identify LP events in emergence from
% general anesthesia (Lewis et al., eLife 2018)
%
% Laura Lewis, Aug. 5 2018
% No warranty provided

% Input:
% chdata: channels x time, raw data
% filtdata: channels x time, should already be bandpassed between 0.2 and 4 Hz.
% fs: sampling rate
% chlabels: cell array of channel labels
% ploton: plotting flag

% output:
% kcevents: kc times in samples for all channels

function allkcs=detectevents(chdata,filtdata,fs,chlabels,ploton)

if nargin<4
    ploton=0;
end

% Set parameters for KC detection
ampthresh=400;
durthresh=0.4;
minthresh=40; % this is set as baseline value
maxthresh=1200; % reject if there are peaks above this
mingap=0.5;

%% then search for peaks and do detection in each channel
p=findpeaks(filtdata',ampthresh);
n=findpeaks(-filtdata',ampthresh);
nc=size(filtdata,1);

%% detect events and select polarity
allkcs=cell(nc,1);
nkcs=zeros(length(chlabels),1);
for c=1:nc
    kcs=[];
    signdata=filtdata(c,:);
    thresh=find(abs(filtdata(c,:))<minthresh);
    ev=[p(c).loc];
    ev=ev(ev>fs*2&ev<length(filtdata)-fs*2);
    for k=1:length(ev)
        currk=ev(k);
        ontime=find(thresh<currk);
        offtime=find(thresh>currk);
        
        if isempty(offtime)|isempty(ontime)
            continue
        end
        ontime=thresh(ontime(end));
        offtime=thresh(offtime(1));
        if offtime-ontime>durthresh*fs&all(abs(chdata(c,currk+(-fs*2:fs*2-1)))<maxthresh) % here, check duration and max amp
            kcs(end+1)=currk;
        end
    end
    
    signdata=-filtdata(c,:);
    ev=[n(c).loc];
    ev=ev(ev>fs*2&ev<length(filtdata)-fs*2);
    for k=1:length(ev)
        currk=ev(k);
        ontime=find(thresh<currk);
        offtime=find(thresh>currk);
        if isempty(offtime)|isempty(ontime)
            continue
        end
        offtime=thresh(offtime(1));
        ontime=thresh(ontime(end));
        
        if offtime-ontime>durthresh*fs&all(abs(chdata(c,currk+(-fs*2:fs*2-1)))<maxthresh) % here, check  duration
            kcs(end+1)=currk;
        end
    end
    
    kcs=sort(kcs);
    polarity=sign(filtdata(c,kcs));
    kcs=kcs(polarity==sign(median(filtdata(c,kcs))));
    if ~isempty(kcs)
        kcs=kcs([1000 diff(kcs)]>mingap*fs);
        allkcs{c}=kcs;
        nkcs(c)=length(kcs);
    end
end

%% plot timeseries
chtime=(1:length(chdata))/fs;
if ploton
    ci=1:nc;
    plotchannels(chdata(ci,:),fs,chlabels(ci),2);
    for i=1:length(ci)
        if ~isempty(allkcs(i))
            plot(chtime(allkcs{i}),i*ones(size(allkcs{i})),'r*');
        end
    end
end

%% triggered spectrogram

minkc=10; % require that there be at least a certain number of events
goodc=[];
for i=1:length(allkcs)
    if length(allkcs{i})>minkc
        goodc(end+1)=i;
    end
end
[~,ci]=max(nkcs); % select single channel of interest

kcevents=allkcs{ci}/fs;

% normalize
params=struct;
params.Fs=fs;
params.fpass=[0 100];
params.tapers=[2 3];
params.trialave=1;
win=[2 2];
movingwin=[0.2 0.05];
[s t f]=mtspecgramtrigc_detrend(chdata(ci,:),kcevents,win,movingwin,params);

normspec=mean(s(t<1,:),1);
normspec=repmat(normspec,length(t),1);

% mean timeseries
trigindx=-win(1)*fs:win(2)*fs;
trigdata=zeros(length(kcevents),length(trigindx));
trigt=trigindx/fs;

for i=1:length(kcevents)
    kcindex=round(kcevents(i)*fs);
    datawin=chdata(ci,kcindex+trigindx);
    trigdata(i,:)=datawin;
end

if ploton
    figure
    axs(1)=subplot(2,1,1);
    imagesc(t-win(1),f,log10(abs((s./normspec)'))); axis xy
    colorbar
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(chlabels(ci));
    caxis([-0.5 0.5]);
    
    axs(2)=subplot(2,1,2);
    plot(trigindx/fs,mean(trigdata)); colorbar
    linkaxes(axs,'x');
end

%% plot identified events

if ploton
    figure; plot(trigindx/fs,trigdata,'LineWidth',1.5);
    xlabel('Time (s)');
    ylabel('amplitude');
end


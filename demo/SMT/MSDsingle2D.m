function s_MSD=MSDsingle2D(trace,pixelsize,frameT,darkT,N)
%% calculate 2dMSD from a singel trajectory
%   Modified from Jie's code MSD_single

% input : trace : the single track (1, frame; 2, x; 3, y; 4, z;
% 5,intensity(or other info, not using)
%        pixelsize : the dimension of single pixel in um
%        frameT : the exposure time in s
%        darkT : the dark time between two frames in s
%        N  : the max number of time interval to calculate, smaller than
%        length of trace/2
%        CMode : use all the steps because of smaller sample size
% output :   s_MSD: MSD data
%                   column 1: time index
%                   column 2: MSD
%                   column 3: std of MSD
%                   column 4: SEM of MSD 


if nargin<4
    error('not enough input!');
    return;
end

L = size(trace,1);

if nargin < 5 | N >L/2
    N = round(L/2);
end
%%
for idxM = 1:N
    % calculate MSD
    MSD_temp=[]; % initial
    Coordinates=trace;
    frames=Coordinates(:,1)-Coordinates(1,1)+1;
    for idxQ=1:min((frames(end)-idxM),length(frames))
        Startf=frames(idxQ);
        Endf=find(frames==(Startf+idxM));
        if Endf
            Startp=Coordinates(idxQ,2:3);
            Endp=Coordinates(Endf,2:3);
            SquareD=sum(((Startp-Endp)*(pixelsize)).^2); % square displacement in um
            MSD_temp=[MSD_temp;SquareD];
        end
    end
    if isempty(MSD_temp)
        MSD(idxM,2:3)=NaN;
    else
        MSD(idxM,2)=mean(MSD_temp);
%         stats=bootstrp(1000, @mean, MSD_temp);
%         MSD(idxM,3)=std(stats);
            MSD(idxM,3)=std(MSD_temp);
            MSD(idxM,4)=std(MSD_temp)/sqrt(length(MSD_temp));
    end
end
Timebin=(1:N)*(frameT+darkT);
MSD(:,1)=Timebin';
MSD(~any(isnan(MSD),2),:);
s_MSD=MSD;

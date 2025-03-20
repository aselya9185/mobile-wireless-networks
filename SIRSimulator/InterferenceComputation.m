function [SIRup SIRdown] = InterferenceComputation(R,gamma,cellLoad,N,sectors,Rue,pcu)

%Simulator of signal to interference ratio
%Author: Andrea Detti

%R=1e3;  % cell radius
%gamma=4; % attenuation exponent
%sectors=3;
%cellLoad=1; % probability of using a resource in the cell
%N=1; % reuse factor (FCA)
%Rue = R; % distance of the reference ue from its bs
%pcu = 1; % power control uplink


BeamWidth=360/sectors; % antenna beam mesuered in degree
prx=10; % prx power set to 10 (not relevant)
trials = 1000; % number of sim trials
M=[[cos(pi/6) 0];[sin(pi/6) 1]]; %hexagonal to cartesian coordinate transformation Matrix
subchannelingPerSector=1; % subchanneling per sector


% create and deploy the (u,v) plane with base stations
[ u v firstTierFilter ] = DeployBS( R,N );

% Fixed Channel Allocation
F=FCA(u,v,N);


% draw cells and bs assigned channels
% figure;
% hold on;
% D=1.5*(sqrt(3*N)*R)+0.1;
% for i=1:length(u)
%     for j=1:length(v)
%         if ~ismember([i j],firstTierFilter,'rows')
%             continue
%         end
%         c=M*[u(i);v(j)];
%         circle(c(1),c(2),R,pi/3);
%         text(c(1),c(2),num2str(F(i,j)));
%     end
% end
% grid
% axis([-1.5*D 1.5*D -1.5*D 1.5*D]);


%
%uplink interference computation
%

ptot=zeros(trials,1);   % matrix storing interference power for each sim trial
for t=1:trials
    
    % deploy the reference ue in the reference cell (0,0) at distance d and
    % random angle (p)
    d=Rue;
    p=rand*(2*pi);
    ue_xr=d*cos(p); % reference ue cartesian coordinate x
    ue_yr=d*sin(p); % reference ue cartesian coordinate y
    bs_uir=find(u==0); % reference bs u coordinate index
    bs_vir=find(v==0); % reference bs v coordinate index
    angleBSrUEr=atan2d(ue_xr,ue_yr); % angle between the reference BS  and the reference ue
    if (angleBSrUEr<0)
        angleBSrUEr=angleBSrUEr+360;
    end
    ue_sr = ceil(angleBSrUEr/BeamWidth);% number of sector of the reference ue;
    %plot(ue_xr,ue_yr,'*r');
    
    %
    %for each cell random deploy a interfering ue (OFDMA/FDMA/TDMA systems)
    %
    for i=1:length(u)
        for j=1:length(v)
            if ~ismember([i j],firstTierFilter,'rows')
                % first tier approximation
                continue
            end
            if (i==bs_uir && j==bs_vir)
                continue;
            end
            
            % random deploy a interfering ue in the cell (i,j)
            d=rand(1,1)*R;
            p=rand(1,1)*2*pi;
            c=M*[u(i);v(j)];
            x=c(1)+d*cos(p);
            y=c(2)+d*sin(p);
            %plot(x,y,'+b');
            
            % computation of distances
            dother=sqrt(x^2+y^2);   % distance interfering ue to the bs(0,0)
            down=d;                 % distance interfering ue to the bs(i,j)
            
            %computation of tx power
            if (pcu==1)
                % tx power with power control, same signal power in bs antenna for each UE
                % same bitrate per UE
                ptx=prx*(down^(gamma));  % tx power so as arriving at the own bs with prx
                %log normal model with d0 equal to 1 and d0 attenuation equal to 1 (not relevant for SIR computation).
            else                
                % tx power without powercontrol, tx power so as achieving prx when at cell edge
                ptx=prx*(R^(gamma));
            end
            
            %worst case, interfering ue located in the direction bs(i,j)-->bs(0,0)
            %and at celle edge, traditional model S/I = sectors/6 (sqrt(3*N)-1)^gamma
            % to be used without powercontrol
            %xws=c(1);
            %yws=c(2);
            %dother=sqrt(xws^2+yws^2)-R;
            
            
            % interfering power at reference bs
            pti=ptx/(dother^(gamma)); % interfering power received at bs(0,0)
            
            % take into account antenna sectoring
            angleBSrUEi=atan2d(x,y); % angle between the reference BS and the interfering ue
            if (angleBSrUEi<0)
                angleBSrUEi=angleBSrUEi+360;
            end
            ue_si = ceil(angleBSrUEi/BeamWidth);% number of the sector of the interfering ue;
            if ue_si~=ue_sr
                %fprintf('different sector, no interference');
                pti=0;
                continue;
            end
            
            % taking into account cell load
            cl=rand(1,1);
            if cl>cellLoad
                %fprintf('no load, no interference');
                pti=0;
                continue;
            end
            
            % taking into account channel reuse
            if (F(i,j)~=F(bs_uir,bs_vir))
                %fprintf('different channels, no interference');
                pti=0;
                continue;
            end
            
            ptot(t)=ptot(t)+pti;
        end
    end
end
% average S/Iother, reference ue reaches reference bs with prx
if (pcu==1)
    % with power control
    ptx_uer=prx*Rue^(gamma);
    % no power control
else
    % without power control
    ptx_uer=prx*R^(gamma);
end

prx_uer=ptx_uer/Rue^(gamma);

SIRup = prx_uer/mean(ptot);
SIupdb=10*log10(SIRup);

%
%dowlink interference computation
%
ptot=zeros(trials,1);

for t=1:trials
    
    % deploy the reference ue in the reference cell (0,0) at distance d and
    % random angle (p)
    d=Rue;
    p=rand*(2*pi);
    ue_xr=d*cos(p); % reference ue cartesian coordinate x
    ue_yr=d*sin(p); % reference ue cartesian coordinate y
    bs_uir=find(u==0); % reference bs u coordinate index
    bs_vir=find(v==0); % reference bs v coordinate index
    angleBSrUEr=atan2d(ue_xr,ue_yr); % angle between the reference BS  and the reference ue
    if (angleBSrUEr<0)
        angleBSrUEr=angleBSrUEr+360;
    end
    ue_sr = ceil(angleBSrUEr/BeamWidth);% number of sector hosting the the reference ue;
    %plot(ue_xr,ue_yr,'*r');
    for i=1:length(u)
        for j=1:length(v)
            if ~ismember([i j],firstTierFilter,'rows')
                continue
            end
            if (i==bs_uir && j==bs_vir)
                continue;
            end
            
            
            c=M*[u(i);v(j)];
            x=c(1);
            y=c(2);
            
            % computation of distances
            dother=sqrt((x-ue_xr)^2+(y-ue_yr)^2);   % distance interfering bs(i,j) to the reference ue
            
            
            %approx
            %dother=sqrt(x^2+y^2); %distance between interfering bs(i,j) and reference bs(0,0) approximates the distance between bs(i,j) and the the reference ue
            % traditional model SIR= sectors/6 * (sqrt(3*N))^gamma
            
            % no downlink power control
            ptx=prx*(R^(gamma));
            
            % interfering power
            pti=ptx/(dother^(gamma)); % power received at reference mobile
            
            % take into account sectoring with different sub-channels
            angleBSiUEr=atan2d(x-ue_xr,y-ue_yr); % angle between the interfering bs(i,j) and the reference ue
            if (angleBSiUEr<0)
                angleBSiUEr=angleBSiUEr+360;
            end
            
            bs_si = ceil(angleBSiUEr/BeamWidth);% number of the sector of the interfering bs(i,j)exposed to the reference ue;
            if (bs_si~=ue_sr && subchannelingPerSector)
                % different sector, no interference if sub-channeling per
                % sector
                %disp('different sector, no interference');
                pti=0;
            end
            
            % taking into account cell load
            cl=rand(1,1);
            if cl>cellLoad
                %disp('no lad, no interference');
                pti=0;
            end
            
            % taking into account channel reuse
            if (F(i,j)~=F(bs_uir,bs_vir))
                % different channels, no interference
                %disp('different channels, no interference');
                pti=0;
            end
            
            ptot(t)=ptot(t)+pti;
        end
    end
end

% average S/Iother, reference ue reaches reference bs with prx
%no power control downlink
ptx_bsr=prx*R^(gamma);
prx_bsr=ptx_bsr/Rue^(gamma);

% average S/Iother
SIRdown = prx_bsr/mean(ptot);
SIdowndb=10*log10(SIRdown);
fprintf('uplink S/I %f dB, downlink S/I %f dB\n',SIupdb,SIdowndb);


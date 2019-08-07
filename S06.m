function S06(summer,current)
%Modified Scherrer blocking index
%summer=1--> summer season, summer=0--> winter season
%current=1-->current climate, current=0--> future time
load('latlon.mat')
year_start=1;
nstep=ceil(15/(latNH(2)-latNH(1)));
%latitudes and longitudes for calculating std
stdNH1=43;
stdNH2=65;
stdSH1=32;
stdSH2=54;

%define the regions latitudes and longitude 
SH_last_index=96;
NH_first_index=97;            
latSH=lat(1:SH_last_index);
latNH=lat(NH_first_index:end);
%North Pacific lon
lonNP_west=101;
lonNP_east=185;
lonNP=lon(lonNP_west:lonNP_east,1);
%North Antantic lon
lonNA_west=241;
lonNA_east=25;
lonNA=[lon(lonNA_west:end,1);lon(1:lonNA_east,1)];
%Russia
lonR_west=26;
lonR_east=100;
lonR=lon(lonR_west:lonR_east,1);
%Southern Hemisphere
lonSH_west=1;
lonSH_east=288;
lonSH=lon(lonSH_west:lonSH_east,1);
%load std
if current==1 && summer==1
    load('/current/25_SDTsummer.mat')
    stdd_Summer_NP=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_NP(:,:,stdNH1:stdNH2,:).*stdd_Summer_NP(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Summer_NA=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_NA(:,:,stdNH1:stdNH2,:).*stdd_Summer_NA(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Summer_R=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_R(:,:,stdNH1:stdNH2,:).*stdd_Summer_R(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Summer_SH=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_SH(:,:,stdSH1:stdSH2,:).*stdd_Summer_SH(:,:,stdSH1:stdSH2,:)),1))),3)),1);
    
    for j=1:length(latSH(stdSH1:stdSH2))
        stdd_Summer_SH(j)=stdd_Summer_SH(j).*abs(sind(45))./(abs(sind(latSH(stdSH1-1+j))));
    end
    for j=1:length(latNH(stdNH1:stdNH2))
        stdd_Summer_NP(j)=stdd_Summer_NP(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Summer_NA(j)=stdd_Summer_NA(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Summer_R(j)=stdd_Summer_R(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
    end
    
    stdd_NP=max(stdd_Summer_NP);
    stdd_NA=max(stdd_Summer_NA);
    stdd_R=max(stdd_Summer_R);
    stdd_SH=max(stdd_Summer_SH);
    clear stdd_Summer_SH stdd_Summer_R stdd_Summer_NA stdd_Summer_NP
elseif current==1 && summer==0
    load('/current/25_SDTwinter.mat')
    stdd_Winter_NP=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_NP(:,:,stdNH1:stdNH2,:).*stdd_Winter_NP(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Winter_NA=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_NA(:,:,stdNH1:stdNH2,:).*stdd_Winter_NA(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Winter_R=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_R(:,:,stdNH1:stdNH2,:).*stdd_Winter_R(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Winter_SH=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_SH(:,:,stdSH1:stdSH2,:).*stdd_Winter_SH(:,:,stdSH1:stdSH2,:)),1))),3)),1);
    for j=1:length(latSH(stdSH1:stdSH2))
        stdd_Winter_SH(j)=stdd_Winter_SH(j).*abs(sind(45))./(abs(sind(latSH(stdSH1-1+j))));
    end
    for j=1:length(latNH(stdNH1:stdNH2))
        stdd_Winter_NP(j)=stdd_Winter_NP(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Winter_NA(j)=stdd_Winter_NA(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Winter_R(j)=stdd_Winter_R(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
    end
    
    stdd_NP=max(stdd_Winter_NP);
    stdd_NA=max(stdd_Winter_NA);
    stdd_R=max(stdd_Winter_R);
    stdd_SH=max(stdd_Winter_SH);
    clear stdd_Winter_SH stdd_Winter_R stdd_Winter_NP stdd_Winter_NA
elseif current==0 && summer==1
    load('/future/25_SDTsummer.mat')
    stdd_Summer_NP=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_NP(:,:,stdNH1:stdNH2,:).*stdd_Summer_NP(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Summer_NA=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_NA(:,:,stdNH1:stdNH2,:).*stdd_Summer_NA(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Summer_R=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_R(:,:,stdNH1:stdNH2,:).*stdd_Summer_R(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Summer_SH=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Summer_SH(:,:,stdSH1:stdSH2,:).*stdd_Summer_SH(:,:,stdSH1:stdSH2,:)),1))),3)),1);
    for j=1:length(latSH(stdSH1:stdSH2))
        stdd_Summer_SH(j)=stdd_Summer_SH(j).*abs(sind(45))./(abs(sind(latSH(stdSH1-1+j))));
    end
    for j=1:length(latNH(stdNH1:stdNH2))
        stdd_Summer_NP(j)=stdd_Summer_NP(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Summer_NA(j)=stdd_Summer_NA(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Summer_R(j)=stdd_Summer_R(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
    end
    
    stdd_NP=max(stdd_Summer_NP);
    stdd_NA=max(stdd_Summer_NA);
    stdd_R=max(stdd_Summer_R);
    stdd_SH=max(stdd_Summer_SH);
    clear stdd_Summer_SH stdd_Summer_R stdd_Summer_NA stdd_Summer_NP
elseif current==0 && summer==0
    load('/future/25_SDTwinter.mat')
    stdd_Winter_NP=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_NP(:,:,stdNH1:stdNH2,:).*stdd_Winter_NP(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Winter_NA=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_NA(:,:,stdNH1:stdNH2,:).*stdd_Winter_NA(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Winter_R=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_R(:,:,stdNH1:stdNH2,:).*stdd_Winter_R(:,:,stdNH1:stdNH2,:)),1))),3)),1);
    stdd_Winter_SH=mean(squeeze(mean(sqrt(squeeze(mean((stdd_Winter_SH(:,:,stdSH1:stdSH2,:).*stdd_Winter_SH(:,:,stdSH1:stdSH2,:)),1))),3)),1);
    for j=1:length(latSH(stdSH1:stdSH2))
        stdd_Winter_SH(j)=stdd_Winter_SH(j).*abs(sind(45))./(abs(sind(latSH(stdSH1-1+j))));
    end
    for j=1:length(latNH(stdNH1:stdNH2))
        stdd_Winter_NP(j)=stdd_Winter_NP(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Winter_NA(j)=stdd_Winter_NA(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
        stdd_Winter_R(j)=stdd_Winter_R(j).*abs(sind(45))./(abs(sind(latNH(stdNH1-1+j))));
    end
    stdd_NP=max(stdd_Winter_NP);
    stdd_NA=max(stdd_Winter_NA);
    stdd_R=max(stdd_Winter_R);
    stdd_SH=max(stdd_Winter_SH);
    clear stdd_Winter_SH stdd_Winter_R stdd_Winter_NP stdd_Winter_NA
end



for m=1:40 %for loop over all the members in large ensembles, 40 for LENS and 20 for GFDL
    posAll={};
    %load anomalies 
    if current==1 && summer==0
        filename = sprintf('%s%d.mat','25_Winter_anomalies',m);
        load(['/current/' filename],'noseasonzSH','noseasonzNH')
        filename = sprintf('%s%d.mat','25_Winter',m);
        load(['/current/' filename])
 allsectors={ZwinterNP(:,:,:,16:end-15),ZwinterNA(:,:,:,16:end-15),ZwinterR(:,:,:,16:end-15),ZwinterSH(:,:,:,16:end-15)};    
    elseif current==1 && summer==1
        filename = sprintf('%s%d.mat','25_Summer_anomalies',m);
        load(['/current/' filename],'noseasonzSH','noseasonzNH')
        filename = sprintf('%s%d.mat','25_Summer',m);
        load(['/current/' filename])
 allsectors={ZsummerNP(:,:,:,16:end-15),ZsummerNA(:,:,:,16:end-15),ZsummerR(:,:,:,16:end-15),ZsummerSH(:,:,:,16:end-15)};    
    elseif current==0 && summer==0
        filename = sprintf('%s%d.mat','25_Winter_anomalies',m);
        load(['/future/' filename],'noseasonzSH','noseasonzNH')
        filename = sprintf('%s%d.mat','25_Winter',m);
        load(['/future/' filename])
 allsectors={ZwinterNP(:,:,:,16:end-15),ZwinterNA(:,:,:,16:end-15),ZwinterR(:,:,:,16:end-15),ZwinterSH(:,:,:,16:end-15)};    
    elseif current==0 && summer==1
        filename = sprintf('%s%d.mat','25_Summer_anomalies',m);
        load(['/future/' filename],'noseasonzSH','noseasonzNH')
        filename = sprintf('%s%d.mat','25_Summer',m);
        load(['/future/' filename])
        
 allsectors={ZsummerNP(:,:,:,16:end-15),ZsummerNA(:,:,:,16:end-15),ZsummerR(:,:,:,16:end-15),ZsummerSH(:,:,:,16:end-15)};    
    end
        for j=1:length(latSH)
            noseasonzSH(:,:,j,:)=noseasonzSH(:,:,j,:).*abs(sind(45))./(abs(sind(latSH(j))));
        end
        for j=1:length(latNH)
            noseasonzNH(:,:,j,:)=noseasonzNH(:,:,j,:).*abs(sind(45))./(abs(sind(latNH(j))));
end
    
    clear ZsummerNP ZsummerNA ZsummerR ZsummerSH
    xx={lonNP,lonNA,lonR,lonSH};
    yy={latNH,latNH,latNH,latSH};
    lonPlotNA=[lon(146:end);lon(1:145)+360];
    XXplot={lon,lonPlotNA,lon,lonSH};
    noseasonzNAplot=cat(2,noseasonzNH(:,146:end,:,:),noseasonzNH(:,1:145,:,:));
    DataPlot={noseasonzNH./stdd_NP,noseasonzNAplot./stdd_NA,noseasonzNH./stdd_R,noseasonzSH./stdd_SH};
    for sector=1:4
        Xplot=XXplot{sector};
        normaldata1=DataPlot{sector};
        pos={};
        Z=allsectors{sector};
        AR=zeros(size(Z));
        AR1=zeros(size(Z));
        AR2=zeros(size(Z));
        x=xx{sector};
        if sector==2
            x=[x(1:48);x(49:end)+360];
        end
        yyy=yy{sector};
        [X,Y]=meshgrid(x,yyy);
        X=X';Y=Y';
        [YY,nx,ny,day]=size(Z);
        if sector==4
            for J=1:size(Z,3)
                if latSH(J)<-75 || latSH(J)>-32
                    continue
                else
                    deltaD=15;
                    AR(:,:,J,:)= (Z(:,:,J-nstep,:)- Z(:,:,J,:))/deltaD;
                    AR1(:,:,J,:)= (Z(:,:,J,:)- Z(:,:,J+nstep,:))/deltaD;
                    AR2(:,:,J,:)= (Z(:,:,J+nstep,:)- Z(:,:,J+2*nstep,:))/deltaD;
                end
            end
        else
            for J=1:size(Z,3)
                if latNH(J)>75 || latNH(J)<32
                    continue
                else
                    deltaD=15;
                    AR(:,:,J,:)= (Z(:,:,J+nstep,:)- Z(:,:,J,:))/deltaD;
                    AR1(:,:,J,:)= (Z(:,:,J,:)- Z(:,:,J-nstep,:))/deltaD;
                    AR2(:,:,J,:)= (Z(:,:,J-nstep,:)- Z(:,:,J-2*nstep,:))/deltaD;
                end
            end
        end
        for y=year_start:YY
            B=zeros((nx*(ny+1)),2,'single');
            B(:,1)=(1:nx*(ny+1));
            B1=B;pp=B;
            Bold=B;Bold1=Bold;Bold2=Bold;
            Bold11=zeros(nx*(ny+1),5,'single');Bold22=Bold11;
            jj=1;ii=jj;
            Posblockinfo=[];
            cc=1;
            for d=1:day
                [d y sector m]
                BB  =squeeze(AR(y,:,:,d));
                BB1  =squeeze(AR1(y,:,:,d));
                BB2  =squeeze(AR2(y,:,:,d));
                index1=find(BB<-10);
                index2=find(BB1>0);
                index3=find(BB2<-5);
                index4=intersect(index1,index2);
                positive=intersect(index4,index3);
                if isempty (positive)
                    disp('P is empty')
                    Bold=zeros(nx*(ny+1),2);
                else
                    [pp,aii,bii]=intersect(Bold(:,1),positive);
                    lr=length(pp);
                    B=[pp,Bold(aii,end)+0.2*ones(lr,1)];
                    pp1=setdiff(positive,Bold(:,1));
                    lr=length(pp1);
                    B1=[pp1,0.2*ones(lr,1)];
                    Bp=[B;B1];
                    Bold=Bp;
                    if ~isempty(find(Bold(:,2)>=1))
                        blockp1=find(Bold(:,2)>=1);
                        lrp=length(blockp1);
                        lrp2=lrp;
                        lrp=ii+lrp-1;
                        bbp=Bold(blockp1,:);
                        [kk11,ai,bi]=intersect(Bold11(:,1),bbp(:,1));
                        ff=0;
                        if ~isempty(kk11)
                            ff=length(ai);
                            for ip=1:ff
                                
                                index1=find(Bold11(:,1)==bbp(bi(ip),1));
                                d1=Bold11(index1(end),end-2);
                                y1=Bold11(index1(end),end-1);
                                m1=Bold11(index1(end),end);
                                if m1==m
                                    if y1==y
                                        if  (d1==d-1)
                                            Bold11(index1(end),:)=[bbp(bi(ip),1),bbp(bi(ip),2),d,y,m];
                                            lrp=lrp-1;
                                        else
                                            Bold11(ii,:)=[bbp(bi(ip),1),bbp(bi(ip),2),d,y,m];
                                            ii=ii+1;
                                        end
                                    end
                                end
                                
                            end
                        end
                        if ~isempty(setdiff(bbp(:,1),Bold11(:,1)))
                            [pp2,aii]=setdiff(bbp(:,1),Bold11(:,1));
                            rr=length(pp2);
                            Bold11(ii:lrp,:)=[bbp(aii,1:2),d*ones(rr,1),y*ones(rr,1),m*ones(rr,1)];
                            ii=ii+lrp2-ff;
                        end
                    end
                end
            end
            Bold11(Bold11==0) = [];
            Bold11=reshape(Bold11,[],5);
            Bold11=sortrows(Bold11,5);
            if isempty (find(Bold11))
                disp('There is no positive block in the system')
            else
                Bold11(:,2)=int8(5*Bold11(:,2));
                [lf,rf]=size(Bold11);%ii=1;
                rr=int8(max(Bold11(:,2)));
                Posblockinfo=zeros(lf,rr+6);
                Posblockinfo(:,1:5)=Bold11;
                i=1;
                 %capturing the contour around the blocking grid point
                while i<=lf
                    gg=Posblockinfo(i,3)-(Posblockinfo(i,2));
                    for jj=1:(Posblockinfo(i,2))
                        c=contourc(double(Xplot),double(yyy),double(squeeze(normaldata1(y,:,:,gg+jj)))',[1 1]);%hh=colorbar;title(['Day' num2str(gg+jj)]);set(hh, 'ylim', [-2 2])
                        if ~isempty(c)
                            [xa,ya,z]=C2xyz(c);
                            for n =1:length(z)
                                    if (inpolygon(X(Posblockinfo(i,1)),Y(Posblockinfo(i,1)),xa{n},ya{n}))
                                         plot((xa{n}),(ya{n}),'r:','linewidth',2);
                                                Posblockinfo(i,jj+5)=polyarea(xa{n},ya{n})*cosd(yC);
                                                [ geom, iner, cpmo ] = polygeom(xa{n},ya{n});
                                                yC = geom(3);
                                    end
                                
                            end
                        else
                            Posblockinfo(i,:)=[];
                            lf=lf-1;
                            i=i-1;
                            break
                        end
                    end
                    i=i+1;
                end
            end
            if isempty (find(Posblockinfo))
                disp('There is no positive block in the system')
            else
                
                [lf,rf]=size(Posblockinfo(:,1));jj=1;
                if lf==1
                    Posblockinfo(1,:)=Posblockinfo(1,:);
                else
                    for uu=1:lf
                        if uu<=lf
                            ia=ismember(Posblockinfo(uu+1:lf,6:10),Posblockinfo(uu,6:10),'rows');
                            ia=find(ia);
                            if ~isempty(ia)
                                as=max(find((Posblockinfo(ia+uu,2))==max(Posblockinfo(ia+uu,2))));
                                Posblockinfo(uu,:)=Posblockinfo(ia(as)+uu,:);
                                Posblockinfo(ia+uu, :) = [];
                                lf=lf-(length(ia));
                            else
                                Posblockinfo(uu,:)=Posblockinfo(uu,:);
                            end
                            
                        else
                            if lf==1
                                Posblockinfo(1,:)=Posblockinfo(1,:);
                            elseif Posblockinfo(lf-1,2:end)==Posblockinfo(lf,2:end)
                                Posblockinfo(lf, :) = [];
                                lf=lf-1;
                            elseif Posblockinfo(lf-1,6:10)==Posblockinfo(lf,6:10)
                                Posblockinfo(lf-1,:)=Posblockinfo(lf,:);
                                Posblockinfo(lf, :) = [];
                                lf=lf-1;
                            else
                                Posblockinfo(lf,:)=Posblockinfo(lf,:);
                            end
                        end
                    end
                end
            end
            if isempty (find(Posblockinfo))
                disp('There is no positive block in the system')
            else
                [lf,rf]=size(Posblockinfo(:,1));ii=1;
                for i=1:lf
                    lf1=Posblockinfo(ii,2);
                    aa=Posblockinfo(ii,6:6+lf1-1);
                    aq=find(~aa);
                    if ~isempty(aq)
                        Posblockinfo(ii,:)=[];
                        ii=ii;
                        lf=lf-1;
                    else
                        ii=ii+1;
                    end
                end
            end
            if isempty (find(Posblockinfo))
                disp('There is no positive block in the system')
            else
                [lf,rf]=size(Posblockinfo(:,1));jj=1;
                if lf==1
                    Posblockinfo(1,:)=Posblockinfo(1,:);
                else
                    for uu=1:lf
                        if uu<=lf
                            [ia,ib]=ismember(Posblockinfo(uu+1:lf,6:end-1),Posblockinfo(uu,6:10));
                            [a1,a2]=find(ib);
                            ia=unique(a1);
                            if ~isempty(ia)
                                if max(Posblockinfo(ia+uu,2))>=Posblockinfo(uu,2)
                                    as=max(find((Posblockinfo(ia+uu,2))==max(Posblockinfo(ia+uu,2))));
                                    Posblockinfo(uu,:)=Posblockinfo(ia(as)+uu,:);
                                    Posblockinfo(ia+uu, :) = [];
                                    lf=lf-(length(ia));
                                elseif max(Posblockinfo(ia+uu,2))<=Posblockinfo(uu,2)
                                    Posblockinfo(ia+uu,:)=[];
                                    lf=lf-(length(ia));
                                end
                            else
                                Posblockinfo(uu,:)=Posblockinfo(uu,:);
                            end
                            
                        else
                            if lf==1
                                Posblockinfo(1,:)=Posblockinfo(1,:);
                            elseif Posblockinfo(lf-1,2:end)==Posblockinfo(lf,2:end)
                                Posblockinfo(lf, :) = [];
                                lf=lf-1;
                            elseif Posblockinfo(lf-1,6:10)==Posblockinfo(lf,6:10)
                                Posblockinfo(lf-1,:)=Posblockinfo(lf,:);
                                Posblockinfo(lf, :) = [];
                                lf=lf-1;
                            else
                                Posblockinfo(lf,:)=Posblockinfo(lf,:);
                            end
                        end
                    end
                end
            end
            if isempty (find(Posblockinfo))
                disp('There is no positive block in the system')
                pos=[pos,nan];
                
            else
                Posblockinfo=sortrows(Posblockinfo,5);
                [lf,rf]=size(Posblockinfo(:,1));
                for i=1:lf; rr=find(Posblockinfo(i,:));rr=rr(end);Posblockinfo(i,end)=mean(Posblockinfo(i,6:rr));end
                areapos(m,sector,y)=mean(Posblockinfo(:,end));
                pos=[pos,Posblockinfo];
            end
        end
        posAll{1,sector}=pos;
    end
    %save the results
    if summer==0 && current==1
        filename = sprintf('%s%d.mat','scherrer_winterC',m);
    elseif  summer==0 && current==0
        filename = sprintf('%s%d.mat','scherrer_winterF',m);
    elseif  summer==1 && current==0
        filename = sprintf('%s%d.mat','scherrer_summerF',m);
    elseif  summer==1 && current==1
        filename = sprintf('%s%d.mat','scherrer_summerC',m);
    end
    save(filename,'posAll','areapos','-v7.3')
    clear posAll areapos
end

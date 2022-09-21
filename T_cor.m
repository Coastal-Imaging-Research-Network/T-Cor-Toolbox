% T-cor : code for a temporal method for celerity detection and bathymetry
% inversion following Almar et al., 2008. Including the use of multiple
% frequencies and error estimates (Bergsma and Almar, 2018;
% Thuan et al, 2019; Abessolo et al., 2019; Abessolo et al, 2020).

close all
clear all
clc

addpath('T_cor_functions\')
f= 'Data_test\'; % Images directory

disp('----------------------------------------------------------------')
disp('TEMPORAL METHOD FOR CELERITY DETECTION AND BATHYMETRY INVERSION')
disp('----------------------------------------------------------------')
disp(' ')

disp('Define stack coordinates')
sk1=[130 520]; % Camera 1 of the video system at Grand Popo, Benin
sk2=[360 1600];
disp('          ... Ok')
disp(' ')

disp('Compute pixel footprint')
[d]=T_cor_funct1(sk1,sk2,[f,'\20140316\I_201403161200.jpg'],'functions_rect\RectGPP.mat');
d=[0 d]; % Pixel footprint
dx=smooth(d,50); % Smooth pixel footprint
DX=cumsum(dx); % Cross-shore distance
figure('Color','w')
hold on
hp(1)=plot(1:length(d),d,'k','LineWidth',1,'MarkerSize',1);
hp(2)=plot(1:length(d),dx,'r','LineWidth',2);
xlabel('Pixels','FontSize',14)
ylabel('Pixel footprint (m)','FontSize',14)
legend(hp,'Initial','Smoothed')
grid on
hold off
box on
disp('          ... Ok')
disp(' ')

disp('Looking for the start date of the processing')
lsdat=ls('Param*');
if isempty(lsdat)
    t_start=datenum(2013,03,1); % Start date
    save('Parameters.mat','t_start')
else
    load('Parameters.mat');
end
t_end = datenum(2016,09,30); % End date
disp('          ... Ok')
disp(' ')

disp('Processing timestack images')
for t_p=t_start:t_end
    
    % Initialization
    files=[];
    Rj=[];
    
    % Image directories to be processed
    lsrep=ls([f,'20*']);
    ind=find(datenum(lsrep,'yyyymmdd')>=t_p&datenum(lsrep,'yyyymmdd')<t_p+1);
    
    % Day by day processing
    if ~isempty(ind==0)
        fprintf('     Day being processed %s \n',datestr(t_p))
        lsrep=lsrep(ind,1:8);
        lsjour=ls([f,lsrep,'/S_3_*']); % Identification des images (timestack) du jour considéré
        
        % Vecteurs de sauvegarde journaliers
        t1=[];
        h1=[];
        T1=[];
        C1=[];
        WL1=[];
        WL2=[];
        compt1=0;
        
        % Process image by image
        for i=1:size(lsjour,1)
            
            % Step 1: Load image
            files=[files [f,lsrep,'/',lsjour(i,:)]'];
            t=datenum(files(end-15:end-4,end)','yyyymmddHHMM'); % Processed image time
            fprintf('          Image in process %s\n',datestr(datenum(lsjour(i,5:16),'yyyymmddHHMM'),'mmmm dd, yyyy HH:MM AM'))
            S=imread(files(:,i)');
            % image(S)
            
            % Step 2: Rectify image
            [S0,resc]=T_cor_funct2(S,dx); % S0 = rectified image, resc = spatial resolution in meters
            [nt,nc,ncol]=size(S0);
            
            % Step 3: Formatting the image
            So=double(S0);
            S_std=nanstd(So);
            clear S
            
            % Step 4: Delimitation of the surf zone
            S_std0=smooth(S_std(:,:,3)); % We consider only the blue component
            clear iS
            iS=find(S_std0>0.5*max(S_std0)); % Indices corresponding to the surf zone
            
            % Displaying the image to be processed
            figure(102)
            clf
            set(gcf,'Color','w')
            image(S0)
            title(sprintf('%s',datestr(datenum(lsjour(i,5:16),'yyyymmddHHMM'),'mmmm dd, yyyy HH:MM AM')))
            box on
            hold on
            plot(smooth(S_std(:,:,1),30).*40,'r') % Blue component
            plot(smooth(S_std(:,:,2),30).*40,'g') % Red component
            plot(smooth(S_std(:,:,3),30).*40,'b') % Green component
            xlim([1 size(S0,2)])
            patch([min(iS) min(iS)],[0 size(S0,1)],'k','Linewidth',2)
            patch([max(iS) max(iS)],[0 size(S0,1)],'k','Linewidth',2)
            hold off
            box on
            xlabel('X (m)','FontSize',14)
            ylabel('Samples','FontSize',14)
            clear S0
            
            % Step 5: Criteria for selecting a good image
            if abs(max(iS)-min(iS))>10&abs(max(iS)-min(iS))<200 % This parameter is site dependent
                disp('               Good image. Continue processing... !')

                S0=double(So(:,:,3)); % We consider only the blue component
                Nlim=400; % Number of samples to consider (to save time) over the 15 min of the image, knowing that 15 min = 1800 samples
                S1=S0(1:min([Nlim size(S0,1)]),:);
                clear S0
                
                % Step 6: Pre-processing of the image (elimination of spatial redundancy and noise, to keep only the motion of the waves)
                icmin=1;
                icmax=nc;
                dt=0.5; % dt=1/freq with freq=2 (acquisition frequency of the camera: 2 Hz)
                % methodPreT=1; % Méthode de prétraitement de l'image.
                % De ce que j'observe, "methodPreT=1" est la meilleure et c'est celle choisie par Almar
                [S2]=T_cor_funct3(S1,icmin,icmax,dt,resc,1);
                clear S
                
                % Step 7: Temporal cross-sorrelation
                dpha=3; % We set the time lag "dpha" which must be less than the wave period
                dc=200; % Max spacing for cross-correlation
                res=10; % Resolution for the cross-correlation
                wdf=0.7; % bandpass filter width (in ratio)
                lx=ones(1,size(S2,2)).*resc; % Cross-shore resolution of the image (in meters)
                
                try
%                     [h,T,R,C,W]=T_cor_funct4(S2,dpha,dt,lx,dc,res,wdf);
                    [h,T,R,C,W]=T_cor_funct4_bis(S2,dpha,dt,lx,dc,res,wdf);
                catch
                    h=[];T=[];R=[];C=[];W=[];
                end
                
%                 % Visualization of the cross-correlation matrix
%                 a=1:size(R,1);
%                 b=1:size(R,2);
%                 [A,B]=meshgrid(a,b);
%                 figure(102)
%                 clf
%                 set(gcf,'Color','w')
%                 contourf(A,B,R');
%                 colorbar;
%                 title(sprintf('Intercorrelation matrix %s',datestr(datenum(lsjour(i,5:16),'yyyymmddHHMM'),'mmmm dd, yyyy HH:MM AM')))
                
                % We add the successive intercorrelations matrices (taking care not to consider the NaN values)
                if ~isempty(Rj)
                    clear in1 in2 in
                    in1=find(~isnan(Rj));
                    in2=find(~isnan(R));
                    in=unique([in1;in2]);
                    for ik=1:length(in)
                        clear i1 i2
                        i1=find(in1==in(ik));
                        i2=find(in2==in(ik));
                        if ~isempty(i1)&~isempty(i2)
                            Rj(in1(i1))=Rj(in1(i1))+R(in2(i2));
                            compt1(in1(i1))=compt1(in1(i1))+1;
                        elseif isempty(i1)&~isempty(i2)
                            Rj(in(ik))=R(in2(i2));
                            compt1(in(ik))=compt1(in(ik))+1;
                        end
                    end
                else
                    Rj=R;
                    compt1=zeros(size(R,1),size(R,2));
                    clear in
                    in=find(~isnan(R));
                    compt1(in)=1;
                end
                
                clear R
                
                % Step 8: Save data
                h1=[h1 h];
                C1=[C1 C];
                WL1=[WL1 W];
                t1=[t1 t];
                T1=[T1 T];
                disp('          ... Ok')
                disp(' ')

            else
                disp('               Process aborted : Blurred image - Image floue !')
            end
            
        end
        
        if ~isempty(T1)
            R2j=NaN(size(Rj,1),size(Rj,2));
            clear iz
            iz=find(compt1~=0);
            R2j(iz)=Rj(iz)./compt1(iz);
            
%             % Visualization of the cross-correlation matrix
%             a=1:size(R2j,1);
%             b=1:size(R2j,2);
%             [A,B]=meshgrid(a,b);
%             figure(103)
%             clf
%             set(gcf,'Color','w')
%             contourf(A,B,R2j');
%             colorbar;
%             title(sprintf('Daily-averaged intercorrelation matrix %s',datestr(datenum(lsjour(i,5:16),'yyyymmddHHMM'),'mmmm dd, yyyy')))
            
            savdir='Output\'; % Output directory
            File_name=strcat('bathy_',datestr(t_p,'yyyymmdd'));
            save([savdir,File_name],'h1','t1','T1','R2j','C1','WL1');
        end
        
        fprintf('     End of treatment of %s OK\n',datestr(t_p))
        t_start = t_p + 1;
        save('Parameters.mat','t_start')
        disp('  ')
        
    end
    
end
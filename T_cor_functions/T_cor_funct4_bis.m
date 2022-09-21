% Cette fonction calcule la corrélation entre une série à l'indice "ic" et
% ses positions précédentes pour estimer la célérité

function [Depthf,Tp,R2M1,C,Wavelengthf]=T_cor_funct4_bis(S,dpha,dt,dx,dc,res,wdf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%S= matrice d'entrée (matrice spatio-temporelle) [ntim ncol]
%Waves should come from up left (1,1) to bottom right (nt,nc)
%dpha : décalage temporel fixé
%dt: résolution temporelle
%dx: résolution spatiale
%dc=Ecartement max pour la cross-correlation
%res: resolution pour la cross-correlation
%wdf=largeur du filtre passe bande (en ratio). Par défaut, wdf=0.7
%Outputs
% Depthf: profondeur
% Tp: periode pic
% R2M1 : matrice d'intercorrelation
% C : Célérité
% Wavelengthf : Longueur d'onde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A2=double(A2(:,:,1));
% A2=detrend(double(abs(diff(A2))));

S=S(end:-1:1,:);

S=detrend(S);

if length(dx)==1
    dx=dx.*ones(size(S,2),1);
end
% [A2]=ImagePreProcessing_BI(A2,1,size(A2,2),1,1,1);

%Spectre
[Pw,ffw]=pwelch(S(:,100),[],[],[],1/dt);
% [Pw,ffw]=pwelch(S,[],[],[],1/dt);
ipicf=find(slidefun(@nanmax,15,Pw)./Pw==1&Pw>nanmean(Pw));%Recherche des pics
% figure;plot(Pw);
% hold on;plot(ipicf,Pw(ipicf),'sr');hold off;
[nt nc cc]=size(S);

Depth=NaN.*ones(1,nc-dc);
Wavelength=NaN.*ones(1,nc-dc);
Celerity=NaN.*ones(1,nc-dc);

% save('C:\Docs_IRD\Codes\SERR1D_DepthInversion/toto')

for ifreq=1:length(ipicf)
    
    % S2=S;
    A2=smooth2(S-smooth2(S,round(1/(dt*ffw(ipicf(ifreq))))+round(wdf/(dt*ffw(ipicf(ifreq)))),1),round(1/(dt*ffw(ipicf(ifreq))))-round(wdf/(dt*ffw(ipicf(ifreq)))),1);
    A2=single(A2(round(4.*wdf/(dt*ffw(ipicf(ifreq)))):end-round(4.*wdf/(dt*ffw(ipicf(ifreq)))),:));
    [nt nc cc]=size(A2);
    
    % n-1 est la valeur du décalage temporel fixé (en pas de temps)
    n=floor(dpha/dt);
    R2(1:dc-1)=0;
    R1M=ones(nc-dc,dc-1,'double')*0;
    
    for ic=1:res:nc-dc
        for lc=1:dc-1
            try
                % variable en temps
                R6 = corrcoef(A2(1:(nt-(n+1)),ic+dc),A2(n+2:nt,ic+dc-lc+1));
                R7 = corrcoef(A2(1:(nt-(n)),ic+dc),A2(n+1:nt,ic+dc-lc+1));
                R8 = corrcoef(A2(1:(nt-(n-1)),ic+dc),A2(n:nt,ic+dc-lc+1));
                % variable en espace
                R9 = corrcoef(A2(1:(nt-n),ic+dc+1),A2(n+1:nt,ic+dc-lc+1));
                R10 = corrcoef(A2(1:(nt-n),ic+dc-1),A2(n+1:nt,ic+dc-lc+1));
                
                % R2(lc)=(dpha/dt-n).*R7(2,1)+(1-((dpha/dt)-n)).*R8(2,1);
                R2(lc)=double(nanmean([R6(2,1) R7(2,1) R8(2,1) R9(2,1) R10(2,1)]));
                
                
                % R2_0(lc)=R6(2,1);
                % R2_1(lc)=R7(2,1);
                % R2_2(lc)=R8(2,1);
            catch
                R2(lc)=double(nanmean([R6(2,1) R7(2,1) R8(2,1)]));
            end
        end % lc
        R1M(ic,1:dc-1)=smooth(R2,3);
    end % ic
    
    clear R6 R7 R8 R9 R10 % Ajout Greg
    
    for irt=1:length(R2)
        R(:,irt)=interp1(1:res:length(R1M(:,irt)),R1M(1:res:length(R1M(:,irt)),irt),1:1:nc-dc);
        R2M(:,irt)=double(smooth(R(:,irt),4));
    end
    
    clear R % Ajout Greg
    
    %Celerity
    R2M=smooth2(R2M,10,2);
    C2M=smooth2(R2M,20,20);
    
    if ifreq==1  % Début ajout Greg
        R2M1=R2M;
    end
    
%     %     Visualisation de la matrice d'intercorrelation
%     a=1:size(R2M,1);
%     b=1:size(R2M,2);
%     [A,B]=meshgrid(a,b);
%     figure;surf(A,B,R2M');view(2);shading('interp');colorbar;
        
    clear R2M % fin ajout Greg
    
    %Wavelength
    % for i=1:size(R2M,1)
    % [hs,htiers,hrms,trms,Hmax,h,Tp,sl,tcrete,t]=Wave_Char(R2M(i,:),1,0,2);
    % L2M(i)=trms ;
    % end
    % L2M=smooth2(L2M,10,1);
    
    %Wave period
    % clear TC2M
    % for i=1:length(1+dc/2:nc-dc/2)
    % [hs,htiers,hrms,trms,Hmax,h,Tp,sl,tcrete,t]=Wave_Char(A2(1:min([500 round(length(A2(:,1))/3)]),i),1,0,2);
    % T2M(i)=trms;
    % end
    %
    % T2M=smooth2(T2M,10,1).*dt;
    Tp(ifreq)=1./ffw(ipicf(ifreq));
    coef(ifreq)=Pw(ipicf(ifreq))/sum(Pw(ipicf)); % Contribution de chaque composante fréquentielle
    
    % NB : Cette interpolation est nécessaire car le calcul des
    % corrélations a été fait avec une résolution de res=10 pour alléger le
    % calcul. Il afut alors revenir à un pas réduit (pas=1).
    [X Y]=meshgrid(1:size(C2M,2),1:size(C2M,1));
    %     [XX YY]=meshgrid(1:0.1:size(C2M,2),1:0.1:size(C2M,1));
    [XX YY]=meshgrid(1:0.5:size(C2M,2),1:0.5:size(C2M,1));
    C2M2=interp2(X,Y,C2M,XX,YY);
    
    %     [ftg indR37]=max(smooth2(C2M2(:,1:round(dc.*10/2))',15,15));
    [ftg,indR37]=max(smooth2(C2M2(:,1:round(dc.*2/2))',15,15));
    indR37=round(indR37(1:2:end)/2);
    
    % Ajout Greg
    xpha=NaN(1,nc-dc);
    for i=1:1:nc-dc
        if ~isnan(indR37(i))
            xpha(i)=abs(dx(i+dc)-dx(i+dc-indR37(i)+1));
        else
            xpha(i)=NaN;
        end
    end
    
%     Depth(ifreq,round(dc/2+(1:length(indR37))))=LinearC(Tp(ifreq),dx(1+dc/2:nc-dc/2)'.*indR37./dpha);
%     Wavelength(ifreq,round(dc/2+(1:length(indR37))))=Tp(ifreq).*dx(1+dc/2:nc-dc/2)'.*indR37./dpha;
%     Celerity(ifreq,round(dc/2+(1:length(indR37))))=dx(1+dc/2:nc-dc/2)'.*indR37./dpha;

    Depth(ifreq,1:length(indR37))=LinearC(Tp(ifreq),xpha./dpha);
    Wavelength(ifreq,1:length(indR37))=Tp(ifreq).*xpha./dpha;
    Celerity(ifreq,1:length(indR37))=xpha./dpha;
end

% Coefficients dependant de la Cdeepwater
Co=1-Celerity'./(9.81*repmat(Tp,size(Celerity,2),1)/(2*pi));
Co(find(Co<0))=0;

% Co=Pw(ipicf)./sum(Pw(ipicf));%coefficient
% Co=repmat(Co',size(Celerity,2),1);
% Co(find(Co<0))=0;

% Depth=median(Depth);
% Co=Pw(ipicf)./sum(Pw(ipicf));%coefficient

Depth(isfinite(Depth)==0)=0;

Depthf=0;
for ifreq=1:length(ipicf)
    if length(ipicf)>1
        Depthf=(Co(:,ifreq)./sum(Co')').*Depth(ifreq,:)'+Depthf;
    else
        Depthf=(Co(:,ifreq)./Co).*Depth(ifreq,:)'+Depthf;
    end
end
C=0;
for ifreq=1:length(ipicf)
    if length(ipicf)>1
        C=(Co(:,ifreq)./sum(Co')').*Celerity(ifreq,:)'+C;
    else
        C=(Co(:,ifreq)./Co).*Celerity(ifreq,:)'+C;
    end
end
Wavelengthf=0;
for ifreq=1:length(ipicf)
    if length(ipicf)>1
        Wavelengthf=(Co(:,ifreq)./sum(Co')').*Wavelength(ifreq,:)'+Wavelengthf;
    else
        Wavelengthf=(Co(:,ifreq)./Co).*Wavelength(ifreq,:)'+Wavelengthf;
    end
end
% Tp=median(Tp); % Moyenne statique
Tp_0=0;
for ifreq=1:length(ipicf)
    Tp_0=Tp_0+(coef(ifreq)*Tp(ifreq)); % Période pondérée
end
clear Tp
Tp=Tp_0;
clear Tp_0
end


function [hs,htiers,hrms,trms,Hmax,h,Tp,sl,tcrete,t]=Wave_Char(d,dt,filt,meth)
%Author: Rafael Almar (rafael.almar@ird.fr)
%Calcul les caractéristiques des vagues à partir d'une série temporelle de
%hauteur d'eau (ou équivalent)
%Inputs
%d : vecteur de serie temporelle
%dt : pas de temps
%filt : filtrage (1) /pas de filtre (0)
%meth : Mean zero crossing
% meth=1 mean zero crossing
% meth=2 up zero crossing
% meth=3 down zero crossing
%z=temps des vagues
%Outputs
%h : hauteurs de vague individuelles (ex: pour faire un histogramme)
%htiers;
%hrms;
%trms;
%Hmax;
%h; hauteurs individuelles des vagues (vecteur)
%Tp; periode peak
%sl: slope (m/s) (vecteur)
%tcrete= temps des cretes individuelles  (vecteur)
%t:periode individuelle (vecteur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%INITIALISATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Chargement du fichier de données :
%d=load('D:\Documents\THESE\DONNEES_INSITU\DGO41001.txt');
% Numero de la colonne du capteur dans le fichier (ex:10)
%nc=5;
% Pas de temps (dt) des mesures des capteurs
%dt = 0.5; % secondes

%Nombre de pas de temps de l'échantillon (nt)
[n1 n2]=size(d);
if n1==1
    d=d';
end




sz=size(d);
nt=sz(1);
zc=meth;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Fin initialisation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chargement de la colonne du capteur à analyser :
% ici colonne 10
% durée de donnée à analyser : environ 10 min (soit 1820 mesures)
k=d(1:nt);
% On detrend le signal brut
%generation d'un vecteur pour les abscisses
%  x=[1:nt];
%  a=reshape(x,nt,1);
%
%  % generation des coefficients du polynome de degré un
%   p=polyfit(a,k,1);
%
%  %Suppression de la partie trendée de la courbe
%
%  k2=k-p(1).*a;
%
%  % Suppression de la partie moyenne
%  km=k2;
%  k2=km-mean(km);

k2=NaNdetrend(d);

clear k4
if(filt==1)
    freq = (1./dt)*(1:nt)/nt;
    %frequences de coupure : %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fréquence de coupure vers les basses fréqences
    % fcb= 1/Tlimit
    fcb=0.05;
    % fréquence de coupure vers les hautes fréquences
    % fch=1/Tlimit
    fch=0.17;
    f1 = fft(k2);
    Indb = find(freq < fcb);
    Indh = find(freq > fch);
    f1(Indb) = 0.0;
    Sinv = ifft(f1);
    k3 = Sinv.*conj(Sinv);
    f1(Indh) = 0.0;
    Sinv = ifft(f1);
    k4 = real(Sinv);
    k2=k4;
    k3=FiltreMean(k2,round(dt));
else%filt==0
    
    
    
    %      k3=k2'-FiltreMean(k2',round(10/(2*dt)));
    k3=FiltreMean(k2,round(dt));
    k4=k2;
    
    
    
end %filt
%  size(k2)
%  figure;plot(k3)
%  round(5/(2*dt))
%k4 est un signal filtré pour les hautes et basses fréquences

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Mesures de la physique des vagues%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Localisation des zero du signal :

%Mean zero crossing
% zc=1 mean zero crossing
% zc=2 up zero crossing
% zc=3 down zero crossing
% for zc=1:2

clear z; j=0;
for i=2:nt-1
    if(k3(i)*k3(i+1)<=0 && k3(i-1)*k3(i)>=0 )
        j=j+1;
        z(j)=i;
    end
end

%  figure(36);plot(k2,'k');hold on;plot(k3,'r');plot(z,0,'r*')

%calcul de la taille du vecteur dans lequel sont stockés les indices des zeros
k=0;
sz=size(z);
for i=1:sz(2)
    if(z(i)~=0)
        k=k+1;
    end
end

if(zc==2)
    clear z2; j=0;
    for i=2:nt-1
        if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])<0 )
            j=j+1;
            z2(j)=i;
        end
    end
    %  figure(36);plot(k3,'k');hold on;plot(z2,0,'r*')
elseif(zc==3)
    clear z2; j=0;
    for i=2:nt-1
        if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])>0 )
            j=j+1;
            z2(j)=i;
        end
    end
    k2=0;
    sz=size(z2);
    for i=1:sz(2)
        if(z2(i)~=0)
            k2=k2+1;
        end
    end
end



% Recherche du max de la hauteur de la vague entre deux zeros

hmax(1:k-1)=0;
hmin(1:k-1)=0;

for i=1:k-1
    if(mean(k4(z(i):z(i+1)))>0.)
        if(max(k3(z(i):z(i+1)))>0)
            hmax(i)=max(k4(z(i):z(i+1)));
        else
            hmax(i)=0;
        end
        hmin(i)=0.;
        
    else
        hmax(i)=0.;
        if(min(k3(z(i):z(i+1)))<0)
            hmin(i)=min(k4(z(i):z(i+1)));
            
        else
            hmin(i)=0.;
        end
    end
end

% calcul de la hauteur de chaque vague passant par le capteur
h=[];
t=[];
kh=0;
%hauteurs: h
% h(1:k-2)=0;
%pentes:sl
sl=[];
for i=1:k-2
    
    if(hmin(i)~=0 & hmax(i+1)~=0)
        if(hmin(i+1)==0 & hmax(i)==0)
            kh=kh+1;
            h(kh)=hmax(i+1)-hmin(i);
            sl(kh)=hmax(i+1)-hmin(i)./((z(i+1)-z(i)));
            try
                tcrete(kh)=z(i+1)+0.1*abs((z(i)-z(i+1)));
                t(kh)=2.*(z(i+1)-z(i))*dt;
            end
        end
    end
end
% Calcul de la hauteur rms moyenne des vagues sur la durée de mesure
%hrms1=0;
clear hrms ord htiers

[sd ord]=sort(h(find(h~=0)),'descend');
%Nombre de vagues pour le calcul
% disp(['Nombre de vagues : ',num2str(length(find(h~=0)))])

%Hauteur max
Hmax=nanmax(h);

%Hauteur 1/3
htiers=nanmean(h(ord(1:round(length(ord)/3))));% on prend le premier tiers pour la hauteur rms
%Hauteur rms
hrms=nanmean(h(1:kh));

%Hauteur significative
hs=4.*std(k4);
% disp(['Hs=',num2str(hs),'; Hrms=',num2str(hrms),'; H1/3=',num2str(htiers),'; Hmax=',num2str(Hmax)])
%Calcul de la periode de la houle pour chaque vague

% clear t
% t(1:k-1)=0;
% for i=1:k-1
%     if(z(i)~=0 | z(i+1)~=0 )
%         t(i)=2*(z(i+1)-z(i))*dt;
%     end
% end

%Calcul de la periode rms de la houle et affichage

trms=nanmean(t);
if(zc~=1)
    trms=nanmean(diff(z2))*dt;
end

%%%%%Periode Peak%%%%%%%%%%%
y=NaNdetrend(d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nombre de points pour la fft
nx=length(y);
%Pas de temps
% x=(1:nx);
%Resolution temporelle (periode d'échantillonnage) (en s)
% dt=1;

%FFT en précisant la taille du vecteur
Pv = pwelch(NaNdetrend(d));
[valmx indmx]=max(Pv);
Tp=1./(indmx./(dt*length(Pv)));

%Affichage de la periode
% disp(['Periode = ',num2str(trms)])

% try
% %Positin des cretes des vagues
% iii=find(tt>0);
% tc=tt(iii)+0.5*diff([tt(iii) max(tt(iii))]);
% id=round(tc);
% tcrete=id(find(k4(id)>0));
% catch
%     tcrete=[];
% end

h=h(1:kh);
% end%zc
try
    length(tcrete)>0;
catch
    tcrete=NaN;
end

end


function [B]=ImagePreProcessing_BI(A,icmin,icmax,dt,resc,methodPreT)
%size(A)= n*m
[nt nc ncol]=size(A);
% %Temporal resolution
% rest=2;
%Spatial resolution
% resc=8;

% tic
% dt=0.5; %BISCA

clear B
%B=A(:,:,1);



if(methodPreT==0)
    B=double(A);
elseif(methodPreT==1)
    % A3(1:nt,icmin:resc:icmax )=zeroavg(double(A(1:nt,icmin:resc:icmax)));
    A3=A;
    % per=PeriodWaves(A3,icmax-50,icmax,dt,2);
    for ic=icmin:resc:icmax
        
        % A2=zeroavg(double(A3(1:nt,ic,1)));
        % interval=[1/(per+2+1) 1/(per+2-1)];
        interval=[1 0.04];
        % %Methode 1 (ancienne) spectrale : filtre ideal filter
        count1=timeseries(A3(:,ic),(1:nt)*dt);
        %Put extrema for unfiltered signal
        y1=idealfilter(count1,interval,'pass');
        y=double(y1.Data);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        B(:,ic)=y;%Norma(y,2);
    end
elseif(methodPreT==2)
    
    for ic=icmin:4*resc:icmax
        [hs,htiers,hrms,trms,Hmax,h]=Wave_Char(smooth(double(A(:,ic,1)),2),dt,0,1);
        T(ic)=trms;
    end
    To=min(T(find(T>3)));
    
    for ic=icmin:resc:icmax
        %Methode 2 : temporelle
        v=double(A(:,ic,1));
        % B(:,ic)=Norma(v'-slidefun (@mean,24,v'));
        B(:,ic)=smooth (v'-smooth(v',round(6*To)),round((2*To)/3));
        % [hs,htiers,hrms,trms,Hmax,h]=Wave_Char(B(:,ic),dt,0,1);trms
        
        % B(:,ic)=Norma(v);
        % B(:,ic)=FiltreMean(v'-FiltreMean(v',10),2);
    end
elseif(methodPreT==3)
    for ic=icmin:resc:icmax
        %Methode 2 : temporelle
        v=detrend(double(A(:,ic,1)));
        %B(:,ic)=(((v'-FiltreMean(v',3)))');
        B(:,ic)=smooth(smooth(v-smooth(v,60),10),3);
    end
    
    %   figure(196);pcolor(B(1:1200,:));shading interp
    
elseif(methodPreT==4)
    clear A2
    
    A2(1:nt,icmin:icmax)=zeroavg(double(A(1:nt,icmin:icmax)));
    per=PeriodWaves(A2(1:nt,icmin:icmax),10,20,dt,2);
    for ic=icmin:icmax
        % filtre 2))
        count1=timeseries(A2(:,ic),(1:nt)*dt);
        %Put extrema for unfiltered signal
        % interval=[1/(per+2+1) 1/(per+2-1)];
        interval=[0.12 0.14];
        
        y1=idealfilter(count1,interval,'pass');
        y=double(y1.Data);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        B(:,ic)=y;
        %  B(:,ic)=A2(:,ic);
        B(:,ic)=Norma(B(:,ic),2);
    end
elseif(methodPreT==5)
    for ic=icmin:resc:icmax
        %Methode 2 : temporelle
        v=double(A(:,ic,1));
        %B(:,ic)=(((v'-FiltreMean(v',3)))');
        
        B(:,ic)=smooth(v'-smooth(v',7),2);
        B(:,ic)=Norma(B(:,ic),2);
    end
    
    
else
    
    
    %Smoothing
    %Facteur horizontal
    facth=1;
    %Facteur vertical
    factv=1;
    
    B=double(A);
    
    B(factv+1:nt-factv-2,facth+1:nc-facth-2)=smoothc(double(A),factv,facth);
    
    
    
end

if  methodPreT~=4
    B(find(isnan(B)==1))=0;
    
    %Interpolation spatiale entre icmin et icmax
    for irt=1:nt
        B(irt,icmin:icmax)=interp1(icmin:resc:icmax,B(irt,icmin:resc:icmax),icmin:icmax);
    end
end



disp('Pre-Traitement ok')
%  figure(52);
%  clf
% subplot(4,1,1);pcolor(B(1:40,:));shading flat
% % subplot(4,1,2);pcolor(A3);shading flat
%  subplot(4,1,3);plot(double(A(:,1600)));set(gca,'ydir','reverse')
%  subplot(4,1,4);plot(B(:,1600));set(gca,'ydir','reverse')
end

function S=smooth2(M,nx,ny)

S=[repmat(M(1,:),nx,1)' M' repmat(M(end,:),nx,1)']';
S=[repmat(S(:,1),1,ny) S repmat(S(:,end),1,ny)];
S=smoothc(S,nx-1,ny-1);

end

function mO = smoothc(mI, Nr, Nc)

% SMOOTHC.M: Smooths matrix data, cosine taper.
% MO=SMOOTHC(MI,Nr,Nc) smooths the data in MI
% using a cosine taper over 2*N+1 successive points, Nr, Nc points on
% each side of the current point.
%
% Inputs: mI - original matrix
% Nr - number of points used to smooth rows
% Nc - number of points to smooth columns
% Outputs:mO - smoothed version of original matrix
%
%
if nargin<2, error('Not enough input arguments!'), end;

% Determine convolution kernel k
Nr=Nr+1;
Nc=Nc+1;
kr=2*Nr+1;
kc=2*Nc+1;
midr=Nr+1;
midc=Nc+1;
maxD=(Nr.^2+Nc.^2).^0.5;
for irow=1:kr;
    for icol=1:kc;
        D=((midr-irow).^2+(midc-icol).^2).^(0.5);
        k(irow,icol)=cos(D*pi/2./maxD);
    end;
end;

k = k./sum(k(:));

% Perform convolution
mO=conv2(mI,k,'valid');
end


function [df]=LinearC(T,c)

precision=0.01;

for i=1:length(c)
    ct=0;%mise à zero du compteur
    w=2*pi/T;
    k=w/c(i);
    g=9.81;%gravite
    do=1000; %valeur quelconque grande
    d=c(i)^2/g;
    while(abs(do-d)>precision)
        ct=ct+1;
        do=d;
        dispe=w^2-g*k*tanh(k*d);
        fdispe=-g*(k^2)./(cosh(k*d)^2);
        d=d-dispe./fdispe;
    end
    df(i)=d;
    
end
end

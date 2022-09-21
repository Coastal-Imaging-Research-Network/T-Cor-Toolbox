function [B]=T_cor_funct3(A,icmin,icmax,dt,resc,methodPreT)
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
A3=A;



if(methodPreT==0)
B=double(A);
elseif(methodPreT==1)
% A3(1:nt,icmin:resc:icmax )=zeroavg(double(A(1:nt,icmin:resc:icmax)));
% per=PeriodWaves(A3,icmax-50,icmax,dt,2);
for ic=icmin:resc:icmax
    
 Tcoupure=1.5; %frecuencia de corte
 fr=1/dt;
 Val=(1/Tcoupure)*2*(1/fr);
 ord=1000; fil=fir1(ord,Val,'low');
 kk1=conv(fil,A3(:,ic));
 S= detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));

 
 Tcoupure=20; %frecuencia de corte
 fr=1/dt;
 Val=(1/Tcoupure)*2*(1/fr);
 ord=1000; fil=fir1(ord,Val,'high');
 kk1=conv(fil,S);
 y= detrend(kk1(  (max(ord)/2)+1: length(kk1)-(max(ord)/2)));



 B(:,ic)=y;%Norma(y,2);
end
    elseif(methodPreT==2)
        
        for ic=icmin:4*resc:icmax   
         [hs,htiers,hrms,trms,Hmax,h]=Wave_Char(smooth(double(A(:,ic,1)),2),dt,0,1);   
            T(ic)=trms;
        end
        To=min(T(find(T>3)))
%         To=min(T(find(T>0))) % Greg a modifié la ligne 51 par ligne 52
        
for ic=icmin:resc:icmax        
%Methode 2 : temporelle
v=double(A(:,ic,1));
% B(:,ic)=Norma(v'-slidefun (@mean,24,v'));
B(:,ic)=smooth (v-smooth(v,round(6*To)),round((2*To)/3));
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
% per=PeriodWaves(A2(1:nt,icmin:icmax),10,20,dt,2);
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

B(:,ic)=smooth(v-smooth(v',7),2);
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

% if  methodPreT~=4
% B(find(isnan(B)==1))=0;
% end

ii=find(abs(mean(B))>0&isnan(abs(mean(B)))==0);
%Interpolation spatiale entre icmin et icmax
for irt=1:nt
    try
     B(irt,icmin:icmax)=interp1(ii,B(irt,ii),icmin:icmax);
    end
end




disp('               Pre-Traitement ok')
%  figure(52);
%  clf
% subplot(4,1,1);pcolor(B(1:40,:));shading flat
% % subplot(4,1,2);pcolor(A3);shading flat
%  subplot(4,1,3);plot(double(A(:,1600)));set(gca,'ydir','reverse')
%  subplot(4,1,4);plot(B(:,1600));set(gca,'ydir','reverse')
end

function [hs,htiers,hrms,trms,Hmax,h,Tp,tindiv]=Wave_Char(d,dt,filt,meth)
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

%Outputs
%h : hauteurs de vague individuelles (ex: pour faire un histogramme)
%tindiv: periodes de vagues individuelles (ex: pour faire un histogramme)
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
 x=[1:nt];
 a=reshape(x,nt,1);
 
 % generation des coefficients du polynome de degré un
  p=polyfit(a,k,1);
 
 %Suppression de la partie trendée de la courbe
 
 k2=k-p(1).*a;
 
 % Suppression de la partie moyenne
 km=k2;
 k2=km-mean(km);
 
 clear k4
if(filt==1)
freq = (1./dt)*(1:nt)/nt;
 %frequences de coupure : %%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %fréquence de coupure vers les basses fréqences
 % fcb= 1/Tlimit
  fcb=0.01;
 % fréquence de coupure vers les hautes fréquences
 % fch=1/Tlimit
  fch=1/0.1;
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
 clear z j
 j=0;
z(1:floor(nt./2))=0;

%Mean zero crossing
% zc=1 mean zero crossing
% zc=2 up zero crossing
% zc=3 down zero crossing
% for zc=1:2

if(zc==1)
 for i=2:nt-1
     if(k3(i)*k3(i+1)<=0 && k3(i-1)*k3(i)>=0 ) 
         j=j+1;
        z(j)=i;
    end
 end
 
%  figure(36);plot(k2,'k');hold on;plot(k3,'r');plot(z,0,'r*')
elseif(zc==2)
    clear z j
     j=0;
 for i=2:nt-1
     if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])<0 ) 
         j=j+1;
        z(j)=i;
    end
 end   
elseif(zc==3)
     clear z j
     j=0;
 for i=2:nt-1
     if(k3(i)*k3(i+1)<=0 && mean([k3(i-1) k3(i)])>0 ) 
         j=j+1;
        z(j)=i;
        mean([k3(i-1) k3(i)])
    end
 end    
end

%calcul de la taille du vecteur dans lequel sont stockés les indices des zeros
k=0;
sz=size(z);
for i=1:sz(2)
    if(z(i)~=0) 
        k=k+1;
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
kh=0;
h(1:k-2)=0;

for i=1:k-2
    
 if(hmin(i)~=0 & hmax(i+1)~=0)
        if(hmin(i+1)==0 & hmax(i)==0)        
        kh=kh+1;
     h(kh)=hmax(i+1)-hmin(i);
        end
 end
end

% Calcul de la hauteur rms moyenne des vagues sur la durée de mesure
%hrms1=0;
clear hrms ord htiers

[sd ord]=sort(h(find(h~=0)),'descend');

%Nombre de vagues pour le calcul
% disp(['Nombre de vagues : ',num2str(length(find(h~=0)))])

% Affichage d'un Histogramme
% figure(85);clf
% bar(h(find(h~=0)));
% set(gcf,'Color','w')
% set(gca,'FontSize',14)
% xlabel('H (m)')
% ylabel ('N°')
% grid on
% [val,nval]=hist(h(find(h~=0)));
% figure(86)
% plot(nval,FiltreMean(val,2),'k')
% set(gcf,'Color','w')
% set(gca,'FontSize',14)
% xlabel('H (m)')
% ylabel ('N°')
% grid on

%Hauteur max
Hmax=max(h);

%Hauteur 1/3
htiers=mean(h(ord(1:round(length(ord)/3))));% on prend le premier tiers pour la hauteur rms
%Hauteur rms
hrms=mean(h(1:kh));

%Hauteur significative
hs=4.*std(k4);
% disp(['Hs=',num2str(hs),'; Hrms=',num2str(hrms),'; H1/3=',num2str(htiers),'; Hmax=',num2str(Hmax)])
%Calcul de la periode de la houle pour chaque vague

clear t
t(1:k-1)=0;
for i=1:k-1
    if(z(i)~=0 | z(i+1)~=0 ) 
        t(i)=2*(z(i+1)-z(i))*dt;
    end
end

tindiv=t;
%Calcul de la periode rms de la houle et affichage

if(zc==1)
trms=mean(t);
elseif(zc==2)
trms=mean(t)/2;
elseif(zc==3)
trms=mean(t)/2;
end    


%%%%%Periode Peak%%%%%%%%%%%
y=detrend(d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Nombre de points pour la fft
nx=length(y);
%Pas de temps
% x=(1:nx);
%Resolution temporelle (periode d'échantillonnage) (en s)
% dt=1;

%FFT en précisant la taille du vecteur
z=fft(y);
[n1 n2]=size(z);
if n1==1
    z=z';
end

t=1./((1:nx/2)./(dt*nx));
in=find(t<6|t>14);
x=1:nx/2;
k=z(x).*conj(z(x)/nx);
k(in)=0;

[sa ind]=lmax(FiltreMean(FiltreMean(k,3),5));


%  disp(ind./(dt*nx))

% io=find(sa>max(sa)/2);
Tp=1./(ind./(dt*nx));
end

function [lmval,indd]=lmax(xx,filt)
%LMAX 	[lmval, indd]=lmax(xx,filt). Find local maxima in vector XX,where
%	LMVAL is the output vector with maxima values, INDD  is the 
%	corresponding indexes, FILT is the number of passes of the small
%	running average filter in order to get rid of small peaks.  Default
%	value FILT =0 (no filtering). FILT in the range from 1 to 3 is 
%	usially sufficient to remove most of a small peaks
%	For example:
%	xx=0:0.01:35; y=sin(xx) + cos(xx ./3); 
%	plot(xx,y); grid; hold on;
%	[b,a]=lmax(y,2)
%	 plot(xx(a),y(a),'r+')
%	see also LMIN, MAX, MIN
	
%**************************************************|
% 	Serge Koptenko, Guigne International Ltd., |
%	phone (709)895-3819, fax (709)895-3822     |
%--------------06/03/97----------------------------|

x=xx;
len_x = length(x);
	fltr=[1 1 1]/3;
  if nargin <2, filt=0; 
	else
x1=x(1); x2=x(len_x); 
	for jj=1:filt,
	c=conv(fltr,x);
	x=c(2:len_x+1);
	x(1)=x1;  
        x(len_x)=x2; 
	end
  end
lmval=[]; indd=[];
i=2;		% start at second data point in time series
    while i < len_x-1,
	if x(i) > x(i-1)
	   if x(i) > x(i+1)	% definite max
lmval =[lmval x(i)];
indd = [ indd i];
	   elseif x(i)==x(i+1)&x(i)==x(i+2)	% 'long' flat spot
%lmval =[lmval x(i)];  	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite max included
i = i + 2;  		% skip 2 points
	   elseif x(i)==x(i+1)	% 'short' flat spot
%lmval =[lmval x(i)];	%1   comment these two lines for strict case
%indd = [ indd i];	%2 when only  definite max included
i = i + 1;		% skip one point
	   end
	end
	i = i + 1;
    end
if filt>0 & ~isempty(indd),
	if (indd(1)<= 3)|(indd(length(indd))+2>length(xx)), 
	   rng=1;	%check if index too close to the edge
	else rng=2;
	end
	  for ii=1:length(indd), 	% Find the real maximum value
	    [val(ii) iind(ii)] = max(xx(indd(ii) -rng:indd(ii) +rng));
	    iind(ii)=indd(ii) + iind(ii)  -rng-1;
	  end
  indd=iind; lmval=val;
else
end
end

function [N,Ma]=zeroavg(M,c)

% function [N,Ma]=zeroavg(M,c)
%
% This function computes the average of each column of a matrix M
% and subtracts it % from every entry in that column. 
% If the column contains a time % serie, the time serie will become 
% zero averaged. 
%
% If an extra argument "c"
% is passed to the function 
% (independent of its value) the averageing will
% not be performed over its columns bit over its rows.
% The zero averaged matrix is returns in matrix N.

[p,q]=size(M);

if nargin == 1
   Ma=sum(M)/p;
   N = M - repmat(Ma,p,1);
end
if nargin>1
   Ma=(sum(M')/q)';
   N=M-repmat(Ma,1,q);
end
end

% Normalise un vecteur presentant 
%une variation d'amplitude par une moyenne glissante
%USAGE [f2]=Norm(f,n)
function Normf=Norma(f,n)

F=f;
version=3;
%n Valeur pour la moyenne glissante
%n=3;

%valeur absolue de f

Normf=F.*0;

if(version==1)
%
clear fct2

for it=1:length(f)
        for i=1:size(F,2)
        f=F(:,i);

    if(it<=n)
  Normft(it)=f(it)./nanmean(abs(f(it:it+n)));
elseif(it>=length(f)-n)
Normft(it)=f(it)./nanmean(abs(f(it-n:it)));
else
  Normft(it)=f(it)./nanmean(abs(f(it-n:it+n)));  
    end
Normf(1:length(Normft),i)=Normft;
        end
end
elseif(version==2)

    for i=1:size(F,2)
        f=F(:,i);
    % Version matricielle
f2=[f(1) ;f(1); f; f(length(f)); f(length(f))];   
absf=[abs(f(1)) ;abs(f(1)); abs(f); abs(f(length(f))); abs(f(length(f)))];
absp2=[abs(f(1)); abs(f(1)); abs(f(1)) ;abs(f(1)); abs(f)];
absp1=[abs(f(1)) ;abs(f(1)); abs(f(1)); abs(f); abs(f(length(f)))];
absm1=[abs(f(1)) ;abs(f) ;abs(f(length(f))) ;abs(f(length(f))); abs(f(length(f)))];
absm2=[abs(f) ;abs(f(length(f))); abs(f(length(f))); abs(f(length(f))); abs(f(length(f)))];

f2=5*f2./((absf+absp2+absp1+absm1+absm2));
Normft=f2(3:(length(f2)-2));
Normf(:,i)=Normft;
    end
    
    
elseif(version==3)
    
for ir=1:size(F,2)
t22=hilbert(F(:,ir));
Normf(:,ir)=F(:,ir)./smooth(sqrt(t22.*conj(t22)),18*2);
end
    
end %if version

% >> figure;plot(f2)
% >> hold on;plot(f,'r')

end
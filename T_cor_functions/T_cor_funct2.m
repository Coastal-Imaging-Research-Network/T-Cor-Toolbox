function [Ii,res]=T_cor_funct2(Img,dxx,res,show)

% Inputs:
%       Img : Image (in uint8)
%       dxx : Pixel footprint vector
%       show = 1 display figure

% outputs:
%       Ii : Image rectifiée
%       res : résolution en mètre de l'image rectifiée

    if nargin<4
        show=0;
        if nargin<3
            res=1;
        end
    end

    % freqt=input('     Entrer la fréquence d''acquisition en Hz (exple GPP: 2 Hz): ');
    freqt=2; % GPP: 2 Hz

    % Compute timestack image in real coordinates
    b=cumsum(dxx(end:-1:1)); % Offshore to the beach
    a=1:size(Img,1); % Time (2 Hz)
    
%     res=input('     Résolution spatiale voulue en mètre (exple 1 m): ');
%     res=1; % 1 mètre de résolution
    
    xi=ceil(min(b)):res:fix(max(b)); % Cross-shore (m)
    yi=a; % Time (s)
    Ii=NaN(length(yi),length(xi),3); % Timestack in real coordinates

    for i=1:3 % For RGB
        for j=1:length(yi)
            clear v
            [bb,ibb]=unique(b);
            v=double(Img(j,ibb,i));
            %         Ii(j,:,i)=interp1(bb,v,xi,'linear');
            Ii(j,:,i)=interp1(bb,v,xi,'nearest');
        end
    end
    
    % Contrôle des limites de bandes
    Ii(Ii>255)=255;
    Ii(Ii<0)=0;
    
    % Passage en uint8
    Ii=uint8(Ii);
    
    if show==1
        figure(99)
        subplot(1,2,1);
        image(1:size(Img,2),1:size(Img,1)/freqt,Img)
        xlabel('Pixels','fontsize',16)
        ylabel('Time (s)','fontsize',16)
        title('Original','FontSize',14)
        subplot(1,2,2);
        image((1:size(Ii,2)).*res,1:size(Ii,1)/freqt,Ii)
        xlabel('X (m)','fontsize',16)
        ylabel('Time (s)','fontsize',16)
        title('Rectified','FontSize',14)
    end
end
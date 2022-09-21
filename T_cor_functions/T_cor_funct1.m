function [dx]=T_cor_funct1(sk2,sk1,path_img,path_rect,p)

% Cette fonction permet de calculer les distances en mètres entre pixels
% consécutifs
% Greg, 04 Mars 2017

lmax=max(abs(sk2(2)-sk1(2)),abs(sk2(1)-sk1(1)));
pas1=(sk2(1)-sk1(1))/lmax;
pas2=(sk2(2)-sk1(2))/lmax;
s1=round(sk1(1):pas1:sk2(1));
s2=round(sk1(2):pas2:sk2(2));
[Vcoord1,Vcoord2]=PixtoCoord2Taller(s2,s1,zeros(1,length(s2)),path_rect);

if nargin==4
    p=1;
end

if p==1
    I=imread(path_img);
    figure;
    image(I)
    hold on
    plot(s2,s1,'.','Linewidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',0.5)
    plot(s2(1),s1(1),'ok','MarkerSize',5)
    text(s2(1)-10,s1(1)-10,'Start')
    plot(s2(end),s1(end),'ok','MarkerSize',5)
    text(s2(end)-10,s1(end)-10,'End')
end

Co1=Vcoord2;Co2=Vcoord1;

dx=abs(sqrt(diff(Co1).^2+diff(Co2).^2));
end

function [Vcoord1,Vcoord2]=PixtoCoord2Taller(Vpix1,Vpix2,z,path_rect)
%Modif le 20 Oct ajout de z en input (coordonnée verticale du plan d'eau)
% Cette fonction permet de calculer les coordonnées GPS correspondantes aux
% coordonnées pixels avec élévation par rapport au plan d'eau.
% Reçue de Rafael Almar

clear Vcoord1 Vcoord2
load(path_rect)

%Version avec altitude
[X_I] = rectifyTaller(double([Vpix1; Vpix2]),Rckk,Tckk,fc,cc,kc,alpha_c,z);

Vcoord1=X_I(1,:);
Vcoord2=X_I(2,:);

end

function [X_r] = rectifyTaller(aaa,Rckk,Tckk,fc,cc,kc,alpha_c,Z)

% aaa contains the pixel coordinates (origin top left, horizontal in the
% first row vertical in the second row).
%
% Z is the level you want to recitfy to (can be
% different for each point, in which case length(Z) must be the same as X_kk))

[n,m]=size(aaa);
if (m==2 && n~=2)
    x_kk=aaa';
end
[n,m]=size(aaa);
if length(Z)==1
    Z=Z.*ones(1,m);
end

%undistort the image coordinates and scale
xn = normalizeTaller(aaa,fc,cc,kc,alpha_c);


xn=[xn;ones(1,length(xn))];

%rotate the normalised coordinates
R=Rckk'*xn;

%rotate the camera location
T=Rckk'*Tckk*ones(1,length(xn));

%figure out what the scale factor is
z=(Z+T(3,:))./R(3,:);

%apply the scale factor to rotated coordinates, and correct for camera
%position.
X_r=[z.*R(1,:)-T(1,:);z.*R(2,:)-T(2,:);Z];
end

function [xn] = normalizeTaller(x_kk,fc,cc,kc,alpha_c)

%normalize
%
%[xn] = normalize(x_kk,fc,cc,kc,alpha_c)
%
%Computes the normalized coordinates xn given the pixel coordinates x_kk
%and the intrinsic camera parameters fc, cc and kc.
%
%INPUT: x_kk: Feature locations on the images
%       fc: Camera focal length
%       cc: Principal point coordinates
%       kc: Distortion coefficients
%       alpha_c: Skew coefficient
%
%OUTPUT: xn: Normalized feature locations on the image plane (a 2XN matrix)
%
%Important functions called within that program:
%
%comp_distortion_oulu: undistort pixel coordinates.

if nargin < 5,
   alpha_c = 0;
   if nargin < 4;
      kc = [0;0;0;0;0];
      if nargin < 3;
         cc = [0;0];
         if nargin < 2,
            fc = [1;1];
         end;
      end;
   end;
end;


% First: Subtract principal point, and divide by the focal length:
x_distort = [(x_kk(1,:) - cc(1))/fc(1);(x_kk(2,:) - cc(2))/fc(2)];

% Second: undo skew
x_distort(1,:) = x_distort(1,:) - alpha_c * x_distort(2,:);

if norm(kc) ~= 0,
	% Third: Compensate for lens distortion:
	xn = comp_distortion_ouluTaller(x_distort,kc);
else
   xn = x_distort;
end
end

function [x] = comp_distortion_ouluTaller(xd,k)

%comp_distortion_oulu.m
%
%[x] = comp_distortion_oulu(xd,k)
%
%Compensates for radial and tangential distortion. Model From Oulu university.
%For more informatino about the distortion model, check the forward projection mapping function:
%project_points.m
%
%INPUT: xd: distorted (normalized) point coordinates in the image plane (2xN matrix)
%       k: Distortion coefficients (radial and tangential) (4x1 vector)
%
%OUTPUT: x: undistorted (normalized) point coordinates in the image plane (2xN matrix)
%
%Method: Iterative method for compensation.
%
%NOTE: This compensation has to be done after the subtraction
%      of the principal point, and division by the focal length.


if length(k) == 1,
    
    [x] = comp_distortionTaller(xd,k);
    
else
    
    k1 = k(1);
    k2 = k(2);
    k3 = k(5);
    p1 = k(3);
    p2 = k(4);
    
    x = xd; 				% initial guess
    
    for kk=1:20,
        
        r_2 = sum(x.^2);
        k_radial =  1 + k1 * r_2 + k2 * r_2.^2 + k3 * r_2.^3;
        delta_x = [2*p1*x(1,:).*x(2,:) + p2*(r_2 + 2*x(1,:).^2);
        p1 * (r_2 + 2*x(2,:).^2)+2*p2*x(1,:).*x(2,:)];
        x = (xd - delta_x)./(ones(2,1)*k_radial);
            
    end
    
end

end

function [x_comp]  = comp_distortionTaller(x_dist,k2)

%       [x_comp] = comp_distortion(x_dist,k2);
%       
%       compensates the radial distortion of the camera
%       on the image plane.
%       
%       x_dist : the image points got without considering the
%                radial distortion.
%       x : The image plane points after correction for the distortion
%       
%       x and x_dist are 2xN arrays
%
%       NOTE : This compensation has to be done after the substraction
%              of the center of projection, and division by the focal
%              length.
%       
%       (do it up to a second order approximation)


[two,N] = size(x_dist);


if (two ~= 2 ), 
    error('ERROR : The dimension of the points should be 2xN');
end;


if length(k2) > 1,
    
    [x_comp]  = comp_distortion_ouluTaller(x_dist,k2);
    
else
    
    radius_2= x_dist(1,:).^2 + x_dist(2,:).^2;
    radial_distortion = 1 + ones(2,1)*(k2 * radius_2);
    radius_2_comp = (x_dist(1,:).^2 + x_dist(2,:).^2) ./ radial_distortion(1,:);
    radial_distortion = 1 + ones(2,1)*(k2 * radius_2_comp);
    x_comp = x_dist ./ radial_distortion;
    
end
end
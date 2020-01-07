clear;
radius=32*10^-3;
focus_distance=63.2*10^-3;
Transducer = get_spherical_shell(radius, focus_distance);
%Transducer=create_spherical_shell_planar_array(1, 1, radius, focus_distance ,.01,.01);
%draw_array(Transducer);
Transducer.amplitude = 0.2;

f0 = 1.16e6;
xmin = 0;
xmax = 32*10^-3;
ymin = 0;
ymax = 32*10^-3;
zmin = 0.4*63.2*10^-3; % Don't capture the pressure field inside the transducer array
zmax = 1.6*63.2*10^-3;

dx = 0.2e-3;
dy = dx;
dz = 0.2e-3;

x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;
delta = [dx dy dz];
coord_grid = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);

medium = set_medium('nonlinearlossy');
medium.attenuationdBcmMHz = 0.022;%.1;%7e-4*689.0284;%1;
medium.nonlinearityparameter = 5.2;
medium.powerlawexponent=2; %http://iopscience.iop.org/article/10.1088/0031-9155/58/17/5833

pistonPts = 200;
RhoMax = 1;
dSigFD = 0.0025;
NHARM = 10; % Number of Harmonics
Exflag = 0; % no effect

% Perform the calculation
tic();
pkzk_cw=kzk_cw(Transducer, coord_grid, medium, f0 , NHARM, dSigFD, pistonPts, RhoMax, Exflag);
t_kzk2 = toc();
fprintf('KZK finished in %f seconds\n', t_kzk2);

%% 
% P=zeros(sz(1:3));
% for ii=1:NHARM
%     P=P+pkzk_cw(:,1,:,ii);
% end
figure;
title('First Harmonic');
%mesh(z,x,squeeze(abs(P(:,1,:))));
mesh(z,x,squeeze(abs(pkzk_cw(:,1,:,8))));
xlabel('axial distance (m)')
ylabel('radial distance (m)')

%% 
nHarm=NHARM;
P=pkzk_cw;
l=2000;
clear A
clear t
clear Ptemp
t=(1/f0:-1/f0/(l-1):0)*10;
f=([1:NHARM]*f0)';
Preal=squeeze(real(P(:,:,:,1)));
[I,J]=find(Preal==max(Preal(:)));

for ii=1:length(f)
    for jj=1:length(t)
        A(ii,jj)=exp(1i*pi*2*f(ii)*t(jj));
        Ptemp(ii,jj)=A(ii,jj)*P(I,:,J,ii);
        
    end
end
figure;
plot(real(sum(Ptemp)));

%% 
clear A
clear t
clear Ptemp
l=100;
t=(1/f0:-1/f0/(l-1):0);
P=squeeze(P);
sz=size(P);
Pall=zeros(sz(1),sz(2),1);
for aa=1:sz(1)
    for bb=1:sz(2)
        for ii=1:length(f)
            for jj=1:length(t)
                A(ii,jj)=exp(1i*pi*2*f(ii)*t(jj));
                Ptemp(aa,bb,ii,jj)=A(ii,jj)*P(aa,bb,ii);
            end
        end        
    end
    disp(aa);
end
Psum = squeeze(sum(Ptemp,3));

%mesh(z,x,abs(Pall));
    
%% rotate the data
Psum = permute(Psum, [2,1,3]);
for jj = 1:length(t)
    P3d(:,:,:,jj) = revolve2D(Psum(:,:,jj));
end
clc;
clear;

%% Size of 3D model.
nx=512;
ny=512;
nz=512;

%% Creat vel & den 3D model.
vp=3500*ones(nx,ny,nz);
 
dens=1200*ones(nx,ny,nz);
 
fid1=fopen('vp512.bin','wb');
fwrite(fid1,vp,'float');

fid2=fopen('rho512.bin','wb');
fwrite(fid2,dens,'float');
fclose(fid1);
fclose(fid2);

%% Display the Marmousi models.
for z=1:nz
    vp_slice(:,:,z)=vp(:,20,z);
    den_slice(:,:,z)=dens(:,20,z);
end

figure(1);
imagesc(vp_slice(:,:)');
title('A slice of velocity model');
figure(2);
imagesc(den_slice(:,:)');
title('A slice of density model');

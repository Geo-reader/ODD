clc;
clear;

%% Size of 3D model.
nx=101;
ny=101;
nz=101;

%% Creat vel & den 3D model.
vp=4000*ones(nx,ny,nz);

vs=2300*ones(nx,ny,nz);
 
dens=1200*ones(nx,ny,nz);
 
fid1=fopen('vp.bin','wb');
fwrite(fid1,vp,'float');

fid2=fopen('vs.bin','wb');
fwrite(fid2,vs,'float');

fid3=fopen('rho.bin','wb');
fwrite(fid3,dens,'float');
fclose(fid1);
fclose(fid2);
fclose(fid3);

%% Display the Marmousi models.
for z=1:nz
    vp_slice(:,:,z)=vp(:,20,z);
    vs_slice(:,:,z)=vs(:,20,z);
    den_slice(:,:,z)=dens(:,20,z);
end

figure(1);
imagesc(vp_slice(:,:)');
title('A slice of p-velocity model');
figure(2);
imagesc(vs_slice(:,:)');
title('A slice of s-velocity model');
figure(3);
imagesc(den_slice(:,:)');
title('A slice of density model');

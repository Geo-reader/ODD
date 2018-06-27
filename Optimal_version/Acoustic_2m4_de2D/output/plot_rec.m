clc;
clear;

%% Size of 3D model.
nx=101;
ny=101;
nz=101;
tmax=1001;

%% Some parameters.
dt=1e-3;
dx=10;
dy=10;
tt=0:dt:(tmax-1)*dt;
xx=0:dx:(nx-1)*dx;
yy=0:dy:(ny-1)*dy;

%% Open & read seismic record data.
fid=fopen('rec_whole.bin','r');

%% Display the seismic record.
for t=1:tmax
     rec=fread(fid,[nx,ny],'float');
     rec_slice(t,:)=rec(:,fix(ny/2));
  %   figure(8);imagesc(xx,yy,rec);colormap gray;title('wave field record');colorbar;
     t
end
figure(12)
imagesc(yy,tt,rec_slice)
title('A common gather record');
xlabel('distance (m)');
ylabel('Time (s)'); 

ma=max(max(rec_slice));
mi=min(min(rec_slice));
caxis([ mi/100, ma/100]); 
fclose(fid);



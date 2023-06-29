clc
close all
[filename,pathname]=uigetfile('*.tif','Select data file');
cd(pathname);
[filepath,name,ext] = fileparts(filename);
im=imread(filename);
im2=im; %e16 midgut
% se=strel("disk",20); % e14 hindgut
% se=strel("disk",8);% e14 midgut

% im2=imtophat(im,se); % e14 midgut, hindgut
% figure, imshow(im2);

% im3=imgaussfilt(im2,3); % e14 hindgut
% im3=imgaussfilt(im2,5); %e14 midgut
im3=imgaussfilt(im2,1); %e16 midgut

% T=adaptthresh(im3,0.5,'ForegroundPolarity','dark'); % e14 hindgut
% T=adaptthresh(im3,0.5,'ForegroundPolarity','dark'); %e14 midgut
T=adaptthresh(im3,0.45,'ForegroundPolarity','dark'); %e16 midgut

imbw=imbinarize(im3,T);
% imbw2=bwmorph(imbw,'thin',10); % e14 midgut
% se2=strel('disk',3); % e14 midgut
% imbw2=imclose(imbw,se2); % e14 midgut
% figure, imshow(imbw2); 

imbw2=imcomplement(bwareaopen(imcomplement(imbw),500)); %e16 midgut
% imbw2=imcomplement(bwareaopen(imcomplement(imbw),200)); % e14 hindgut
% imbw3=imcomplement(bwareaopen(imcomplement(imbw2),500)); % e14 midgut
imbw2=bwareaopen(imbw2,100); % e14 midgut
figure, imshow(im)
figure, imshow(imbw2)
% saveas(gcf,'e16_midgut_binary.tif')

%%

% b=imresize(imbw2,0.5);
% % b=imcomplement(imbw2);
% % b=imbw4;
% % b=imbinarize(imcomplement(i)); 
% % b=imbinarize(i);
% % b2=bwskel(b); % e14 hindgut
% b2=bwskel(b,'MinBranchLength',20); % e14 midgut
% b3=bwmorph(b2,'branchpoints'); 
% figure, imshow(b2);
% % saveas(gcf,'e14_hindgut_skel.tif')
% % figure, imshow(b3);
% branches=sum(b3(:)==1);

 %%
 %%
% b=imbw2;
b=imresize(imbw2,1.66);
f=im2double(b);
f=imgaussfilt(f,1);
[x,y]=meshgrid(-(length(f)/2):((length(f)/2)-1),-(length(f)/2):((length(f)/2)-1));
z=sqrt(x.^2+y.^2);
c=z<(length(f)/2);
f=f.*c;
% figure, imshow(f);


df=1/length(f); %frequency resolution
sampleIndex = (-length(f)/2):(length(f)/2)-1; %raw index for FFT plot
x=sampleIndex*df; %x-axis index converted to frequencies
y=sampleIndex*df;


F=fftshift(fft2(f)/length(f));
I=mat2gray(log(abs(F)+5));

figure, imagesc(x,y,I);

axis equal; axis tight; colormap(jet); 
c=colorbar;
ylabel(c,'Normalized Amplitude','Rotation',360);

set(gca,'TickDir','out')
set(gca,'TickLength',[0.02,0.02])
set(gca,'LineWidth',1.5)
set(gca,'FontSize',20)
xlim([-0.05 0.05]);
xticks([-0.025 0 0.025]);
ylim([-0.05 0.05]);
yticks([-0.025 0 0.025]);
xlabel('Frequency (1/um)'); ylabel('Frequency (1/um)');
% saveas(gcf,'mgctl_2_FFT.pdf')




%%
% I=im2uint8(I);
% figure, imshow(I);
% colorbar
% imshow(f)
% F = fft(f);
% F2 = log(abs(F));
% figure, imagesc(abs(fftshift(F)))
% imshow(F);
% colormap(jet); colorbar
% F = fft2(f,256,256);
% imshow(log(abs(F)),[-1 5]); colormap(jet); colorbar
% F = fft2(f,256,256);F2 = fftshift(F);
% imshow(log(abs(F2)),[-1 5]); colormap(jet); colorbar
% 
% P = peaks(20);
% X = repmat(P,[5 10]);
% imagesc(X)
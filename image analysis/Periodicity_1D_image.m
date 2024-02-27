clear -pathname;
format long;
[filename,pathname]=uigetfile('*.png','Select data file');
cd(pathname);
[filepath,name,ext] = fileparts(filename);
calib = 1.66;     % in um
RGB = imread(filename);
Full_list = rgb2gray(RGB);
figure;
imshow(Full_list);
[d1,d2] = size(Full_list);
Length = d2*calib;
X_right = d2*calib;
X_left = 1*calib;
Delta_x = 2;                         % Bin size in um
Bins = floor(Length./Delta_x);
dx = floor(Delta_x./calib);
Periodicity = zeros(Bins,1);
Interval = zeros(Bins, 1);

for i = 1:d2
    Sum_y(i,1) = sum(Full_list(:,i));
end;

for i = 1:Bins
    Periodicity(i,1) = Periodicity(i,1) + sum(Sum_y((i-1)*dx+1:i*dx));
    Interval(i,1) = i*Delta_x;
end;

%Get intensity-position plot
Norm_factor = mean(Periodicity);
norm_periodicity = Periodicity(:,1)./Norm_factor;
% Periodicity_2 = smooth(Periodicity);
figure; 
plot(Interval(:,1),norm_periodicity,'k','LineWidth',3)
%axis([0,Length,0,1]);
xlabel('Position (um)','FontSize',14,'FontWeight','bold','Color','black')
ylabel('Normalized Intensity (a.u.)','FontSize',14,'FontWeight','bold','Color','black')
set(gca,'XTick',(0:10*Delta_x:Length))
%set(gca,'YTick',(0:0.1:1));
% saveas(gcf,'hind_intensity_long4.png') %********************************************************

xlswrite([name,'.xls'],Interval(:,1),1,'A1')
xlswrite([name,'.xls'],norm_periodicity,1,'B1')

%Get Fourier transform
%n = 2^nextpow2(Bins)
Y = fft(Periodicity);
P2 = abs(Y/Bins);
P1 = P2(1:Bins/2+1);
P1(2:end-1) = 2*P1(2:end-1);
%P1 = P1./max(P1);
f = 1/Delta_x*(0:Bins/2)/Bins;
figure;
plot(f,P1,'k','LineWidth',1.5)
%title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('Frequency (1/um)','FontSize',14,'FontWeight','bold','Color','black')
ylabel('Norm. Amplitude','FontSize',14,'FontWeight','bold','Color','black')

% saveas(gcf,'hind_fourier_long4.png') %*************************************************

xlswrite([name,'.xls'],f.',2,'A1')
xlswrite([name,'.xls'],P1,2,'B1')

%Get most probable frequency
b = max(P1(2:end));
c = max(max(P1(2:find(P1 == b)-1)),max(P1((find(P1 == b)+1):end)));
F_max = (find(P1 == b)-1)/Delta_x/Bins;
T = 1./F_max
F_goodness = b/c;
xlswrite([name,'.xls'],T,4,'A1')

%Get autocorrelation plot and amplitude
acf = autocorr(Periodicity,Bins-1);
int=Interval(:,1)-Delta_x;
figure;
plot(int,acf,'k','LineWidth', 3);
xlabel('Lag (um)','FontSize',14,'FontWeight','bold','Color','black')
ylabel('Normalized Autocorrelation Amplitude','FontSize',14,'FontWeight','bold','Color','black')
set(gca,'XTick',(0:10*Delta_x:Length))
set(gca,'YTick',(-1:0.1:1));

% saveas(gcf,'hind_autocorr_long4.png') %********************************************************************

xlswrite([name,'.xls'],Interval(:,1)-Delta_x,3,'A1')
xlswrite([name,'.xls'],acf,3,'B1')

maxpks=islocalmax(acf);
minpks=islocalmin(acf);
maxamps=acf(maxpks);
minamps=acf(minpks);
figure,plot(int,acf,int(minpks),acf(minpks),'r*',int(maxpks),acf(maxpks),'b*')
autoamp=(maxamps(1)-mean([minamps(1);minamps(2)],1))/2
xlswrite([name,'.xls'],autoamp,5,'A1')


% i = 1;
% while(acf(i)>acf(i+1))
%     i = i+1;
% end;
% Auto_amp = max(acf(i:end));  %need to verify according to the plot!!

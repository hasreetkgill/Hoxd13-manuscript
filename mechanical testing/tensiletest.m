function [name, datastruct] = tensiletest(scale, stepsize, l, kb, w)
% for j = 1:nvids % multiple videos
    clearvars -except scale stepsize l kb ri ro w name datastruct j %nvids 

    [filename, pathname] = uigetfile({'*.avi;','All Image Files';...
              '*.*','All Files' },'mytitle');
%     name = strsplit(pathname,filesep);
%     name = name{end-1};
    name = filename(1:(end-4));
    vidObj = VideoReader(filename); 
    vidHeight = vidObj.Height;
    vidWidth = vidObj.Width;
  
    mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),... 
        'colormap',[]);
    k = 1;
    while hasFrame(vidObj);
        mov(k).cdata = readFrame(vidObj);
        k = k+1;
    end

    [temp, framenum]=size(mov);

    % calculate stress and strain
    videoFileReader = vision.VideoFileReader(filename); 
    videoPlayer = vision.VideoPlayer();
    objectFrame = videoFileReader();

    figure; imshow(objectFrame);
%     if j==1
%     msgbox('click on: A. 3 pairs of points on tissue and B. 1 point to track on lever');
%     else
%     end
    points = ginput(7); 
    pointImage = insertMarker(objectFrame,points,'+','Color','white');
    figure; imshow(pointImage);

    % Track points
    tracker = vision.PointTracker('MaxBidirectionalError',1);
    initialize(tracker,points,objectFrame);

    % for regular tests
%   initialpos = points(7,2);
%     dist(1,1) = sqrt((points(2,1)-points(1,1))^2+(points(2,2)-points(1,2))^2);
%     dist(2,1) = sqrt((points(4,1)-points(3,1))^2+(points(4,2)-points(3,2))^2);
%     dist(3,1) = sqrt((points(6,1)-points(5,1))^2+(points(6,2)-points(5,2))^2);

    i = 1;
    while ~isDone(videoFileReader);
      frame = videoFileReader();
      [points,validity] = tracker(frame); % Tracker determines new location of points per frame
      disp(i+1) = points(7,2) - initialpos;
      dist(1,i+1) = sqrt((points(2,1)-points(1,1))^2+(points(2,2)-points(1,2))^2); % Store distance between points per frame
      dist(2,i+1) = sqrt((points(4,1)-points(3,1))^2+(points(4,2)-points(3,2))^2);
      dist(3,i+1) = sqrt((points(6,1)-points(5,1))^2+(points(6,2)-points(5,2))^2);
      out = insertMarker(frame,points(validity, :),'+');
      videoPlayer(out); % Displays video with positions of points marked in green
      i = i+1;
    end

%     applied(1) = 0;
%     for i =2:framenum;
%         applied(i)=applied(i-1)+stepsize;
%     end

    disp = disp*scale;
    deflection = abs(disp-applied'); 
    force = kb*deflection; % (Pa*mm^2)
%     datastruct.stress = force/(2*w*l); % for ring
%     datastruct.stress = force/(pi*(w*l)^2); % for cylinder eg. foregut
    datastruct.stress = force/((pi*(ro)^2)-(pi*(ri)^2)); % for whole tube
    datastruct.force = force;
%     datastruct.stress{j} = force/(2*w*l); 
%     datastruct.force{j} = force;

    dist = dist*scale;
    strain(1,:) = (dist(1,:)-dist(1,1))./dist(1,1);
    strain(2,:) = (dist(2,:)-dist(2,1))./dist(2,1);
    strain(3,:) = (dist(3,:)-dist(3,1))./dist(3,1);
    avgdist = mean(dist,1);
    datastruct.avgnormdist = avgdist - (avgdist(1));
    datastruct.avgstrain = mean(strain,1);
    datastruct.time = (0.5:0.5:(0.5*framenum));
%     datastruct.avgnormdist{j} = avgdist - (avgdist(1));
%     datastruct.avgstrain{j} = mean(strain,1);
    stdvstrain = std(strain,1);
    sestrain = stdvstrain./sqrt(length(datastruct.avgstrain));
    extracteddata = [strain;datastruct.avgstrain;sestrain;datastruct.stress];
    % extracteddata = [strain;datastruct.avgstrain{j};sestrain;datastruct.stress{j}];

%     titles = {'strain 1';'strain 2';'strain 3';'average strain';'standard error strain';'stress'};
%     xlswrite([name,'.xls'],titles,1,'A1')
%     xlswrite([name,'.xls'],extracteddata,1,'B1')
%     xlswrite([name,'.xls'],titles,j,'A1')
%     xlswrite([name,'.xls'],extracteddata,j,'B1')
% end

end


% calculate lever velocity from a video of lever pulling without sample
function disp = levercal(scale)

[filename, pathname] = uigetfile({'*.avi;','All Image Files';... % choose lever calibration video
          '*.*','All Files' },'mytitle');
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

videoFileReader = vision.VideoFileReader(filename); 
videoPlayer = vision.VideoPlayer();
objectFrame = videoFileReader();

figure; imshow(objectFrame);
msgbox('click on point to track on lever');
points = ginput(1); 
pointImage = insertMarker(objectFrame,points,'+','Color','white');
figure; imshow(pointImage);

tracker = vision.PointTracker('MaxBidirectionalError',1); % click on tip of lever
initialize(tracker,points,objectFrame);
initialpos = points(1,2); 

disp(1) = points(1,2) - initialpos;

i = 1;
while ~isDone(videoFileReader);
  frame = videoFileReader();
  [points,validity] = tracker(frame); % Tracker determines new location of points per frame
  disp(i+1) = points(1,2) - initialpos;
  out = insertMarker(frame,points(validity, :),'+');
  videoPlayer(out); % Displays video with positions of points marked in green
  i = i+1;
end

disp = disp*scale;

end


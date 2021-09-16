%This is the dynamic object demonstration
%09152021
% Clear the window and the buffle:
sca;
close all;
clearvars;
commandwindow
%% ----------------Set up screen parameters--------------------%
ScreenWidthInCM = 29;% in cm; 
ScreenHeightInCM = 18;% in cm; 
ViewingDistance = 30;% in cm

%% ----------------Some user control for debug mode
%displaywindow = [0 0 720 450]; %Smaller window for debugging
displaywindow = []; %uncomment it if want to use the full screen

%% ----------------Set up Fixation Point parameters--------------------%
FixationPointColor = [255 0 0];%Fixation point depicted in red;
FixationPointSizePix = 20;%Fixation point diameter,in pixel;

%% ----------------Set up Object parameters--------------------%

%Object locations (in degree),represented in polar coordinate;
ObjectViewAngle = 180;%degree %Define the righward direction zero;
ObjectEccentricityDeg = 10;%degree  

ObjectSizeDeg = [7, 5];%in degree; This is the initial width and height;
%ObjectSizeCM = [3.5, 2.5];%in cm; This is the initial width and height;

%Object moving parameters
ObjectMotionVelocty = 0.04; %cm/s; positive: zooming; negative: shrinking
TotalMotionTime = 1000; %in ms
PureZooming = 1; %Flag to control motion perspective: 1: zooming on the original corr, 0: towards subject

%Load Object images
ObjectPath = 'Z:\ProgramsInOHLab\DynamicObject\img\';%Path in windows; Surface pro4;
%ObjectImage = dir([ObjectPath '*.png']);
ObjectImage = dir([ObjectPath '*.jpg']);


NumberOfObjects = length(ObjectImage); 
ObjectName = {ObjectImage(:).name};
RandomSelection = randperm(NumberOfObjects);  

NumberOfObjectsInOneTrial = 1;


try
%% -----------------Start setting for the display-------------------%
% Default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
%Skip the synchronize checking temporally for illustration; just temporally....
Screen('Preference', 'SkipSyncTests', 1)
% Seed the random number generator. 
rng('shuffle')
% Get the screen numbers. 
screens = Screen('Screens');
% Always draw on the remote screen;
screenNumber = max(screens);

% Define black and white according to the display of the screen(white will be 1 and black 0).
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
windowRect = Screen('Rect', screenNumber);
screenXpixels = windowRect(3);% Get the length of the on screen window in pixels.
screenYpixels = windowRect(4);% Get the width of the on screen window in pixels.
[xCenterpixels, yCenterpixels] = RectCenter(windowRect);% Get the centre coordinate of the window in pixels
%Set up fixation point location
%if isempty(displaywindow)
    FixationPointLoc_X = xCenterpixels;%Set the fixation point location to the center of the screen;
    FixationPointLoc_Y = yCenterpixels;%Set the fixation point location to the center of the screen;
    %{
else
    FixationPointLoc_X = displaywindow(3)/2;%Set the fixation point location to the center of the display window;
    FixationPointLoc_Y = displaywindow(4)/2;%Set the fixation point location to the center of the display window;
end    
%}
%Coordinate transform for object locations
%Transform cm to pixel;
PixPerCM_W = screenXpixels / ScreenWidthInCM;%How many cm does one pixel occupy in width;pix/cm
PixPerCM_H = screenYpixels / ScreenHeightInCM;%How many cm does one pixel occupy in height;pix/cm

%Transform object location from degree to cm;
ObjectEccentricityCM = tan(ObjectEccentricityDeg * pi / 180) * ViewingDistance;
%Transform object location representation from polor coordinate to cartesian coordinate;
ObjectLoc_XCM = ObjectEccentricityCM * cos(ObjectViewAngle * pi/180);
ObjectLoc_YCM = ObjectEccentricityCM * sin(ObjectViewAngle * pi/180);

%Transform cm to pixel;
ObjectLoc_X = FixationPointLoc_X + ObjectLoc_XCM * PixPerCM_W;
ObjectLoc_Y = FixationPointLoc_Y - ObjectLoc_YCM * PixPerCM_H;

%Transfrom degree to pixel
WholeViewingAngle_W = 2 * atan(ScreenWidthInCM / 2 / ViewingDistance) / pi * 180;
WholeViewingAngle_H = 2 * atan(ScreenHeightInCM / 2 / ViewingDistance) / pi * 180;

PixPerDeg_W = screenXpixels / WholeViewingAngle_W;
PixPerDeg_H = screenYpixels / WholeViewingAngle_H;


ObjectSize_W = ObjectSizeDeg(1) * PixPerDeg_W; %in pixel
ObjectSize_H = ObjectSizeDeg(2) * PixPerDeg_H; %in pixel

ObjectSize_WCM = ObjectSize_W / PixPerCM_W; %in cm
ObjectSize_HCM = ObjectSize_H / PixPerCM_H; %in cm



%Set up the initial location for the presentation of objects
baseRect = [0 0 ObjectSize_W ObjectSize_H];
centeredRect = CenterRectOnPointd(baseRect, ObjectLoc_X, ObjectLoc_Y);

%% ------------Basic set up for display window----------------------------%
% Open an on screen window and color the backgrounds black.
[window, ~] = Screen('OpenWindow', screenNumber, black, displaywindow);
% Enable alpha blending for anti-aliasing
Screen('BlendFunction', window, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
% Measure the vertical refresh rate of the monitor
ifi = Screen('GetFlipInterval', window);
% Retreive the maximum priority number
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

%% -------------Prepare the stimulus---------------------------------------%
%Prepare the object sets;
for i=1:NumberOfObjects
    ObjectSeq{i}=imread([ObjectPath ObjectName{i}]);
    texObjectSeq{i} = Screen('MakeTexture', window, ObjectSeq{i});
end 

%% ----------------Set up Object Motion parameters--------------------%
fr=Screen('NominalFrameRate', screenNumber);%Frame rate, in hz
SecPerFrame = 1 / fr; %in s
TotalMotionFrame = ceil(TotalMotionTime / 1000 / SecPerFrame);%total frame
ObjectMotionVeloctyInFrame_mm = ObjectMotionVelocty * TotalMotionTime / TotalMotionFrame;% mm/frame, unit change to preserve precision

%Calculate perspective change; Now assume screen center as zero for
%convenience
ObjCoor = [ObjectLoc_XCM * 1000, ObjectLoc_YCM * 1000, 0];%Set screen plane as zero; transfer to mm to preserve precision
EyeCoor = [0, 0, ViewingDistance * 1000];% transfer to mm
FPCoor = [0,0,0];

%First, translational component, projection of motion vector onto object-FP
%coordinate
if PureZooming
    ObjMoveVecUnit = [0,0,1]; %Unit Motion Vector
else
    ObjMoveVecUnit = (EyeCoor - ObjCoor)/norm((ObjCoor - EyeCoor)); %Unit Motion Vector
end
ObjMoveVecFull = ObjMoveVecUnit * ObjectMotionVeloctyInFrame_mm ; % Full Motion Vector

ObjMoveVecFull_Translation = dot(ObjMoveVecFull,(FPCoor - ObjCoor));
%Decompose into x and y 
ObjMoveVecFull_Translation_X = ObjMoveVecFull_Translation * cos(ObjectViewAngle * pi/180);
ObjMoveVecFull_Translation_Y = ObjMoveVecFull_Translation * sin(ObjectViewAngle * pi/180);

ObjMoveVecFull_X_Vec = linspace(0,ObjMoveVecFull_Translation_X,TotalMotionFrame) / 1000 * PixPerCM_W;%change into pixel
ObjMoveVecFull_Y_Vec = linspace(0,ObjMoveVecFull_Translation_Y,TotalMotionFrame) / 1000 * PixPerCM_H;%change into pixel

ObjectLoc_X_ChangeVec = ObjectLoc_X - ObjMoveVecFull_X_Vec;
ObjectLoc_Y_ChangeVec = ObjectLoc_Y + ObjMoveVecFull_Y_Vec;
%Shift to the real coordinate
ObjectLoc_X_ChangeVec_FP = ObjectLoc_X_ChangeVec;
ObjectLoc_Y_ChangeVec_FP = ObjectLoc_Y_ChangeVec;



%Second, expending/contraction component, projection of motion vector onto Eye-FP
%coordinate
ObjMoveVecFull_Expanding = dot(ObjMoveVecFull,(EyeCoor - FPCoor));
%The end vertical distance 
ObjMoveVecFull_ExpandingEnd = ViewingDistance * 1000 - ObjMoveVecFull_Expanding; %into mm

ObjMoveVecFull_ExpandingEnd_Deg_W = 2 * atan(ObjectSize_WCM / 2 * 1000 / ObjMoveVecFull_ExpandingEnd) * 180 / pi; %in deg
ObjMoveVecFull_ExpandingEnd_Deg_H = 2 * atan(ObjectSize_HCM / 2 * 1000 / ObjMoveVecFull_ExpandingEnd) * 180 / pi; %in deg

ObjMoveVecFull_ExpandingEnd_W = ObjMoveVecFull_ExpandingEnd_Deg_W * PixPerDeg_W; % in pixel
ObjMoveVecFull_ExpandingEnd_H = ObjMoveVecFull_ExpandingEnd_Deg_H * PixPerDeg_H; % in pixel

ObjMoveVecFull_ExpandingVec_W = linspace(ObjectSize_W,ObjMoveVecFull_ExpandingEnd_W, TotalMotionFrame);
ObjMoveVecFull_ExpandingVec_H = linspace(ObjectSize_H,ObjMoveVecFull_ExpandingEnd_H, TotalMotionFrame);

%Now make the base rectangle 
for i = 1 : TotalMotionFrame
    baseRectM = [0 0 ObjMoveVecFull_ExpandingVec_W(i) ObjMoveVecFull_ExpandingVec_H(i)];
    centeredRectM(i,:) = CenterRectOnPointd(baseRectM, ObjectLoc_X_ChangeVec_FP(i), ObjectLoc_Y_ChangeVec_FP(i));
end


%% ----------------Draw the stimulus in sequence--------------------------%
% First: Fixation Points
Screen('DrawDots', window, [FixationPointLoc_X FixationPointLoc_Y], FixationPointSizePix, FixationPointColor, [], 2);
Screen('Flip', window);
WaitSecs(1);%Wait for 1 second;

% Next: Stimulus
for i=1:NumberOfObjectsInOneTrial
    %Initial presentation
    Screen('DrawDots', window, [FixationPointLoc_X FixationPointLoc_Y], FixationPointSizePix, FixationPointColor, [], 2);
    Screen('DrawTexture', window, texObjectSeq{i}, [], centeredRect, 0);
    Screen('Flip', window);
    WaitSecs(0.5);
    
    %Start to zooming or shrinking
    for j = 1 : TotalMotionFrame
        WaitSecs(SecPerFrame/2);
        Screen('DrawDots', window, [FixationPointLoc_X FixationPointLoc_Y], FixationPointSizePix, FixationPointColor, [], 2);
        Screen('DrawTexture', window, texObjectSeq{i}, [], centeredRectM(j,:), 0);
        Screen('Flip', window);
        WaitSecs(SecPerFrame/2);  
    end
  
end
 WaitSecs(0.5);
Screen('DrawDots', window, [FixationPointLoc_X FixationPointLoc_Y], FixationPointSizePix, FixationPointColor, [], 2);
Screen('Flip', window);

%% ------------------End of the trial, close the window--------------------%
Priority(0);
sca
commandwindow

catch
    Priority(0);
    sca
    commandwindow
    % Output the error message that describes the error:
    psychrethrow(psychlasterror);
    
    keyboard
end




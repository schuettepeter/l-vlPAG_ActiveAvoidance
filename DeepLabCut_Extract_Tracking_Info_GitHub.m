%% Convert DeepLabCut output (.h5 file) to Tracking file with angle, position, vel, accel, etc. as well as behavior categorization

clearvars -except RTPP folderzzz framesToDelete truncateTracking lastFrame hotplate_xp_CNN moving_shockgrid open_field folders folderzz jjj CO2 deleteUserDefinedZone setByParentScript miniscope other_vid fear_cond_chamber opto_xp fear_cond_chamber ...
    simpleRat_xp ratClimbUp_xp hotplate_xp moving_shockgrid mouse_toyrat burrow_xp burrow_norat_xp shockgrid_xp EPM toyratClimbUp_xp

setByParentScript = 0; %select '1' if xp type chosen by parent "Workflow" script.

isHD = 0;
dontUseH5 = 0; %just use existing Tracking.mat file

% TO READ IN h5 FILE:

%Tracking smooth window width
windowWidthSmooth = 1;

truncateTracking = 0;
lastFrame = 1600000;

framesToDelete = [];
miniscope = 1; 
other_vid = 0; 

%select assay here:
simpleRat_xp = 1;
mouse_toyrat = 0;
burrow_xp = 0;
burrow_norat_xp = 0;
shockgrid_xp = 0;
EPM = 0;
open_field = 0;
RTPP = 0;

cropOutsideRearingConstraint = 0;
deleteUserDefinedZone = 0; 

if mouse_toyrat == 1 | shockgrid_xp == 1 | fear_cond_chamber == 1 | moving_shockgrid == 1 | EPM == 1 | burrow_norat_xp == 1 | toyratClimbUp_xp == 1 | CO2 == 1 | hotplate_xp == 1 | hotplate_xp_CNN == 1 |open_field == 1 | RTPP == 1 | CO2_miniscope == 1 | cricket == 1;
    mouse_only = 1; mouse_rat = 0;
elseif simpleRat_xp == 1 | burrow_xp == 1 | ratClimbUp_xp == 1
    mouse_only = 0; mouse_rat = 1;
elseif artificial_prey==1
    mouse_only = 0; mouse_rat = 0;    
end

crop_tracking = 0;
area_occupancy = 1; %determine how much time mouse occupies user-selected area
flipRatSide = 0; %1 if threat (any threat) on left side, 0 if threat on right side

basedir = pwd;
data = [];
data_grid = [];

    S = dir('*.h5');  %dir('*354000.h5');
    dataAdd = h5read(S.name,'/df_with_missing/table');
    dataAddMat = dataAdd.values_block_0(:,:);
    data = [data, dataAddMat];

%Create video object
        if exist('behav.mp4', 'file')
            obj = VideoReader('behav.mp4');
            vidname = 'behav.mp4';
        elseif exist('behav.avi', 'file')
            obj = VideoReader('behav.avi');
            vidname = 'behav.avi';
        elseif exist('behavCam.avi', 'file')
            obj = VideoReader('behavCam.avi');
            vidname = 'behavCam.avi';
        elseif exist('saline.avi','file')
            obj = VideoReader('saline.avi');
            vidname = 'behav.avi';
        elseif exist('cno.avi','file')
            obj = VideoReader('cno.avi');
            vidname = 'behav.avi';
        elseif exist('behav.mov','file')
            obj = VideoReader('behav.mov');
            vidname = 'behav.mov';
        else
            temp = dir('*.avi'); temp = temp.name;
            obj = VideoReader(temp);
            vidname = temp; clearvars temp
        end
    
if mouse_toyrat == 1
%read in toy rat location
if ~exist([basedir '/toyRat_coords.mat'], 'file')
        frame = read(obj,200);
        image(frame)
        disp('choose location of toy rat')
        toyRat_coords = ginput(1);
        toyRat_coords = round(toyRat_coords);
        save([basedir '/toyRat_coords'], 'toyRat_coords') 
else
    load([basedir, '/toyRat_coords.mat'])
end 
end

%read in bounds for rearing classification
if EPM == 0 && ratClimbUp_xp == 0 && toyratClimbUp_xp == 0 && hotplate_xp_CNN == 0 && CO2_miniscope == 0 && fear_cond_chamber == 0 && open_field==0
if ~exist([basedir '/rearing_constraint.mat'], 'file') & ~exist([basedir '/bounds.mat'], 'file')
        %if miniscope == 0
            frame = read(obj,200);
            image(frame); axis image
        %elseif miniscope == 1
        %    frame = double(msReadFrame(obj,i,false,false,false))/255;
        %    image(frame); axis image
        %end
        disp('choose bounds for rearing classification')
        rearing_constraint = ginput(2);
        rearing_constraint = round(rearing_constraint);
        save([basedir '/rearing_constraint'], 'rearing_constraint') 
else
    if exist([basedir '/rearing_constraint.mat'], 'file')
        load([basedir, '/rearing_constraint.mat'])
    end
    if exist([basedir '/bounds.mat'], 'file')
        load([basedir, '/bounds.mat'])
        %rearing_constraint = Bounds.hose;
        rearing_constraint = Bounds;
    end
end
end

%read in bounds for elevated plus maze closed and open walls
if EPM == 1
if ~exist([basedir '/PlusArmN.mat'], 'file')
        if miniscope == 0 | miniscope == 1
        frame = read(obj,200);
        image(frame); axis image;
        end
        
        disp('choose bounds for North Arm')
        PlusArmN = ginput(2);
        PlusArmN = round(PlusArmN);
        save([basedir '/PlusArmN'], 'PlusArmN') 
        
        image(frame)
        disp('choose bounds for East Arm')
        PlusArmE = ginput(2);
        PlusArmE = round(PlusArmE);
        save([basedir '/PlusArmE'], 'PlusArmE') 
        
        image(frame)
        disp('choose bounds for South Arm')
        PlusArmS = ginput(2);
        PlusArmS = round(PlusArmS);
        save([basedir '/PlusArmS'], 'PlusArmS') 
        
        image(frame)
        disp('choose bounds for West Arm')
        PlusArmW = ginput(2);
        PlusArmW = round(PlusArmW);
        save([basedir '/PlusArmW'], 'PlusArmW') 

else
    load([basedir, '/PlusArmN.mat']); load([basedir, '/PlusArmE.mat']); load([basedir, '/PlusArmS.mat']); load([basedir, '/PlusArmW.mat']);
end
end

%read in bounds for head dips
if EPM == 1
if ~exist([basedir '/PlusArmNW.mat'], 'file')
        if miniscope == 0 | miniscope == 1
        frame = read(obj,200);
        image(frame); axis image;
        end
        %if miniscope == 1
        %frame = double(msReadFrame(obj,i,false,false,false))/255;
        %image(frame); axis image;
        %end
        
        disp('choose bounds for North East Arm')
        PlusArmNE = ginput(2);
        PlusArmNE = round(PlusArmNE);
        save([basedir '/PlusArmNE'], 'PlusArmNE') 
        
        image(frame)
        disp('choose bounds for South East Arm')
        PlusArmSE = ginput(2);
        PlusArmSE = round(PlusArmSE);
        save([basedir '/PlusArmSE'], 'PlusArmSE') 
        
        image(frame)
        disp('choose bounds for South West Arm')
        PlusArmSW = ginput(2);
        PlusArmSW = round(PlusArmSW);
        save([basedir '/PlusArmSW'], 'PlusArmSW') 
        
        image(frame)
        disp('choose bounds for North West Arm')
        PlusArmNW = ginput(2);
        PlusArmNW = round(PlusArmNW);
        save([basedir '/PlusArmNW'], 'PlusArmNW') 

else
    load([basedir, '/PlusArmNE.mat']); load([basedir, '/PlusArmNW.mat']); load([basedir, '/PlusArmSE.mat']); load([basedir, '/PlusArmSW.mat']);
end
end

if open_field == 1
    if ~exist('of_center.mat', 'file')
        reader = VideoReader(vidname);
        frame = readFrame(reader);
        image(frame)
        text(20,30, 'choose open field center box (top left and bottom right corner)', 'Color', 'green')
        of_center = ginput(2);
        of_center = round(of_center);
        save('of_center', 'of_center');
    else
            load('of_center.mat')
    end
end

        frame = read(obj,100);
        image(frame); axis image; 
        
    if (burrow_xp == 1 | burrow_norat_xp == 1) & ~exist('burrow_boundary.mat', 'file')
        disp('choose bounds for burrow classification')
        burrow_boundary = ginput(2);
        burrow_boundary = round(burrow_boundary);
        save([basedir '/burrow_boundary'], 'burrow_boundary') 
        
    elseif (burrow_xp == 1 | burrow_norat_xp == 1) & exist('burrow_boundary.mat', 'file')
        load([basedir, '/burrow_boundary.mat'])  
        
    elseif shockgrid_xp == 1 & ~exist('shockgrid_boundary.mat', 'file')
        disp('choose bounds for shockgrid location')
        shockgrid_boundary = ginput(1);
        shockgrid_boundary = round(shockgrid_boundary);
        save([basedir '/shockgrid_boundary'], 'shockgrid_boundary') 
        
    elseif shockgrid_xp == 1 & exist('shockgrid_boundary.mat', 'file')
        load([basedir, '/shockgrid_boundary.mat'])  
    end
    
if burrow_xp == 100
%read in bounds for non-burrow corner
if ~exist([basedir '/nonburrow_boundary.mat'], 'file')
        obj = VideoReader('rat.avi');
        frame = read(obj,100);
        image(frame)
        disp('choose bounds for nonburrow corner classification')
        nonburrow_boundary = ginput(2);
        nonburrow_boundary = round(nonburrow_boundary);
        save([basedir '/nonburrow_boundary'], 'nonburrow_boundary') 
else
    load([basedir, '/nonburrow_boundary.mat'])  
end
end

if burrow_xp == 1 | burrow_norat_xp == 1
%calculate burrow center for distance measures
burrow_center = [(burrow_boundary(1,1) + burrow_boundary(2,1)) / 2, (burrow_boundary(1,2) + burrow_boundary(2,2)) ./ 2];
end

if EPM == 0 && ratClimbUp_xp == 0 && toyratClimbUp_xp == 0 && hotplate_xp_CNN == 0 && RTPP == 0 && hotplate_xp == 0 && CO2 == 0 && fear_cond_chamber == 0 && open_field==0
%calculate corner coordinates
corners.NW_corner = [rearing_constraint(1,1), rearing_constraint(1,2)];
corners.NE_corner = [rearing_constraint(2,1), rearing_constraint(1,2)];
corners.SW_corner = [rearing_constraint(1,1), rearing_constraint(2,2)];
corners.SE_corner = [rearing_constraint(2,1), rearing_constraint(2,2)];
end

if exist('Bounds.mat', 'file')
    load('Bounds.mat')
end

% then convert to struct:
if CO2 == 0 && hotplate_xp==0
    
if cricket == 0 && artificial_prey==0
        Tracking.mouseNose = hampel(data(1:2, :))';
        Tracking.mouseLeftEar = hampel(data(4:5, :))';
        Tracking.mouseRightEar = hampel(data(7:8, :))';
        Tracking.mouseTailbase = hampel(data(10:11, :))';
end

if artificial_prey == 1
        Tracking.artificialPrey = hampel(data(1:2, :))';
end

if cricket == 1
        Tracking.mouseNose = hampel(data(4:5, :))';
        Tracking.mouseLeftEar = hampel(data(7:8, :))';
        Tracking.mouseRightEar = hampel(data(10:11, :))';
        Tracking.mouseTailbase = hampel(data(13:14, :))';
        Tracking.cricket = hampel(data(1:2, :))';
        
        if exist('cricketInFrame')
            Tracking.cricket(~cricketInFrame,:) = nan;
        end
end

if mouse_rat == 1;
    Tracking.ratNose = hampel(data(13:14, :))';
    Tracking.ratLeftNose = hampel(data(16:17, :))';
    Tracking.ratRightNose = hampel(data(19:20, :))';
    Tracking.ratTailbase = hampel(data(22:23, :))';
end

% and remove frames if necessary:
Tracking.mouseNose(framesToDelete,:) = nan;
Tracking.mouseLeftEar(framesToDelete,:) = nan;
Tracking.mouseRightEar(framesToDelete,:) = nan;
Tracking.mouseTailbase(framesToDelete,:) = nan;

if mouse_rat == 1;
Tracking.ratNose(framesToDelete,:) = nan;
Tracking.ratLeftNose(framesToDelete,:) = nan;
Tracking.ratRightNose(framesToDelete,:) = nan;
Tracking.ratTailbase(framesToDelete,:) = nan;
end

% change points to nans that fall outside y bounds
if crop_tracking == 1
offset = 0;

if ratClimbUp_xp==0 & toyratClimbUp_xp==0
if exist('bounds.mat', 'file')
    load('bounds.mat')
    %rearing_constraint = Bounds.enclosure;
    rearing_constraint = bounds;
end
end
end
end

%% and add important metrics

if CO2 == 0 && hotplate_xp==0 && artificial_prey==0
%Angle of head direction -- first find point between ears -- then find
%angle made between that point and nose

for i = 1:length(Tracking.mouseNose)
    Tracking.mouseBetwEars(i,:) = (Tracking.mouseLeftEar(i,:) + Tracking.mouseRightEar(i,:)).'/2;
end
Tracking.mouse_position = Tracking.mouseBetwEars;

% center all points of mouse
temp_x = nanmean([Tracking.mouseLeftEar(:,1), Tracking.mouseRightEar(:,1), Tracking.mouseTailbase(:,1)],2);
temp_y = nanmean([Tracking.mouseLeftEar(:,2), Tracking.mouseRightEar(:,2), Tracking.mouseTailbase(:,2)],2);

Tracking.centerMouse = [temp_x, temp_y];

if mouse_rat == 1
for i = 1:length(Tracking.ratNose)
    Tracking.ratBetwEars(i,:) = (Tracking.ratLeftNose(i,:) + Tracking.ratRightNose(i,:)).'/2;
end
Tracking.rat_position = Tracking.ratBetwEars;
end

for i = 1:length(Tracking.mouseNose)
    Tracking.mouseAngle(i) = atan2(Tracking.mouseBetwEars(i,2) - Tracking.mouseNose(i,2), Tracking.mouseBetwEars(i,1) - Tracking.mouseNose(i,1));
end
Tracking.mouseAngle = Tracking.mouseAngle';

if mouse_rat == 1
for i = 1:length(Tracking.ratNose)
    Tracking.ratAngle(i) = atan2(Tracking.ratBetwEars(i,2) - Tracking.ratNose(i,2), Tracking.ratBetwEars(i,1) - Tracking.ratNose(i,1));
end
Tracking.ratAngle = Tracking.ratAngle';

%find tracking between direction of mouse head and rat, between ears
for i = 1:length(Tracking.mouseBetwEars)
    Tracking.angleMouseBetwEarsRatBetwEars(i,:) = atan2(Tracking.mouseBetwEars(i,2) - Tracking.ratBetwEars(i,2), Tracking.mouseBetwEars(i,1) - Tracking.ratBetwEars(i,1));
end

%find difference of mouse head direction and angle between mouse and rat
Tracking.angleDiffMouseHeadDirRat = abs(abs(Tracking.mouseAngle) - abs(Tracking.angleMouseBetwEarsRatBetwEars));

end

if mouse_toyrat==1

for i = 1:length(Tracking.mouse_position)
    Tracking.ratAngle(i) = atan2(toyRat_coords(1,2) - toyRat_coords(1,2), toyRat_coords(1,1) - (toyRat_coords(1,1)-10));
end
Tracking.ratAngle = Tracking.ratAngle';

%find tracking between direction of mouse head and rat, between ears
for i = 1:length(Tracking.mouseBetwEars)
    Tracking.angleMouseBetwEarsRatBetwEars(i,:) = atan2(Tracking.mouseBetwEars(i,2) - toyRat_coords(1,2), Tracking.mouseBetwEars(i,1) - toyRat_coords(1,1));
end

%find difference of mouse head direction and angle between mouse and rat
Tracking.angleDiffMouseHeadDirRat = abs(abs(Tracking.mouseAngle) - abs(Tracking.angleMouseBetwEarsRatBetwEars));

%find points on toy rat and make angle zero.
tempNose = Tracking.mouseNose(:,1) > toyRat_coords(1)-5 &  Tracking.mouseNose(:,1) < toyRat_coords(1)+90 & ...
    Tracking.mouseNose(:,2) > toyRat_coords(2)-20 &  Tracking.mouseNose(:,2) < toyRat_coords(2)+20;
tempTail = Tracking.mouseTailbase(:,1) > toyRat_coords(1)-5 &  Tracking.mouseTailbase(:,1) < toyRat_coords(1)+90 & ...
    Tracking.mouseTailbase(:,2) > toyRat_coords(2)-20 &  Tracking.mouseTailbase(:,2) < toyRat_coords(2)+20;
temp = tempNose & tempTail;
temp = find(temp);
Tracking.angleDiffMouseHeadDirRat(temp) = 0;
end

%and take cosine and sine of rat angles
Tracking.sinMouseAngle = sin(Tracking.mouseAngle); 
Tracking.cosMouseAngle = cos(Tracking.mouseAngle); 

if mouse_rat == 1
sinRatAngle = sin(Tracking.ratAngle)';
cosRatAngle = cos(Tracking.ratAngle)';
Tracking.sinMouseRatAngle = sin(Tracking.mouseAngle - Tracking.ratAngle)'; 
Tracking.cosMouseRatAngle = cos(Tracking.mouseAngle-Tracking.ratAngle)';
end


%% Positions of rat mouse and related
%mouse velocity and acceleration
w = gausswin(5); %for smoothing motion parameters

Tracking.directionalVelocity = diff(Tracking.mouseBetwEars);
Tracking.directionalVelocity = [Tracking.directionalVelocity(1,:); Tracking.directionalVelocity];

positiondiffMouse = diff(Tracking.mouseBetwEars);
positiondiffMouse = [positiondiffMouse(1,:); positiondiffMouse]; %keep same length

for i = 1:size(positiondiffMouse, 1)
   Tracking.mouseVel(i) = sqrt(positiondiffMouse(i,1)^2 + positiondiffMouse(i,2)^2); %find hypotenuse
end
Tracking.mouseVel = Tracking.mouseVel';

Tracking.mouseAccel = diff(Tracking.mouseVel);
Tracking.mouseAccel = [Tracking.mouseAccel(1); Tracking.mouseAccel];

%velocity mouse tailbase
positiondiffMouseTailbase = diff(Tracking.mouseTailbase);
positiondiffMouseTailbase = [positiondiffMouseTailbase(1,:); positiondiffMouseTailbase]; %keep same length

%change in angle mouse head
Tracking.positiondiffMouseAngle = diff(Tracking.mouseAngle);
Tracking.positiondiffMouseAngle = [Tracking.positiondiffMouseAngle(1); Tracking.positiondiffMouseAngle]; %keep same length

%tailbase velocity of mouse
for i = 1:size(positiondiffMouseTailbase, 1)
   Tracking.mouseTailbaseVel(i) = sqrt(positiondiffMouseTailbase(i,1)^2 + positiondiffMouseTailbase(i,2)^2); %find hypotenuse
end
Tracking.mouseTailbaseVel = Tracking.mouseTailbaseVel';

%nose velocity of mouse
positiondiffMouseNose = diff(Tracking.mouseNose);
positiondiffMouseNose = [positiondiffMouseNose(1,:); positiondiffMouseNose]; %keep same length

for i = 1:size(positiondiffMouseNose, 1)
   Tracking.mouseNoseVel(i) = sqrt(positiondiffMouseNose(i,1)^2 + positiondiffMouseNose(i,2)^2); %find hypotenuse
end
Tracking.mouseNoseVel = Tracking.mouseNoseVel';

if mouse_rat == 1
%rat
positiondiffRat = filter(w,1,diff(Tracking.ratBetwEars));
positiondiffRat = [positiondiffRat(1,:); positiondiffRat]; %keep same length

for i = 1:size(positiondiffRat, 1)
   Tracking.ratVel(i) = sqrt(positiondiffRat(i,1)^2 + positiondiffRat(i,2)^2); %find hypotenuse
end
Tracking.ratVel = Tracking.ratVel';

Tracking.ratAccel = diff(Tracking.ratVel);
Tracking.ratAccel = [Tracking.ratAccel(1); Tracking.ratAccel];
end

if cricket == 1
%cricket
positiondiffCricket = filter(w,1,diff(Tracking.cricket));
positiondiffCricket = [positiondiffCricket(1,:); positiondiffCricket]; %keep same length

for i = 1:size(positiondiffCricket, 1)
   Tracking.cricketVel(i) = sqrt(positiondiffCricket(i,1)^2 + positiondiffCricket(i,2)^2); %find hypotenuse
end
Tracking.cricketVel = Tracking.cricketVel';

Tracking.cricketAccel = diff(Tracking.cricketVel);
Tracking.cricketAccel = [Tracking.cricketAccel(1); Tracking.cricketAccel];
end

%x coord position of mouse and rat
Tracking.mouseXPosition = Tracking.mouseBetwEars(:,1);
if mouse_rat == 1
Tracking.ratXPosition = Tracking.ratBetwEars(:,1);
end
if cricket == 1
Tracking.cricketXPosition = Tracking.cricket(:,1);
end

if mouse_rat == 1
%distance between mouse and rat at each frame
for i = 1:length(Tracking.mouseBetwEars)
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2);Tracking.ratBetwEars(i,1), Tracking.ratBetwEars(i,2)];
    Tracking.distanceMouseRat(i) = pdist(temp,'euclidean');
end
Tracking.distanceMouseRat = Tracking.distanceMouseRat';

%change in distance between mouse and rat at each frame
diffDistanceMouseRat = diff(Tracking.distanceMouseRat);
Tracking.diffDistanceMouseRat = [diffDistanceMouseRat(1); diffDistanceMouseRat];
end

if mouse_toyrat == 1 | toyratClimbUp_xp == 1 %toyRat_coords
%distance between mouse and rat at each frame
for i = 1:length(Tracking.mouseBetwEars)
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2); toyRat_coords(1), toyRat_coords(2)];
    Tracking.distanceMouseToyRat(i) = pdist(temp,'euclidean');
end
Tracking.distanceMouseToyRat = Tracking.distanceMouseToyRat';

%change in distance between mouse and rat at each frame
diffDistanceMouseToyRat = diff(Tracking.distanceMouseToyRat);
Tracking.diffDistanceMouseToyRat = [diffDistanceMouseToyRat(1); diffDistanceMouseToyRat];
end

if shockgrid_xp == 1
    
%distance between mouse and shockgrid boundary
for i = 1:length(Tracking.mouseBetwEars)
    Tracking.distanceShockgrid(i) = abs(shockgrid_boundary(1) - Tracking.mouseBetwEars(i,1));
end
Tracking.distanceShockgrid = Tracking.distanceShockgrid';
    
%change in distance between mouse and shockgrid at each frame
diffDistanceShockgrid = diff(Tracking.distanceShockgrid);
Tracking.diffDistanceShockgrid = [diffDistanceShockgrid(1); diffDistanceShockgrid];
end

%distance between mouse nose and tailbase at each frame (i.e. for
%calculating stretch)
for i = 1:size(Tracking.mouseNose, 1)
   Tracking.distNoseTailbase(i) = sqrt((Tracking.mouseNose(i,1) - Tracking.mouseTailbase(i,1))^2 + (Tracking.mouseNose(i,2) - Tracking.mouseTailbase(i,2))^2); %find hypotenuse
end
Tracking.distNoseTailbase = Tracking.distNoseTailbase';
%idx2delStr = find(Tracking.distNoseTailbase > 90);
%Tracking.distNoseTailbase(idx2delStr) = nan;

%and for distance from between ears to tailbase
for i = 1:size(Tracking.mouseNose, 1)
   Tracking.distEarsTailbase(i) = sqrt((Tracking.mouse_position(i,1) - Tracking.mouseTailbase(i,1))^2 + (Tracking.mouse_position(i,2) - Tracking.mouseTailbase(i,2))^2); %find hypotenuse
end
Tracking.distEarsTailbase = Tracking.distEarsTailbase';
idx2delStr = find(Tracking.distEarsTailbase > 90);
Tracking.distEarsTailbase(idx2delStr) = nan;


%then calculate zscore --seems like possible cutoff for stretch is 1.8
stdDistNoseTailbase = nanstd(Tracking.distNoseTailbase);
meanDistNoseTailbase = nanmean(Tracking.distNoseTailbase);
Tracking.distNoseTailBaseZ = (Tracking.distNoseTailbase - meanDistNoseTailbase) ./ stdDistNoseTailbase;

%calculate distance form center of burrow for each frame
if burrow_xp == 1 | burrow_norat_xp == 1
for i = 1:length(Tracking.mouseBetwEars)
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2);burrow_center(1,1), burrow_center(1,2)];
    Tracking.distanceFromBurrow(i) = pdist(temp,'euclidean');
end
Tracking.distanceFromBurrow = Tracking.distanceFromBurrow';
end

%calculate distance from corners of enclosure
if burrow_xp == 1 | burrow_norat_xp == 1
for i = 1:length(Tracking.mouseBetwEars)
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2);corners.NW_corner(1,1), corners.NW_corner(1,2)];
    Tracking.distance_NW_Corner(i) = (pdist(temp,'euclidean'))';
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2);corners.NE_corner(1,1), corners.NE_corner(1,2)];
    Tracking.distance_NE_Corner(i) = (pdist(temp,'euclidean'))';
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2);corners.SW_corner(1,1), corners.SW_corner(1,2)];
    Tracking.distance_SW_Corner(i) = (pdist(temp,'euclidean'))';
    temp = [Tracking.mouseBetwEars(i,1), Tracking.mouseBetwEars(i,2);corners.SE_corner(1,1), corners.SE_corner(1,2)];
    Tracking.distance_SE_Corner(i) = (pdist(temp,'euclidean'))';
end
Tracking.distance_NW_Corner = Tracking.distance_NW_Corner'; Tracking.distance_NE_Corner = Tracking.distance_NE_Corner';
Tracking.distance_SW_Corner = Tracking.distance_SW_Corner'; Tracking.distance_SE_Corner = Tracking.distance_SE_Corner';
end
end

%% Automatic behavior scoring

if other_vid == 1
    sampleRate = 30; %was 15 fps
elseif miniscope == 1
    sampleRate = 30; %fps
end

if exist('fps', 'file')
    load('fps.mat')
    sampleRate = fps;
end

%% stretch -- find indices where distance from nose to tailbase is greater than threshold and tailbase stationary
if CO2 == 0 && hotplate_xp == 0 && artificial_prey==0
    
if EPM == 1 | RTPP == 1 | open_field==1
    strThresh = 54; 
    strThreshEars = 51; 
    strMax = 120; 

elseif EPM == 0 && ratClimbUp_xp == 0 && toyratClimbUp_xp == 0
    strThresh = 69; 
    strThreshEars = 65; 
    strMax = 100; 
    if isHD==1
        strThresh = 75;
        strThreshEars = 70;
        strMax = 130;
    end
    
elseif ratClimbUp_xp == 1 | toyratClimbUp_xp == 1
    strThresh = 106;
end
    
tailbaseThreshVel = 1; %was 1 -- max movement allowed for stretch tailbase

if isHD==1
    tailbaseThreshVel = 2; %was 1 -- max movement allowed for stretch tailbase
end    

if miniscope == 0
stretchIndices = Tracking.distNoseTailbase > strThresh & Tracking.distNoseTailbase < strMax & Tracking.mouseTailbaseVel < tailbaseThreshVel;
elseif miniscope == 1 && artificial_prey==0
    stretchIndices = Tracking.distEarsTailbase > strThreshEars & Tracking.distNoseTailbase < strMax & Tracking.mouseTailbaseVel < tailbaseThreshVel;
end

windowWidth = round(sampleRate .* .75);
counts = conv(stretchIndices, ones(1, windowWidth), 'same');
stretchIndices = counts >= round(sampleRate * .5);

stretchSwitch = diff(stretchIndices);
stretchSwitchIndices = find(stretchSwitch);
    
jj = 1; kk = 1;
for i=1:length(stretchSwitchIndices)
    if stretchSwitch(stretchSwitchIndices(i)) == 1
        stretchStart(jj) = stretchSwitchIndices(i);
        jj = jj+1;
    elseif stretchSwitch(stretchSwitchIndices(i)) == -1
        stretchEnd(kk) = stretchSwitchIndices(i);
        kk = kk+1;
    end
end

%remove indices if start & end points not complete bookends (acc by start or end)
if exist('stretchStart')
if stretchStart(1) > stretchEnd(1)
    stretchIndices(1:stretchStart(1)) = 0;
end
if stretchStart(end) > stretchEnd(end)
    stretchIndices(stretchStart(end):length(stretchIndices)) = 0;
end

%update start and end indices
 clearvars stretchSwitch stretchSwitchIndices stretchStart stretchEnd;
% 
 stretchSwitch = diff(stretchIndices);
 stretchSwitchIndices = find(stretchSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(stretchSwitchIndices)
     if stretchSwitch(stretchSwitchIndices(i)) == 1
         stretchStart(jj) = stretchSwitchIndices(i);
         jj = jj+1;
     elseif stretchSwitch(stretchSwitchIndices(i)) == -1
         stretchEnd(kk) = stretchSwitchIndices(i);
         kk = kk+1;
     end
 end

 for i = 1:length(stretchStart)
     if stretchEnd(i) - stretchStart(i) < round(sampleRate .* .3) %was .5 -- set the minimum duration of stretch here
         stretchIndices(stretchStart(i):stretchEnd(i)) = 0;
     end
 end 
 %make sure no stretches counted on platform
 if ratClimbUp_xp == 1 | toyratClimbUp_xp == 1
     for i = 1:length(stretchStart)
     if Tracking.mouse_position(stretchStart(i),2) < 350  %set the lowest y value to be considered a stretch (don't count if on platform)
         stretchIndices(stretchStart(i):stretchEnd(i)) = 0;
     end
 end   
 end
 
%update start and end indices
 clearvars stretchSwitch stretchSwitchIndices stretchStart stretchEnd;
% 
 stretchSwitch = diff(stretchIndices);
 stretchSwitchIndices = find(stretchSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(stretchSwitchIndices)
     if stretchSwitch(stretchSwitchIndices(i)) == 1
         stretchStart(jj) = stretchSwitchIndices(i);
         jj = jj+1;
     elseif stretchSwitch(stretchSwitchIndices(i)) == -1
         stretchEnd(kk) = stretchSwitchIndices(i);
         kk = kk+1;
     end
 end
end
end

%% escape
if mouse_rat == 1 | shockgrid_xp == 1 | mouse_toyrat == 1 | ratClimbUp_xp == 1 | toyratClimbUp_xp == 1 | moving_shockgrid == 1 | cricket == 1 | artificial_prey == 1

if flipRatSide == 0 && ratClimbUp_xp == 0 && toyratClimbUp_xp == 0 && moving_shockgrid == 0
if mouse_toyrat == 1 
    aboveEscThreshold = Tracking.diffDistanceMouseToyRat(:,1) > 2; %pixels per frame
elseif mouse_rat == 1
    aboveEscThreshold = Tracking.directionalVelocity(:,1) < -2; %pixels per frame
elseif shockgrid_xp == 1
    aboveEscThreshold = Tracking.diffDistanceShockgrid(:,1) > 2 & Tracking.mouse_position(:,1) < shockgrid_boundary(1); %pixels per frame
end
end

if flipRatSide == 0 && ratClimbUp_xp == 1 | toyratClimbUp_xp == 1
if toyratClimbUp_xp == 1 
    aboveEscThreshold = Tracking.diffDistanceMouseToyRat(:,1) > .5 & Tracking.mouse_position(:,2) > Bounds.lowerTrack(1,2); %pixels per frame
elseif ratClimbUp_xp == 1
    aboveEscThreshold = Tracking.directionalVelocity(:,1) > 2 & Tracking.mouse_position(:,2) > Bounds.lowerTrack(1,2); %pixels per frame
end
end

if flipRatSide == 1
if mouse_toyrat == 1 
    aboveEscThreshold = Tracking.diffDistanceMouseToyRat(:,1) > 2; %pixels per frame
elseif mouse_rat == 1
    aboveEscThreshold = Tracking.directionalVelocity(:,1) > 2; %pixels per frame
elseif shockgrid_xp == 1
    aboveEscThreshold = Tracking.diffDistanceShockgrid(:) > 2 & Tracking.mouse_position(:,1) > shockgrid_boundary(1); %pixels per frame
end
end

% Make a moving window to count the number of frames minimum velocity.
windowWidth = round(sampleRate .* .7); %was * 1.5
counts = conv(aboveEscThreshold, ones(1, windowWidth), 'same');
escapeIndices = counts >= round(sampleRate * .4); %was .7

escapeSwitch = diff(escapeIndices);
escapeSwitchIndices = find(escapeSwitch);

if escapeIndices(length(escapeIndices)) == 1
    escapeIndices(escapeSwitchIndices(end):length(escapeIndices)) = 0;
elseif escapeIndices(1) == 1
    escapeIndices(1:escapeSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(escapeSwitchIndices)
    if escapeSwitch(escapeSwitchIndices(i)) == 1
        escapeStart(jj) = escapeSwitchIndices(i);
        jj = jj+1;
    elseif escapeSwitch(escapeSwitchIndices(i)) == -1
        escapeEnd(kk) = escapeSwitchIndices(i);
        kk = kk+1;
    end
end

%remove indices if start & end points not complete bookends (acc by start or end)
if exist('escapeStart')
if escapeStart(1) > escapeEnd(1)
    escapeIndices(1:escapeStart(1)) = 0;
end
if escapeStart(end) > escapeEnd(end)
    escapeIndices(escapeStart(end):length(escapeIndices)) = 0;
end

% %update start and end indices
 clearvars escapeSwitch escapeSwitchIndices escapeStart escapeEnd;
% 
 escapeSwitch = diff(escapeIndices);
 escapeSwitchIndices = find(escapeSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(escapeSwitchIndices)
     if escapeSwitch(escapeSwitchIndices(i)) == 1
         escapeStart(jj) = escapeSwitchIndices(i);
         jj = jj+1;
     elseif escapeSwitch(escapeSwitchIndices(i)) == -1
         escapeEnd(kk) = escapeSwitchIndices(i);
         kk = kk+1;
     end
 end
 
%adjust initial distance from rat
if mouse_rat == 1
    if exist('escapeStart')
 for i = 1:length(escapeStart)
     if Tracking.distanceMouseRat(escapeStart(i)) > 250 %was 425
         escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
     end
 end
 
 %delete instance if mouse doesn't escape rat by a certain amount
  for i = 1:length(escapeStart)
     if Tracking.distanceMouseRat(escapeEnd(i)) - Tracking.distanceMouseRat(escapeStart(i)) < 100 
         escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
     end
  end
    end
end

%adjust initial distance from  toy rat
if mouse_toyrat == 1
    if exist('escapeStart')
 for i = 1:length(escapeStart)
     if Tracking.distanceMouseToyRat(escapeStart(i)) > 250 %was 425
         escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
     end
 end
 
 %delete instance if mouse doesn't escape rat by a certain amount
  for i = 1:length(escapeStart)
     if Tracking.distanceMouseToyRat(escapeEnd(i)) - Tracking.distanceMouseToyRat(escapeStart(i)) < 100 
         escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
     end
  end
    end
end

if shockgrid_xp == 1
    if exist('escapeStart')
        %initial minimum distance of escape
         for i = 1:length(escapeStart)
             if abs(Tracking.distanceShockgrid(escapeStart(i))) > 250 %was 250
                 escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
             end
         end
         %distance run during escape
    for i = 1:length(escapeStart)
            if abs(Tracking.distanceShockgrid(escapeEnd(i)) - Tracking.distanceShockgrid(escapeStart(i))) < 100 
                escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
            end
    end

    end
end

%adjust minimum escape duration
    if length(escapeStart) > length(escapeEnd)
        escapeStart = escapeStart(1:end-1)
    end
 for i = 1:length(escapeStart)
     if escapeEnd(i) - escapeStart(i) < round(sampleRate .* .5) %set the minimum duration here
         escapeIndices(escapeStart(i):escapeEnd(i)) = 0;
     end
 end 

% %update start and end indices
 clearvars escapeSwitch escapeSwitchIndices escapeStart escapeEnd;
% 
 escapeSwitch = diff(escapeIndices);
 escapeSwitchIndices = find(escapeSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(escapeSwitchIndices)
     if escapeSwitch(escapeSwitchIndices(i)) == 1
         escapeStart(jj) = escapeSwitchIndices(i);
         jj = jj+1;
     elseif escapeSwitch(escapeSwitchIndices(i)) == -1
         escapeEnd(kk) = escapeSwitchIndices(i);
         kk = kk+1;
     end
 end

end
end

 
%% approach

if mouse_rat == 1 | shockgrid_xp == 1 | mouse_toyrat == 1 | ratClimbUp_xp == 1 | toyratClimbUp_xp == 1 | moving_shockgrid == 1 | cricket == 1 | artificial_prey==1

    minPixCloserApp = 100;  
    if isHD==1
    minPixCloserApp = 300;   
    end

if flipRatSide == 0
    appDistThresh = 200; %set minimum approach end distance to count as approach towards rat (must pass a certain threshold)
elseif flipRatSide == 1
    appDistThresh = 200; %set minimum approach end distance to count as approach towards rat (must pass a certain threshold)
end
    
minVelApp = -2; 

if flipRatSide == 0
if mouse_rat == 1 && ratClimbUp_xp == 0 && toyratClimbUp_xp == 0
    aboveAppThreshold = Tracking.directionalVelocity(:,1) > 2; %pixels per frame
elseif mouse_toyrat == 1 | toyratClimbUp_xp == 1
    aboveAppThreshold = Tracking.diffDistanceMouseToyRat(:,1) < -2; %pixels per frame    
elseif shockgrid_xp == 1
    aboveAppThreshold = Tracking.diffDistanceShockgrid(:,1) < -2 & Tracking.mouse_position(:,1) < shockgrid_boundary(1); %pixels per frame
elseif moving_shockgrid == 1 
%separate rules for when shockgrid is moving and when it is stationary 
    aboveAppThreshold_gridMove = gridMoveIndices == 1 & abs(Tracking.mouseVel(:,1)) > 2 & (abs(Tracking.mouseNose(:,1)-Tracking.gridCenter) < abs(Tracking.mouseTailbase(:,1)-Tracking.gridCenter)); %was 5 for HD; if grid is moving, approach defined as mouse moving, and facing towards the grid.
    aboveAppThreshold_gridNoMoveRight = gridMoveIndices == 0 & Tracking.diffDistMouseShockgrid > 2 & Tracking.gridCenter > 300; %was 1 for SD, 5 for HD, otherwise treat is as you would escape from stationary rat or shockgrid with low threshold
    aboveAppThreshold_gridNoMoveLeft = gridMoveIndices == 0 & Tracking.diffDistMouseShockgrid < -2 & Tracking.gridCenter < 300; %was 1 for SD, 5 for HD, otherwise treat is as you would escape from stationary rat or shockgrid with low threshold
    aboveAppThreshold = aboveAppThreshold_gridMove + aboveAppThreshold_gridNoMoveLeft + aboveAppThreshold_gridNoMoveRight;
if isHD==1
    aboveAppThreshold_gridMove = gridMoveIndices == 1 & abs(Tracking.mouseVel(:,1)) > 5 & (abs(Tracking.mouseNose(:,1)-Tracking.gridCenter) < abs(Tracking.mouseTailbase(:,1)-Tracking.gridCenter)); %was 5 for HD; if grid is moving, approach defined as mouse moving, and facing towards the grid.
    aboveAppThreshold_gridNoMoveRight = gridMoveIndices == 0 & Tracking.diffDistMouseShockgrid > 5 & Tracking.gridCenter > 500; %was 1 for SD, 5 for HD, otherwise treat is as you would escape from stationary rat or shockgrid with low threshold
    aboveAppThreshold_gridNoMoveLeft = gridMoveIndices == 0 & Tracking.diffDistMouseShockgrid < -5 & Tracking.gridCenter < 500; %was 1 for SD, 5 for HD, otherwise treat is as you would escape from stationary rat or shockgrid with low threshold
    aboveAppThreshold = aboveAppThreshold_gridMove + aboveAppThreshold_gridNoMoveLeft + aboveAppThreshold_gridNoMoveRight;    
end
end
end

if flipRatSide == 1
if mouse_toyrat == 1 
    aboveAppThreshold = Tracking.diffDistanceMouseToyRat(:,1) < -2; %pixels per frame
elseif mouse_rat == 1
    aboveAppThreshold = Tracking.directionalVelocity(:,1) < -2; %pixels per frame
elseif shockgrid_xp == 1
    aboveAppThreshold = Tracking.diffDistanceShockgrid(:) < -2 & Tracking.mouse_position(:,1) > shockgrid_boundary(1); %pixels per frame
end
end
 
% Make a moving window to count the number of frames minimum velocity.
windowWidth = round(sampleRate .* 1);
counts = conv(aboveAppThreshold, ones(1, windowWidth), 'same');
approachIndices = counts >= round(sampleRate * .5); %changed this to .75 from .25

approachSwitch = diff(approachIndices);
approachSwitchIndices = find(approachSwitch);

if approachIndices(length(approachIndices)) == 1
    approachIndices(approachSwitchIndices(end):length(approachIndices)) = 0;
elseif approachIndices(1) == 1
    approachIndices(1:approachSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(approachSwitchIndices)
    if approachSwitch(approachSwitchIndices(i)) == 1
        approachStart(jj) = approachSwitchIndices(i);
        jj = jj+1;
    elseif approachSwitch(approachSwitchIndices(i)) == -1
        approachEnd(kk) = approachSwitchIndices(i);
        kk = kk+1;
    end
end

%remove indices if start & end points not complete bookends (acc by start or end)
if exist('approachStart')
if approachStart(1) > approachEnd(1)
    approachIndices(1:approachStart(1)) = 0;
end
if approachStart(end) > approachEnd(end)
    approachIndices(approachStart(end):length(approachIndices)) = 0;
end
end

 clearvars approachSwitch approachSwitchIndices approachStart approachEnd;
 approachSwitch = diff(approachIndices);
 approachSwitchIndices = find(approachSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(approachSwitchIndices)
     if approachSwitch(approachSwitchIndices(i)) == 1
         approachStart(jj) = approachSwitchIndices(i);
         jj = jj+1;
     elseif approachSwitch(approachSwitchIndices(i)) == -1
         approachEnd(kk) = approachSwitchIndices(i);
         kk = kk+1;
     end
 end

 if exist('approachStart')
%adjust minimum duration
 for i = 1:length(approachStart)
     if approachEnd(i) - approachStart(i) < round(sampleRate .* .5) %set the minimum duration here -- fraction of a second
         approachIndices(approachStart(i):approachEnd(i)) = 0;
     end
 end
 
% %update start and end indices
 clearvars approachSwitch approachSwitchIndices approachStart approachEnd;
 approachSwitch = diff(approachIndices);
 approachSwitchIndices = find(approachSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(approachSwitchIndices)
     if approachSwitch(approachSwitchIndices(i)) == 1
         approachStart(jj) = approachSwitchIndices(i);
         jj = jj+1;
     elseif approachSwitch(approachSwitchIndices(i)) == -1
         approachEnd(kk) = approachSwitchIndices(i);
         kk = kk+1;
     end
 end
 end

%fix if does not meet minimum distance criterion
if mouse_rat == 1
    
if exist('approachEnd')        
        if length(approachEnd) > length(approachStart)
            approachEnd = approachEnd(2:end);
        end
for i = 1:length(approachEnd)
    if Tracking.distanceMouseRat(approachStart(i)) - Tracking.distanceMouseRat(approachEnd(i)) < minPixCloserApp
        approachIndices(approachStart(i):approachEnd(i)) = 0;
    end
end
end
end

if mouse_toyrat == 1
    
if exist('approachEnd')        
        if length(approachEnd) > length(approachStart)
            approachEnd = approachEnd(2:end);
        end
for i = 1:length(approachEnd)
    if Tracking.distanceMouseToyRat(approachStart(i)) - Tracking.distanceMouseToyRat(approachEnd(i)) < minPixCloserApp
        approachIndices(approachStart(i):approachEnd(i)) = 0;
    end
end
end
end
end

%% approach burrow
if burrow_xp == 1 | burrow_norat_xp == 1
minAppStartDistance = 100; 
minAppEndDistance = 60; 
minPixCloserApp = 40;  
minVelApp = -3;
% Find indices where distance is increasing between mouse and rat above a threshold:
velocityBurrow = diff(Tracking.distanceFromBurrow);
aboveAppBurrThreshold = velocityBurrow < minVelApp; %pixels per frame
% Make a moving window to count the number of frames minimum velocity.
windowWidth = round(sampleRate .* 1.5);
counts = conv(aboveAppBurrThreshold, ones(1, windowWidth), 'same');
approachBurrIndices = counts >= round(sampleRate * .25);

approachBurrSwitch = diff(approachBurrIndices);
approachBurrSwitchIndices = find(approachBurrSwitch);

if approachBurrIndices(length(approachBurrIndices)) == 1
    approachBurrIndices(approachBurrSwitchIndices(end):length(approachBurrIndices)) = 0;
elseif approachBurrIndices(1) == 1
    approachBurrIndices(1:approachBurrSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(approachBurrSwitchIndices)
    if approachBurrSwitch(approachBurrSwitchIndices(i)) == 1
        approachBurrStart(jj) = approachBurrSwitchIndices(i);
        jj = jj+1;
    elseif approachBurrSwitch(approachBurrSwitchIndices(i)) == -1
        approachBurrEnd(kk) = approachBurrSwitchIndices(i);
        kk = kk+1;
    end
end

%remove indices if start & end points not complete bookends (acc by start or end)
if exist('approachBurrStart')
if approachBurrStart(1) > approachBurrEnd(1)
    approachBurrIndices(1:approachBurrStart(1)) = 0;
end
if approachBurrStart(end) > approachBurrEnd(end)
    approachBurrIndices(approachBurrStart(end):length(approachBurrIndices)) = 0;
end

% %update start and end indices
 clearvars approachBurrSwitch approachBurrSwitchIndices approachBurrStart approachBurrEnd;
 approachBurrSwitch = diff(approachBurrIndices);
 approachBurrSwitchIndices = find(approachBurrSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(approachBurrSwitchIndices)
     if approachBurrSwitch(approachBurrSwitchIndices(i)) == 1
         approachBurrStart(jj) = approachBurrSwitchIndices(i);
         jj = jj+1;
     elseif approachBurrSwitch(approachBurrSwitchIndices(i)) == -1
         approachBurrEnd(kk) = approachBurrSwitchIndices(i);
         kk = kk+1;
     end
 end

% %adjust minimum duration
 for i = 1:length(approachBurrStart)
     if approachBurrEnd(i) - approachBurrStart(i) < round(sampleRate .* .5) %set the minimum duration here -- fraction of a second
         approachBurrIndices(approachBurrStart(i):approachBurrEnd(i)) = 0;
     end
 end
 
% remove if preceded by escape by 3 seconds
if exist('escapeStart')
 for i = 1:length(approachBurrStart)
     if approachBurrStart(i) > sampleRate .* 3;
        if sum(escapeIndices(approachBurrStart(i) - sampleRate .* 3:approachBurrStart(i))) > 0;
          approachBurrIndices(approachBurrStart(i):approachBurrEnd(i)) = 0;
        end
     end
 end
end

% %update start and end indices
 clearvars approachBurrSwitch approachBurrSwitchIndices approachBurrStart approachBurrEnd;
 approachBurrSwitch = diff(approachBurrIndices);
 approachBurrSwitchIndices = find(approachBurrSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(approachBurrSwitchIndices)
     if approachBurrSwitch(approachBurrSwitchIndices(i)) == 1
         approachBurrStart(jj) = approachBurrSwitchIndices(i);
         jj = jj+1;
     elseif approachBurrSwitch(approachBurrSwitchIndices(i)) == -1
         approachBurrEnd(kk) = approachBurrSwitchIndices(i);
         kk = kk+1;
     end
 end
 
 % %remove approaches that don't start or end within a certain distance of the
 % burrow
 for i = 1:length(approachBurrStart)
     if Tracking.distanceFromBurrow(approachBurrEnd(i)) > minAppEndDistance %|| Tracking.distanceFromBurrow(approachBurrStart(i)) < minAppStartDistance  %set the minimum duration here -- fraction of a second
         approachBurrIndices(approachBurrStart(i):approachBurrEnd(i)) = 0;
     end
 end
% 
% %update start and end indices
 clearvars approachBurrSwitch approachBurrSwitchIndices approachBurrStart approachBurrEnd;
 approachBurrSwitch = diff(approachBurrIndices);
 approachBurrSwitchIndices = find(approachBurrSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(approachBurrSwitchIndices)
     if approachBurrSwitch(approachBurrSwitchIndices(i)) == 1
         approachBurrStart(jj) = approachBurrSwitchIndices(i);
         jj = jj+1;
     elseif approachBurrSwitch(approachBurrSwitchIndices(i)) == -1
         approachBurrEnd(kk) = approachBurrSwitchIndices(i);
         kk = kk+1;
     end
 end
 
%fix if does not meet minimum distance criterion
for i = 1:length(approachBurrStart)
    if Tracking.distanceFromBurrow(approachBurrStart(i)) - Tracking.distanceFromBurrow(approachBurrEnd(i)) < minPixCloserApp
        approachBurrIndices(approachBurrStart(i):approachBurrEnd(i)) = 0;
    end
end

% %update start and end indices
 clearvars approachBurrSwitch approachBurrSwitchIndices approachBurrStart approachBurrEnd;
 approachBurrSwitch = diff(approachBurrIndices);
 approachBurrSwitchIndices = find(approachBurrSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(approachBurrSwitchIndices)
     if approachBurrSwitch(approachBurrSwitchIndices(i)) == 1
         approachBurrStart(jj) = approachBurrSwitchIndices(i);
         jj = jj+1;
     elseif approachBurrSwitch(approachBurrSwitchIndices(i)) == -1
         approachBurrEnd(kk) = approachBurrSwitchIndices(i);
         kk = kk+1;
     end
 end
end
end

%% depart burrow
if burrow_xp == 1 | burrow_norat_xp == 1
maxDepStartDistance = 60; %set minimum depart start distance
minDepEndDistance = 80; %set the minimum depart end distance
minPixCloserApp = 40;  %set the minimum increase in distance from rat
%appDistThresh = 175; %set minimum approach end distance to count as approach towards rat (must pass a certain threshold)
minVelDep = 2;
% Find indices where distance is increasing between mouse and rat above a threshold:
velocityBurrow = diff(Tracking.distanceFromBurrow);
aboveDepBurrThreshold = velocityBurrow > minVelDep; %pixels per frame
% Make a moving window to count the number of frames minimum velocity.
windowWidth = round(sampleRate .* 1.5);
counts = conv(aboveDepBurrThreshold, ones(1, windowWidth), 'same');
departBurrIndices = counts >= round(sampleRate * .25);

departBurrSwitch = diff(departBurrIndices);
departBurrSwitchIndices = find(departBurrSwitch);

if departBurrIndices(length(departBurrIndices)) == 1
    departBurrIndices(departBurrSwitchIndices(end):length(departBurrIndices)) = 0;
elseif departBurrIndices(1) == 1
    departBurrIndices(1:departBurrSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(departBurrSwitchIndices)
    if departBurrSwitch(departBurrSwitchIndices(i)) == 1
        departBurrStart(jj) = departBurrSwitchIndices(i);
        jj = jj+1;
    elseif departBurrSwitch(departBurrSwitchIndices(i)) == -1
        departBurrEnd(kk) = departBurrSwitchIndices(i);
        kk = kk+1;
    end
end

 %remove indices if start & end points not complete bookends (acc by start or end)
if exist('departBurrStart')
if departBurrStart(1) > departBurrEnd(1)
    departBurrIndices(1:departBurrStart(1)) = 0;
end
if departBurrStart(end) > departBurrEnd(end)
    departBurrIndices(departBurrStart(end):length(departBurrIndices)) = 0;
end

%update start and end indices
 clearvars departBurrSwitch departBurrSwitchIndices departBurrStart departBurrEnd;
 departBurrSwitch = diff(departBurrIndices);
 departBurrSwitchIndices = find(departBurrSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(departBurrSwitchIndices)
     if departBurrSwitch(departBurrSwitchIndices(i)) == 1
         departBurrStart(jj) = departBurrSwitchIndices(i);
         jj = jj+1;
     elseif departBurrSwitch(departBurrSwitchIndices(i)) == -1
         departBurrEnd(kk) = departBurrSwitchIndices(i);
         kk = kk+1;
     end
 end

% %adjust minimum duration
 for i = 1:length(departBurrStart)
     if departBurrEnd(i) - departBurrStart(i) < round(sampleRate .* .25) %set the minimum duration here -- fraction of a second
         departBurrIndices(departBurrStart(i):departBurrEnd(i)) = 0;
     end
 end

%remove departs that don't start or end within a certain distance of the burrow
 for i = 1:length(departBurrStart)
     if Tracking.distanceFromBurrow(departBurrEnd(i)) < minDepEndDistance | Tracking.distanceFromBurrow(departBurrStart(i)) > maxDepStartDistance  %set the minimum duration here -- fraction of a second
         departBurrIndices(departBurrStart(i):departBurrEnd(i)) = 0;
     end
 end
  
% %update start and end indices
 clearvars departBurrSwitch departBurrSwitchIndices departBurrStart departBurrEnd;
 departBurrSwitch = diff(departBurrIndices);
 departBurrSwitchIndices = find(departBurrSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(departBurrSwitchIndices)
     if departBurrSwitch(departBurrSwitchIndices(i)) == 1
         departBurrStart(jj) = departBurrSwitchIndices(i);
         jj = jj+1;
     elseif departBurrSwitch(departBurrSwitchIndices(i)) == -1
         departBurrEnd(kk) = departBurrSwitchIndices(i);
         kk = kk+1;
     end
 end
 
end
end

%% freeze

if artificial_prey==0

freezeMoveThresh = .2; 
freezeAngleThresh = .03; 

if EPM == 0 && CO2 == 0 && hotplate_xp == 0 && ratClimbUp_xp == 0 && moving_shockgrid == 0 && toyratClimbUp_xp==0
if flipRatSide == 0
withinFreezeThreshold = Tracking.mouseTailbaseVel(:) < freezeMoveThresh ...
    & Tracking.mouseNoseVel(:) < freezeMoveThresh ...
    & abs(Tracking.positiondiffMouseAngle(:)) < freezeAngleThresh %& ...
    %Tracking.distEarsTailbase(:) < (strThreshEars-0) %& Tracking.mouseNose(:,1) > Tracking.mouseTailbase(:,1); %and must be facing towards rat if on the right and not stretching
elseif flipRatSide == 1
withinFreezeThreshold = Tracking.mouseTailbaseVel(:) < freezeMoveThresh & Tracking.mouseNoseVel(:) < freezeMoveThresh & abs(Tracking.positiondiffMouseAngle(:)) < freezeAngleThresh %& ...
    %Tracking.distNoseTailbase(:) < (strThresh-0);% & Tracking.mouseNose(:,1) < Tracking.mouseTailbase(:,1); %and must be facing towards rat if on the left and not stretching
end
end

if moving_shockgrid == 1
  if miniscope == 0
    withinFreezeThreshold = Tracking.mouseTailbaseVel(:) < freezeMoveThresh & Tracking.mouseNoseVel(:) < freezeMoveThresh & abs(Tracking.positiondiffMouseAngle(:)) < freezeAngleThresh & Tracking.distNoseTailbase(:) < ((strThresh)-5); %and not stretching    
  elseif miniscope == 1
    withinFreezeThreshold = Tracking.mouseTailbaseVel(:) < freezeMoveThresh & Tracking.mouseNoseVel(:) < freezeMoveThresh & abs(Tracking.positiondiffMouseAngle(:)) < freezeAngleThresh & Tracking.distEarsTailbase(:) < (strThreshEars); %and not stretching 
  end
end

if ratClimbUp_xp == 1 | toyratClimbUp_xp == 1
freezeMoveThresh = .4; %was .3 -- set movement maximum for freeze
freezeAngleThresh = .05; % was .05 , set angular movement max for freeze (of head)

    withinFreezeThreshold = Tracking.mouseTailbaseVel(:) < freezeMoveThresh & Tracking.mouseNoseVel(:) < freezeMoveThresh & ...
        Tracking.mouseNose(:,1) < Tracking.mouseTailbase(:,1) %& Tracking.distNoseTailbase(:) < (strThresh); %must be facing the rat and not stretching
end

if hotplate_xp_CNN == 1
freezeMoveThresh = .8;
%withinFreezeThreshold = Tracking.mouseVel < freezeMoveThresh & Tracking.mouse_position(1:end,2) > Bounds.lowerTrack(1,2);

withinFreezeThreshold = Tracking.mouseNoseVel < freezeMoveThresh & Tracking.mouse_position(1:end,2) > Bounds.lowerTrack(1,2);

end

if CO2 == 1 | hotplate_xp == 1
freezeMoveThresh = .8;
%withinFreezeThreshold = Tracking.mouseVel < freezeMoveThresh & Tracking.mouse_position(1:end,2) > Bounds.lowerTrack(1,2);

withinFreezeThreshold = Tracking.mouseVel < freezeMoveThresh & Tracking.mouse_position(1:end,2) > Bounds.lowerTrack(1,2);

end


if EPM == 1
withinFreezeThreshold = Tracking.mouseTailbaseVel(:) < freezeMoveThresh & Tracking.positiondiffMouseAngle(:) < freezeAngleThresh & ...
   Tracking.distNoseTailbase(:) < (strThresh-2); % and must not be stretching, or even borderline stretching
end

windowWidth = round(sampleRate .* 1);
counts = conv(withinFreezeThreshold, ones(1, windowWidth), 'same');
freezeIndices = counts >= round(sampleRate .* .5);

%find freeze start and end indices
freezeSwitch = diff(freezeIndices);
freezeSwitchIndices = find(freezeSwitch);

if freezeIndices(length(freezeIndices)) == 1
    freezeIndices(freezeSwitchIndices(end):length(freezeIndices)) = 0;
elseif freezeIndices(1) == 1
    freezeIndices(1:freezeSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(freezeSwitchIndices)
    if freezeSwitch(freezeSwitchIndices(i)) == 1
        freezeStart(jj) = freezeSwitchIndices(i);
        jj = jj+1;
    elseif freezeSwitch(freezeSwitchIndices(i)) == -1
        freezeEnd(kk) = freezeSwitchIndices(i);
        kk = kk+1;
    end
end
     if freezeIndices(length(freezeIndices)) == 1
         freezeEnd(kk) = length(freezeIndices);
     end
     
%remove indices if start & end points not complete bookends (acc by start or end)
if exist('freezeStart')
if freezeStart(1) > freezeEnd(1)
    freezeIndices(1:freezeStart(1)) = 0;
end
if freezeStart(end) > freezeEnd(end)
    freezeIndices(freezeStart(end):length(freezeIndices)) = 0;
end

% %update approach start and end indices
 clearvars freezeSwitch freezeSwitchIndices freezeStart freezeEnd;
% 
 freezeSwitch = diff(freezeIndices);
 freezeSwitchIndices = find(freezeSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(freezeSwitchIndices)
     if freezeSwitch(freezeSwitchIndices(i)) == 1
         freezeStart(jj) = freezeSwitchIndices(i);
         jj = jj+1;
     elseif freezeSwitch(freezeSwitchIndices(i)) == -1
         freezeEnd(kk) = freezeSwitchIndices(i);
         kk = kk+1;
     end
 end

%adjust minimum freeze duration
    if length(freezeStart) > length(freezeEnd)
        freezeStart = freezeStart(1:end-1);
    end
    
 for i = 1:length(freezeStart)
     if freezeEnd(i) - freezeStart(i) < round(sampleRate .* .5) %set the minimum duration here
         freezeIndices(freezeStart(i):freezeEnd(i)) = 0;
     end
 end
 
 %specify for climb up that freeze can't happen on the platform
if ratClimbUp_xp == 1
    for i = 1:length(freezeStart)
       if Tracking.mouse_position(freezeStart(i),2) < 350
           freezeIndices(freezeStart(i):freezeEnd(i)) = 0;
       end
    end
     
end

%
% %update approach start and end indices
 clearvars freezeSwitch freezeSwitchIndices freezeStart freezeEnd;
% 
 freezeSwitch = diff(freezeIndices);
 freezeSwitchIndices = find(freezeSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(freezeSwitchIndices)
     if freezeSwitch(freezeSwitchIndices(i)) == 1
         freezeStart(jj) = freezeSwitchIndices(i);
         jj = jj+1;
     elseif freezeSwitch(freezeSwitchIndices(i)) == -1
         freezeEnd(kk) = freezeSwitchIndices(i);
         kk = kk+1;
     end
 end
end

Tracking.freezeFrac = sum(freezeIndices) / length(Tracking.mouse_position);

end

%% REAR UP WALL

if EPM == 0 && ratClimbUp_xp == 0 && toyratClimbUp_xp == 0 && hotplate_xp == 0 && CO2 == 0 && RTPP == 0 && hotplate_xp_CNN == 0 && fear_cond_chamber == 0  && artificial_prey==0 && open_field==0
% Find indices where mouse nose passes beyond constraint rectangle:
rearingIndices = Tracking.mouseNose(:,1) > rearing_constraint(2,1) | Tracking.mouseNose(:,1) < rearing_constraint(1,1) | ...
    Tracking.mouseNose(:,2) < rearing_constraint(1,2) | Tracking.mouseNose(:,2) > rearing_constraint(2,2);
% Make a moving window to count the number of frames minimum velocity.
windowWidth = round(sampleRate .* 1.0);
counts = conv(rearingIndices, ones(1, windowWidth), 'same');
rearingIndices = counts >= round(sampleRate * .35);

rearingSwitch = diff(rearingIndices);
rearingSwitchIndices = find(rearingSwitch);

if rearingIndices(length(rearingIndices)) == 1
    rearingIndices(rearingSwitchIndices(end):length(rearingIndices)) = 0;
elseif rearingIndices(1) == 1
    rearingIndices(1:rearingSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(rearingSwitchIndices)
    if rearingSwitch(rearingSwitchIndices(i)) == 1
        rearingStart(jj) = rearingSwitchIndices(i);
        jj = jj+1;
    elseif rearingSwitch(rearingSwitchIndices(i)) == -1
        rearingEnd(kk) = rearingSwitchIndices(i);
        kk = kk+1;
    end
end

%remove single indices of behavior
if sum(rearingIndices)==0
    rearingIndices(find(rearingIndices==1)) = 0;
    clearvars rearingStart rearingEnd
end

%remove indices if start & end points not complete bookends (acc by start or end)
if exist('rearingStart')
if rearingStart(1) > rearingEnd(1)
    rearingIndices(1:rearingStart(1)) = 0;
end
if rearingStart(end) > rearingEnd(end)
    rearingIndices(rearingStart(end):length(rearingIndices)) = 0;
end

% %update approach start and end indices
 clearvars rearingSwitch rearingSwitchIndices rearingStart rearingEnd;
% 
 rearingSwitch = diff(rearingIndices);
 rearingSwitchIndices = find(rearingSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(rearingSwitchIndices)
     if rearingSwitch(rearingSwitchIndices(i)) == 1
         rearingStart(jj) = rearingSwitchIndices(i);
         jj = jj+1;
     elseif rearingSwitch(rearingSwitchIndices(i)) == -1
         rearingEnd(kk) = rearingSwitchIndices(i);
         kk = kk+1;
     end
 end
 
     if rearingIndices(length(rearingIndices)) == 1
         rearingEnd(kk) = length(rearingIndices);
     end

 %adjust minimum duration
 if exist('rearingEnd', 'var')
 for i = 1:length(rearingStart)
     if rearingEnd(i) - rearingStart(i) < round(sampleRate .* .25) %set the minimum duration here
         rearingIndices(rearingStart(i):rearingEnd(i)) = 0;
     end
 end
 end
 
% %update approach start and end indices
 clearvars rearingSwitch rearingSwitchIndices rearingStart rearingEnd;
 rearingSwitch = diff(rearingIndices);
 rearingSwitchIndices = find(rearingSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(rearingSwitchIndices)
     if rearingSwitch(rearingSwitchIndices(i)) == 1
         rearingStart(jj) = rearingSwitchIndices(i);
         jj = jj+1;
     elseif rearingSwitch(rearingSwitchIndices(i)) == -1
         rearingEnd(kk) = rearingSwitchIndices(i);
         kk = kk+1;
     end
 end
 
end
end

%% HIDE INSIDE BURROW AND NON-BURROW
%burrow
if burrow_xp == 1 | burrow_norat_xp == 1
    %count indices that are within bounding box of burrow-- mouse nose and
    %tailbase
    insideBurrowIndices = Tracking.mouseBetwEars(:,1) > burrow_boundary(1,1) & Tracking.mouseBetwEars(:,1) < burrow_boundary(2,1) & ...
        Tracking.mouseBetwEars(:,2) > burrow_boundary(1,2) & Tracking.mouseBetwEars(:,2) < burrow_boundary(2,2)
        %& Tracking.mouseTailbase(:,1) > burrow_boundary(1,1) & Tracking.mouseTailbase(:,1) < burrow_boundary(2,1) & ...
        %Tracking.mouseTailbase(:,2) > burrow_boundary(1,2) & Tracking.mouseTailbase(:,2) < burrow_boundary(2,2);

    % Make a moving window to count the number of frames.
windowWidth = round(sampleRate .* 1.5);
counts = conv(insideBurrowIndices, ones(1, windowWidth), 'same');
insideBurrowIndices = counts >= round(sampleRate * .5);

insideBurrowSwitch = diff(insideBurrowIndices);
insideBurrowSwitchIndices = find(insideBurrowSwitch);

if insideBurrowIndices(length(insideBurrowIndices)) == 1
    insideBurrowIndices(insideBurrowSwitchIndices(end):length(insideBurrowIndices)) = 0;
elseif insideBurrowIndices(1) == 1
    insideBurrowIndices(1:insideBurrowSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(insideBurrowSwitchIndices)
    if insideBurrowSwitch(insideBurrowSwitchIndices(i)) == 1
        insideBurrowStart(jj) = insideBurrowSwitchIndices(i);
        jj = jj+1;
    elseif insideBurrowSwitch(insideBurrowSwitchIndices(i)) == -1
        insideBurrowEnd(kk) = insideBurrowSwitchIndices(i);
        kk = kk+1;
    end
end

%remove indices if start & end points not complete bookends (acc by start or end)
if exist('insideBurrowStart')
if insideBurrowStart(1) > insideBurrowEnd(1)
    insideBurrowIndices(1:insideBurrowStart(1)) = 0;
end
if insideBurrowStart(end) > insideBurrowEnd(end)
    insideBurrowIndices(insideBurrowStart(end):length(insideBurrowIndices)) = 0;
end

% %update approach start and end indices
 clearvars insideBurrowSwitch insideBurrowSwitchIndices insideBurrowStart insideBurrowEnd;
% 
 insideBurrowSwitch = diff(insideBurrowIndices);
 insideBurrowSwitchIndices = find(insideBurrowSwitch);
% 
 jj = 1; kk = 1;
 for i=1:length(insideBurrowSwitchIndices)
     if insideBurrowSwitch(insideBurrowSwitchIndices(i)) == 1
         insideBurrowStart(jj) = insideBurrowSwitchIndices(i);
         jj = jj+1;
     elseif insideBurrowSwitch(insideBurrowSwitchIndices(i)) == -1
         insideBurrowEnd(kk) = insideBurrowSwitchIndices(i);
         kk = kk+1;
     end
 end

if 0 == 1
%non-burrow

    %count indices that are within bounding box of burrow-- mouse nose and
    %tailbase
    insideNonBurrowIndices = Tracking.mouseBetwEars(:,1) > nonburrow_boundary(1,1) & Tracking.mouseBetwEars(:,1) < nonburrow_boundary(2,1) & ...
        Tracking.mouseBetwEars(:,2) > nonburrow_boundary(1,2) & Tracking.mouseBetwEars(:,2) < nonburrow_boundary(2,2)
        %& Tracking.mouseTailbase(:,1) > burrow_boundary(1,1) & Tracking.mouseTailbase(:,1) < burrow_boundary(2,1) & ...
        %Tracking.mouseTailbase(:,2) > burrow_boundary(1,2) & Tracking.mouseTailbase(:,2) < burrow_boundary(2,2);

    % Make a moving window to count the number of frames.
windowWidth = round(sampleRate .* 1.5);
counts = conv(insideNonBurrowIndices, ones(1, windowWidth), 'same');
insideNonBurrowIndices = counts >= round(sampleRate * .5);

insideNonBurrowSwitch = diff(insideNonBurrowIndices);
insideNonBurrowSwitchIndices = find(insideNonBurrowSwitch);

if insideNonBurrowIndices(length(insideNonBurrowIndices)) == 1
    insideNonBurrowIndices(insideNonBurrowSwitchIndices(end):length(insideNonBurrowIndices)) = 0;
elseif insideNonBurrowIndices(1) == 1
    insideNonBurrowIndices(1:insideNonBurrowSwitchIndices(1)) = 0;
end

jj = 1; kk = 1;
for i=1:length(insideNonBurrowSwitchIndices)
    if insideNonBurrowSwitch(insideNonBurrowSwitchIndices(i)) == 1
        insideNonBurrowStart(jj) = insideNonBurrowSwitchIndices(i);
        jj = jj+1;
    elseif insideNonBurrowSwitch(insideNonBurrowSwitchIndices(i)) == -1
        insideNonBurrowEnd(kk) = insideNonBurrowSwitchIndices(i);
        kk = kk+1;
    end
end
end
end
end

%% Metrics for Elevated Plus Maze
if EPM == 1
Tracking.EPMRegion = [];
load('PlusArmN.mat'); load('PlusArmS.mat'); load('PlusArmE.mat'); load('PlusArmW.mat');
% for each frame, determine the arm of occupancy
for i = 1:length(Tracking.mouse_position)
    %North
   if Tracking.mouse_position(i,1) > PlusArmN(1,1) && Tracking.mouse_position(i,1) < PlusArmN(2,1) && Tracking.mouse_position(i,2) > PlusArmN(1,2) && Tracking.mouse_position(i,2) < PlusArmN(2,2) && ...
      Tracking.mouseTailbase(i,1) > PlusArmN(1,1) && Tracking.mouseTailbase(i,1) < PlusArmN(2,1) && Tracking.mouseTailbase(i,2) > PlusArmN(1,2) && Tracking.mouseTailbase(i,2) < PlusArmN(2,2);     
        Tracking.EPMRegion(i) = 1;
   elseif Tracking.mouse_position(i,1) > PlusArmE(1,1) && Tracking.mouse_position(i,1) < PlusArmE(2,1) && Tracking.mouse_position(i,2) > PlusArmE(1,2) && Tracking.mouse_position(i,2) < PlusArmE(2,2) && ...
      Tracking.mouseTailbase(i,1) > PlusArmE(1,1) && Tracking.mouseTailbase(i,1) < PlusArmE(2,1) && Tracking.mouseTailbase(i,2) > PlusArmE(1,2) && Tracking.mouseTailbase(i,2) < PlusArmE(2,2);     
        Tracking.EPMRegion(i) = 2;
   elseif Tracking.mouse_position(i,1) > PlusArmS(1,1) && Tracking.mouse_position(i,1) < PlusArmS(2,1) && Tracking.mouse_position(i,2) > PlusArmS(1,2) && Tracking.mouse_position(i,2) < PlusArmS(2,2) && ...
      Tracking.mouseTailbase(i,1) > PlusArmS(1,1) && Tracking.mouseTailbase(i,1) < PlusArmS(2,1) && Tracking.mouseTailbase(i,2) > PlusArmS(1,2) && Tracking.mouseTailbase(i,2) < PlusArmS(2,2);     
        Tracking.EPMRegion(i) = 3;
   elseif Tracking.mouse_position(i,1) > PlusArmW(1,1) && Tracking.mouse_position(i,1) < PlusArmW(2,1) && Tracking.mouse_position(i,2) > PlusArmW(1,2) && Tracking.mouse_position(i,2) < PlusArmW(2,2) && ...
      Tracking.mouseTailbase(i,1) > PlusArmW(1,1) && Tracking.mouseTailbase(i,1) < PlusArmW(2,1) && Tracking.mouseTailbase(i,2) > PlusArmW(1,2) && Tracking.mouseTailbase(i,2) < PlusArmW(2,2);     
        Tracking.EPMRegion(i) = 4; 
   elseif Tracking.mouse_position(i,1) > PlusArmW(2,1) && Tracking.mouse_position(i,1) < PlusArmE(1,1) && Tracking.mouse_position(i,2) > PlusArmN(2,2) && Tracking.mouse_position(i,2) < PlusArmS(2,1) && ...
      Tracking.mouseTailbase(i,1) > PlusArmW(2,1) && Tracking.mouseTailbase(i,1) < PlusArmE(1,1) && Tracking.mouseTailbase(i,2) > PlusArmN(2,2) && Tracking.mouseTailbase(i,2) < PlusArmS(2,1);     
        Tracking.EPMRegion(i) = 5;
   else Tracking.EPMRegion(i) = 10;
   end
end

openArmNorthSouth = 0;
openArmEastWest = 1;

%find percent time in open and closed arms
%open arm
if openArmNorthSouth == 1
openCount = length(find(Tracking.EPMRegion == 1 | Tracking.EPMRegion == 3));
%closed arm
closedCount = length(find(Tracking.EPMRegion == 2 | Tracking.EPMRegion == 4));
%center
centerCount = length(find(Tracking.EPMRegion == 5));

%calculate closed and open arm indices
closedArmIndices = zeros(length(Tracking.mouse_position),1);
closedArmIndices(find(Tracking.EPMRegion == 2 | Tracking.EPMRegion == 4)) = 1;
openArmIndices = zeros(length(Tracking.mouse_position),1);
openArmIndices(find(Tracking.EPMRegion == 1 | Tracking.EPMRegion == 3)) = 1;
end

if openArmEastWest == 1
openCount = length(find(Tracking.EPMRegion == 2 | Tracking.EPMRegion == 4));
%closed arm
closedCount = length(find(Tracking.EPMRegion == 1 | Tracking.EPMRegion == 3));
%center
centerCount = length(find(Tracking.EPMRegion == 5));

%calculate closed and open arm indices
openArmIndices = zeros(length(Tracking.mouse_position),1);
openArmIndices(find(Tracking.EPMRegion == 2 | Tracking.EPMRegion == 4)) = 1;
closedArmIndices = zeros(length(Tracking.mouse_position),1);
closedArmIndices(find(Tracking.EPMRegion == 1 | Tracking.EPMRegion == 3)) = 1;
end

%calc percent occupancy
Tracking.PercentTimeClosedEPM = closedCount / length(Tracking.mouse_position);
Tracking.PercentTimeOpenEPM = openCount / length(Tracking.mouse_position);
Tracking.PercentTimeCenterEPM = centerCount / length(Tracking.mouse_position);

%find # entries into open arms
temp = diff(openArmIndices);
openArmEntryStart = find(temp == 1);

    %to avoid tracking jitter at the boundary or dropout false positives,
    %remove any entry preceded by 3 seconds by another entry
    idxToDelete = zeros(1,length(openArmEntryStart));
    for i = 1:length(openArmEntryStart)-1
        if openArmEntryStart(i+1) - openArmEntryStart(i) < 70;
            idxToDelete(i+1) = 1;
        end
    end
    openArmEntryStart(find(idxToDelete==1)) = [];
    
openArmEntryEnd = find(temp == -1);

if openArmEntryEnd(1) < openArmEntryStart(1)
    openArmEntryEnd(1) = [];
end

if length(idxToDelete) > length(openArmEntryEnd)
    idxToDelete = idxToDelete(1:end-1);
end
    openArmEntryEnd(find(idxToDelete==1)) = [];
    
if length(openArmEntryStart) > length(openArmEntryEnd)
    openArmEntryStart(end) = [];
end
    
clearvars idxToDelete  

temp = diff(closedArmIndices);
closedArmEntryStart = find(temp == 1);

    %to avoid tracking jitter at the boundary or dropout false positives,
    %remove any entry preceded by 3 seconds by another entry
    idxToDelete = zeros(1,length(closedArmEntryStart));
    for i = 1:length(closedArmEntryStart)-1
        if closedArmEntryStart(i+1) - closedArmEntryStart(i) < 70;
            idxToDelete(i+1) = 1;
        end
    end
    closedArmEntryStart(find(idxToDelete==1)) = [];

closedArmEntryEnd = find(temp == -1);

    
if closedArmEntryEnd(1) < closedArmEntryStart(1)
    closedArmEntryEnd(1) = [];
end

if length(idxToDelete) > length(openArmEntryEnd)
    idxToDelete = idxToDelete(1:end-1);
end

    closedArmEntryEnd(find(idxToDelete==1)) = [];

clearvars idxToDelete
   
    minSamples = 15;
    temp = openArmEntryEnd-openArmEntryStart;
    temp = find(temp < minSamples);
    openArmEntryStart(temp) = []; openArmEntryEnd(temp) = [];    
end

%% Metrics for Open Field
if open_field == 1
Tracking.OF_Region = [];
load('of_center.mat');
% for each frame, determine the arm of occupancy
for i = 1:length(Tracking.mouse_position)
    %Center of open field
   if Tracking.mouse_position(i,1) > of_center(1,1) && Tracking.mouse_position(i,1) < of_center(2,1) && Tracking.mouse_position(i,2) > of_center(1,2) && Tracking.mouse_position(i,2) < of_center(2,2) && ...
           Tracking.mouseTailbase(i,1) > of_center(1,1) && Tracking.mouseTailbase(i,1) < of_center(2,1) && Tracking.mouseTailbase(i,2) > of_center(1,2) && Tracking.mouseTailbase(i,2) < of_center(2,2);
        Tracking.OF_Region(i) = 1;
   else
        Tracking.OF_Region(i) = 0;
   end
end
Tracking.OF_Region = Tracking.OF_Region';
Tracking.OF_Region_FracTime = sum(Tracking.OF_Region) ./ length(Tracking.OF_Region); 
end

%% RTPP 

if RTPP == 1
load('Bounds.mat')
boundsLeftThreshold = Tracking.mouseNose(:,2) < Bounds.left(2,2) & Tracking.mouseNose(:,2) > Bounds.left(1,2) & Tracking.mouseNose(:,1) > Bounds.left(1,1) & Tracking.mouseNose(:,1) < Bounds.left(2,1) & ...
    Tracking.mouseTailbase(:,2) < Bounds.left(2,2) & Tracking.mouseTailbase(:,2) > Bounds.left(1,2) & Tracking.mouseTailbase(:,1) > Bounds.left(1,1) & Tracking.mouseTailbase(:,1) < Bounds.left(2,1);
boundsRightThreshold = Tracking.mouseNose(:,2) < Bounds.right(2,2) & Tracking.mouseNose(:,2) > Bounds.right(1,2) & Tracking.mouseNose(:,1) > Bounds.right(1,1) & Tracking.mouseNose(:,1) < Bounds.right(2,1) & ...
    Tracking.mouseTailbase(:,2) < Bounds.right(2,2) & Tracking.mouseTailbase(:,2) > Bounds.right(1,2) & Tracking.mouseTailbase(:,1) > Bounds.right(1,1) & Tracking.mouseTailbase(:,1) < Bounds.right(2,1);

Tracking.leftBoxSecsOcc = sum(boundsLeftThreshold)/25; %divide by framerate
Tracking.rightBoxSecsOcc = sum(boundsRightThreshold)/25;
Tracking.leftBoxFracOcc = sum(boundsLeftThreshold) ./ length(Tracking.mouse_position);
Tracking.rightBoxFracOcc = sum(boundsRightThreshold) ./ length(Tracking.mouse_position);
end

%% find head dips in EPM
if EPM == 1
   
    Tracking.mouse_position = Tracking.mouseNose;
    
    headDip = Tracking.mouse_position(:,1) > PlusArmNE(1,1) & Tracking.mouse_position(:,2) < PlusArmNE(2,2) | ...
        Tracking.mouse_position(:,1) > PlusArmSE(1,1) & Tracking.mouse_position(:,2) > PlusArmSE(1,2) | ...
        Tracking.mouse_position(:,1) < PlusArmNW(2,1) & Tracking.mouse_position(:,2) < PlusArmNW(2,2) | ...
        Tracking.mouse_position(:,1) < PlusArmSW(2,1) & Tracking.mouse_position(:,2) > PlusArmSW(1,2) | ...
        Tracking.mouse_position(:,2) < PlusArmNE(1,2) | ...
        Tracking.mouse_position(:,2) > PlusArmSE(2,2) | ...
        Tracking.mouse_position(:,1) < PlusArmNW(1,1) | ...
        Tracking.mouse_position(:,1) > PlusArmNE(2,1);
        %Tracking.mouse_position(:,1) > PlusArmNE(1,1) & Tracking.mouse_position(:,2) < PlusArmNE(2,2) | ...
    
    
    % Make a moving window to count the number of frames minimum velocity.
    windowWidth = round(sampleRate .* .4); %was * 1.5
    counts = conv(headDip, ones(1, windowWidth), 'same');
    headDipIndices = counts >= round(sampleRate * .3); %was .7

    headDipSwitch = diff(headDipIndices);
    headDipSwitchIndices = find(headDipSwitch);

    if headDipIndices(length(headDipIndices)) == 1
        headDipIndices(headDipSwitchIndices(end):length(headDipIndices)) = 0;
    elseif headDipIndices(1) == 1
        headDipIndices(1:headDipSwitchIndices(1)) = 0;
    end

    jj = 1; kk = 1;
    for i=1:length(headDipSwitchIndices)
        if headDipSwitch(headDipSwitchIndices(i)) == 1
            headDipStart(jj) = headDipSwitchIndices(i);
            jj = jj+1;
        elseif headDipSwitch(headDipSwitchIndices(i)) == -1
            headDipEnd(kk) = headDipSwitchIndices(i);
            kk = kk+1;
        end
    end

    %remove indices if start & end points not complete bookends (acc by start or end)
    if exist('headDipStart')
    if headDipStart(1) > headDipEnd(1)
        headDipIndices(1:headDipStart(1)) = 0;
    end
    if headDipStart(end) > headDipEnd(end)
        headDipIndices(headDipStart(end):length(headDipIndices)) = 0;
    end

    % %update start and end indices
     clearvars headDipSwitch headDipSwitchIndices headDipStart headDipEnd;
    % 
     headDipSwitch = diff(headDipIndices);
     headDipSwitchIndices = find(headDipSwitch);
     
    if headDipIndices(length(headDipIndices)) == 1
        headDipIndices(headDipSwitchIndices(end):length(headDipIndices)) = 0;
    elseif headDipIndices(1) == 1
        headDipIndices(1:headDipSwitchIndices(1)) = 0;
    end

    jj = 1; kk = 1;
    for i=1:length(headDipSwitchIndices)
        if headDipSwitch(headDipSwitchIndices(i)) == 1
            headDipStart(jj) = headDipSwitchIndices(i);
            jj = jj+1;
        elseif headDipSwitch(headDipSwitchIndices(i)) == -1
            headDipEnd(kk) = headDipSwitchIndices(i);
            kk = kk+1;
        end
    end
    
    end
end
    
%% save output in current folder

    clearvars -except gridCrossNumFrames idxToRight_All idxToLeft_All escapeSucc escTrialDistMouseGrid escTrialFracEnclosure escSuccTrialIndices escSuccTrialStart escSuccTrialEnd escUnsuccTrialIndices escUnsuccTrialStart escUnsuccTrialEnd rightwardsGridMove escVelMax escVelMean escapeSuccFrac escapeSucc retreatStart retreatEnd retreatIndices escapeStart_gridNoMove escapeEnd_gridNoMove escapeIndices_gridNoMove artificial_prey folderzzz folderzz framesToDelete hotplate_xp_CNN open_field opto_xp toyratClimbUp_xp hotplate_xp CO2 fear_cond_chamber ratClimbUp_xp deleteUserDefinedZone setByParentScript miniscope other_vid fear_cond_chamber ...
    simpleRat_xp moving_shockgrid mouse_toyrat burrow_xp RTPP shockgrid_xp EPM session_name folders nn jj burrow_norat_xp boundsApproach boundsStart boundsEnd boundsIndices closedArmEntryStart closedArmEntryEnd closedArmExitStart closedArmIndices openArmIndices openArmEntryStart openArmEntryEnd openArmExitStart crouchsniffStart crouchsniffEnd crouchsniffIndices departBurrStart departBurrEnd departBurrIndices insideBurrowStart insideBurrowEnd insideBurrowIndices approachBurrStart approachBurrEnd approachBurrIndices vertStart vertEnd vertIndices approachIndices approachStart approachEnd freezeIndices freezeStart ...
    allTimeUpHoseIndices boxStart boxEnd boxIndices headDipStart headDipEnd headDipIndices gridMoveStart gridMoveEnd moveStart moveEnd moveIndices gridMoveIndices freezeEnd climbUStart climbUEnd climbUIndices jumpStart jumpEnd jumpIndices climbHoseStart climbHoseEnd climbHoseIndices stretchIndices stretchStart stretchEnd escapeIndices escapeStart escapeEnd rearingStart rearingEnd rearingIndices Tracking insideNonBurrowEnd insideNonBurrowStart insideNonBurrowIndices sessions h

    save('Tracking', 'Tracking'); 
    clearvars Tracking
    save('Behavior', '-regexp', '^(?!(artificial_prey|moving_shockgrid|open_field|folders|jjj|nn|ans|hotplate_xp_CNN|RTPP|setByParentScript|EPM|jj|i|deleteUserDefinedZone|burrow_xp|fear_cond_chamber|burrow_norat_xp|hotplate_xp|CO2|miniscope|mouse_toyrat|other_vid|ratClimbUp_xp|shockgrid_xp|simpleRat_xp|fear_cond_chamber|toyratClimbUp_xp)$).')


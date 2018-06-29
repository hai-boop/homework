%%%demo of matlab%%%


function multiObjectTracking()

% create system objects used for reading video, detecting moving objects,
% and displaying the results
obj = setupSystemObjects(); %初始化函数
tracks = initializeTracks(); % create an empty array of tracks  %初始化轨迹对象

nextId = 1; % ID of the next track

% detect moving objects, and track them across video frames
while ~isDone(obj.reader)
    frame = readFrame();  %读取一帧
    [centroids, bboxes, mask] = detectObjects(frame); %前景检测
    predictNewLocationsOfTracks();  %根据位置进行卡尔曼预测
    [assignments, unassignedTracks, unassignedDetections] = ...
        detectionToTrackAssignment(); %匈牙利匹配算法进行匹配
    
    updateAssignedTracks();%分配好的轨迹更新
    updateUnassignedTracks();%未分配的轨迹更新
    deleteLostTracks();%删除丢掉的轨迹
    createNewTracks();%创建新轨迹
    
    displayTrackingResults();%结果展示
end


%% Create System Objects
% Create System objects used for reading the video frames, detecting
% foreground objects, and displaying results.

    function obj = setupSystemObjects()
        % Initialize Video I/O
        % Create objects for reading a video from a file, drawing the tracked
        % objects in each frame, and playing the video.
        
        % create a video file reader
        obj.reader = vision.VideoFileReader('yourvideo.avi');         %读入视频
        
        % create two video players, one to display the video,
        % and one to display the foreground mask
        obj.videoPlayer = vision.VideoPlayer('Position', [20, 400, 700, 400]);   %创建两个窗口
        obj.maskPlayer = vision.VideoPlayer('Position', [740, 400, 700, 400]);
        
        % Create system objects for foreground detection and blob analysis
        
        obj.detector = vision.ForegroundDetector('NumGaussians', 3, ...   %GMM进行前景检测，高斯核数目为3，前40帧为背景帧，域值为0.7
            'NumTrainingFrames', 40, 'MinimumBackgroundRatio', 0.7);   
        
        % Connected groups of foreground pixels are likely to correspond to moving
        % objects.  The blob analysis system object is used to find such groups
        % (called 'blobs' or 'connected components'), and compute their
        % characteristics, such as area, centroid, and the bounding box.
        
        obj.blobAnalyser = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...  %输出质心和外接矩形
            'AreaOutputPort', true, 'CentroidOutputPort', true, ...
            'MinimumBlobArea', 400);
    end

%% Initialize Tracks
   

    function tracks = initializeTracks()
        % create an empty array of tracks
        tracks = struct(...
            'id', {}, ...  %轨迹ID
            'bbox', {}, ... %外接矩形
            'kalmanFilter', {}, ...%轨迹的卡尔曼滤波器
            'age', {}, ...%总数量
            'totalVisibleCount', {}, ...%可视数量
            'consecutiveInvisibleCount', {});%不可视数量
    end

%% Read a Video Frame
% Read the next video frame from the video file.
    function frame = readFrame()
        frame = obj.reader.step();%激活读图函数
    end

%% Detect Objects 

    function [centroids, bboxes, mask] = detectObjects(frame)
        
        % detect foreground
        mask = obj.detector.step(frame);
        
        % apply morphological operations to remove noise and fill in holes
        mask = imopen(mask, strel('rectangle', [3,3]));%开运算
        mask = imclose(mask, strel('rectangle', [15, 15])); %闭运算
        mask = imfill(mask, 'holes');%填洞
        
        % perform blob analysis to find connected components
        [~, centroids, bboxes] = obj.blobAnalyser.step(mask);
    end

%% Predict New Locations of Existing Tracks
% Use the Kalman filter to predict the centroid of each track in the
% current frame, and update its bounding box accordingly.

    function predictNewLocationsOfTracks()
        for i = 1:length(tracks)
            bbox = tracks(i).bbox;
            
            % predict the current location of the track
            predictedCentroid = predict(tracks(i).kalmanFilter);%根据以前的轨迹，预测当前位置
            
            % shift the bounding box so that its center is at 
            % the predicted location
            predictedCentroid = int32(predictedCentroid) - bbox(3:4) / 2;
            tracks(i).bbox = [predictedCentroid, bbox(3:4)];%真正的当前位置
        end
    end

%% Assign Detections to Tracks

    function [assignments, unassignedTracks, unassignedDetections] = ...
            detectionToTrackAssignment()
        
        nTracks = length(tracks);
        nDetections = size(centroids, 1);
        
        % compute the cost of assigning each detection to each track
        cost = zeros(nTracks, nDetections);
        for i = 1:nTracks
            cost(i, :) = distance(tracks(i).kalmanFilter, centroids);%损失矩阵计算
        end
        
        % solve the assignment problem
        costOfNonAssignment = 20;
        [assignments, unassignedTracks, unassignedDetections] = ...
            assignDetectionsToTracks(cost, costOfNonAssignment);%匈牙利算法匹配
    end

%% Update Assigned Tracks


    function updateAssignedTracks()
        numAssignedTracks = size(assignments, 1);
        for i = 1:numAssignedTracks
            trackIdx = assignments(i, 1);
            detectionIdx = assignments(i, 2);
            centroid = centroids(detectionIdx, :);
            bbox = bboxes(detectionIdx, :);
            
            % correct the estimate of the object's location
            % using the new detection
            correct(tracks(trackIdx).kalmanFilter, centroid);
            
            % replace predicted bounding box with detected
            % bounding box
            tracks(trackIdx).bbox = bbox;
            
            % update track's age
            tracks(trackIdx).age = tracks(trackIdx).age + 1;
            
            % update visibility
            tracks(trackIdx).totalVisibleCount = ...
                tracks(trackIdx).totalVisibleCount + 1;
            tracks(trackIdx).consecutiveInvisibleCount = 0;
        end
    end

%% Update Unassigned Tracks
% Mark each unassigned track as invisible, and increase its age by 1.

    function updateUnassignedTracks()
        for i = 1:length(unassignedTracks)
            ind = unassignedTracks(i);
            tracks(ind).age = tracks(ind).age + 1;
            tracks(ind).consecutiveInvisibleCount = ...
                tracks(ind).consecutiveInvisibleCount + 1;
        end
    end

%% Delete Lost Tracks

    function deleteLostTracks()
        if isempty(tracks)
            return;
        end
        
        invisibleForTooLong = 10;
        ageThreshold = 8;
        
        % compute the fraction of the track's age for which it was visible
        ages = [tracks(:).age];
        totalVisibleCounts = [tracks(:).totalVisibleCount];
        visibility = totalVisibleCounts ./ ages;
        
        % find the indices of 'lost' tracks
        lostInds = (ages < ageThreshold & visibility < 0.6) | ...
            [tracks(:).consecutiveInvisibleCount] >= invisibleForTooLong;
        
        % delete lost tracks
        tracks = tracks(~lostInds);
    end

%% Create New Tracks

    function createNewTracks()
        centroids = centroids(unassignedDetections, :);
        bboxes = bboxes(unassignedDetections, :);
        
        for i = 1:size(centroids, 1)
            
            centroid = centroids(i,:);
            bbox = bboxes(i, :);
            
            % create a Kalman filter object
            kalmanFilter = configureKalmanFilter('ConstantVelocity', ...
                centroid, [200, 50], [100, 25], 100);
            
            % create a new track
            newTrack = struct(...
                'id', nextId, ...
                'bbox', bbox, ...
                'kalmanFilter', kalmanFilter, ...
                'age', 1, ...
                'totalVisibleCount', 1, ...
                'consecutiveInvisibleCount', 0);
            
            % add it to the array of tracks
            tracks(end + 1) = newTrack;
            
            % increment the next id
            nextId = nextId + 1;
        end
    end

%% Display Tracking Results

    function displayTrackingResults()
        % convert the frame and the mask to uint8 RGB 
        frame = im2uint8(frame);
        mask = uint8(repmat(mask, [1, 1, 3])) .* 255;
        
        minVisibleCount = 8;
        if ~isempty(tracks)
              
            % noisy detections tend to result in short-lived tracks
            % only display tracks that have been visible for more than 
            % a minimum number of frames.
            reliableTrackInds = ...
                [tracks(:).totalVisibleCount] > minVisibleCount;
            reliableTracks = tracks(reliableTrackInds);
            
            % display the objects. If an object has not been detected
            % in this frame, display its predicted bounding box.
            if ~isempty(reliableTracks)
                % get bounding boxes
                bboxes = cat(1, reliableTracks.bbox);
                
                % get ids
                ids = int32([reliableTracks(:).id]);
                
                % create labels for objects indicating the ones for 
                % which we display the predicted rather than the actual 
                % location
                labels = cellstr(int2str(ids'));
                predictedTrackInds = ...
                    [reliableTracks(:).consecutiveInvisibleCount] > 0;
                isPredicted = cell(size(labels));
                isPredicted(predictedTrackInds) = {' predicted'};
                labels = strcat(labels, isPredicted);
                
                % draw on the frame
                frame = insertObjectAnnotation(frame, 'rectangle', ...
                    bboxes, labels);
                
                % draw on the mask
                mask = insertObjectAnnotation(mask, 'rectangle', ...
                    bboxes, labels);
            end
        end
        
        % display the mask and the frame
        obj.maskPlayer.step(mask);        
        obj.videoPlayer.step(frame);
    end

displayEndOfDemoMessage(mfilename)
end

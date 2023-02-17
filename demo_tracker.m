function demo_tracker()

% TODO: put name oy four tracker here
tracker_name = 'ms';
% TODO: select a sequence you want to test on
sequence = 'jogging';
% TODO: give path to the dataset folder
dataset_path = 'vot2014\';

use_reinitialization = true;
skip_after_fail = 5;

% specify initialize and update function
initialize = str2func(sprintf('%s_initialize', tracker_name));
update = str2func(sprintf('%s_update', tracker_name));

% read all frames in the folder
base_path = fullfile(dataset_path, sequence);
img_dir = dir(fullfile(base_path, '*.jpg'));

% read ground-truth
% bounding box format: [x,y,width, height]
gt = dlmread(fullfile(base_path, 'groundtruth.txt'));
if size(gt,2) > 4
    % ground-truth in format: [x0,y0,x1,y1,x2,y2,x3,y3], convert:
    X = gt(:,1:2:end);
    Y = gt(:,2:2:end);
    X0 = min(X,[],2);
    Y0 = min(Y,[],2);
    W = max(X,[],2) - min(X,[],2) + 1;
    H = max(Y,[],2) - min(Y,[],2) + 1;
    gt = [X0, Y0, W, H];
end

start_frame = 1; %580
n_failures = 0;

figure(1); clf;
frame = 1;
while frame <= numel(img_dir)

    % read frame
    img = imread(fullfile(base_path, img_dir(frame).name));

    if frame == start_frame
        % initialize tracker
        tracker = initialize(img, gt(frame,:));
        bbox = gt(frame, :);
    else
        % update tracker (target localization + model update)
        %tic()
        [tracker, bbox] = update(tracker, img);
        %toc()
    end
    
    % show image
    imshow(img);
    hold on;
    rectangle('Position',bbox, 'LineWidth',2, 'EdgeColor','y');
    % show current number of failures
    text(12, 15, sprintf('Failures: %d', n_failures), 'Color','w', ...
        'FontSize',10, 'FontWeight','normal', ...
        'BackgroundColor','k', 'Margin',1);
    % show frame number
    text(12, 35, sprintf('Frame: #%d', frame), 'Color','w', ...
        'FontSize',10, 'FontWeight','normal', ...
        'BackgroundColor','k', 'Margin',1);
    
    hold off;
    drawnow;
    
    % detect failures and reinit
    if use_reinitialization
        area = rectint(bbox, gt(frame,:));
        if area < eps && use_reinitialization
            disp('Failure detected. Reinitializing tracker...');
            frame = frame + skip_after_fail - 1;  % skip 5 frames at reinit (like VOT)
            start_frame = frame + 1;
            n_failures = n_failures + 1;
        end
    end
    
    frame = frame + 1;
    
end

end  % endfunction


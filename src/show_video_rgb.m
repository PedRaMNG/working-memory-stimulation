function [handle] = show_video(frames, properties)
    % show_video Calls the implay function, adjusts colormaps/ranges for sections,
    % and sets the frame rate (without attempting to close the error dialog).
    % Call it with 2 parameters:
    % show_video(frames, properties)
    % frames - 4D array of images (combined sections) [Height x Width x 3(RGB) x Time]
    % properties - struct containing video properties and section properties
    %   properties.video_fps - frame rate for the video
    %   properties.section_properties - cell array of structs, each defining
    %                                colormap and limits for a section.
    %                                Example:
    %                                properties.section_properties = {
    %                                    struct('cmap', cmap1, 'limits', limits1),
    %                                    struct('cmap', cmap2, 'limits', limits2),
    %                                    ...
    %                                };
    % Returns a handle to the player

    if nargin < 2
        properties = [];
    end
    if ~isfield(properties, 'file_name')
        properties.file_name = 'video';
    end
    if ~isfield(properties, 'video_fps')
        fps = 10;
    else
        fps = properties.video_fps;
    end
    if ~isfield(properties, 'frames_range')
        frames_range = [1, size(frames, 4)];
    else
        if (properties.frames_range(2) > size(frames, 4))
            frames_range(1) = 1;
            frames_range(2) = size(frames, 4);
        else
            frames_range = properties.frames_range;
        end
    end
    if ~isfield(properties, 'fit_to_window')
        properties.fit_to_window = true;
    end

    if ~isfield(properties, 'section_properties')
        % If section_properties is not provided, use default colormap and limits
        default_cmap = hot(256);
        default_cmap = default_cmap(end:-1:1, :);
        default_limits = [double(min(frames(:))), double(max(frames(:)))];
        properties.section_properties = cell(1, 6); % 6 sections
        for i = 1:6
            properties.section_properties{i} = struct('cmap', default_cmap, 'limits', default_limits);
        end
    end

    frames_video = frames(:, :, :, frames_range(1):frames_range(2));
    handle = implay(frames_video, fps);

    % Define manual window size and position
    if ~isfield(properties, 'figure_position')
        properties.figure_position = [100, 100, 800, 600]; % Default position [left, bottom, width, height]
    end

    % Apply the manual position unless fullscreen mode is enabled
    if properties.full_screen
        handle.Parent.WindowState = 'maximized';
    else
        handle.Parent.Position = properties.figure_position;
    end

    if properties.fit_to_window
        showHiddenHandles = get(0, 'showHiddenHandles');
        set(0, 'showHiddenHandles', 'on');
        ftw = handle.Parent.findobj('TooltipString', 'Maintain fit to window');
        if ~isempty(ftw)
            ftw.ClickedCallback();
        end
        set(0, 'showHiddenHandles', showHiddenHandles);
    end

    handle.Parent.Name = properties.file_name;
end

function toggle_GUI(handles, action, custom_CloseRequestFcn)
    if ~exist('custom_CloseRequestFcn', 'var')
        custom_CloseRequestFcn = @closereq;
    end

    % Read INFO
    INFO = getappdata(handles.figure,'INFO');
    
    if strcmp(action, 'list')
        INFO.ActionableInterfaceObjects = findobj(handles.figure, 'Enable', 'on');
        setappdata(handles.figure, 'INFO',INFO);
        return
    else
        set(INFO.ActionableInterfaceObjects, 'Enable', action);
    end
    
    % Adjust the close policy of the window
    if strcmp(action, 'on')
        handles.figure.CloseRequestFcn = custom_CloseRequestFcn;
    elseif strcmp(action, 'off')
        handles.figure.CloseRequestFcn = '';
    end

function brainomatic
% A simple tool for single channel EEG recording / display / comparison
%
% Peter Zeidman

% Number of frames to process at a time
frameLength = 20000;

% Size of buffer for recording
bufferLengthMins = 7;

% Default AM demodulation frequency (off by default)
am_freq = 1000;

reader        = []; % Reader object
Fs            = []; % Samples per second
signal_buffer = []; % All data received so far
fttx          = {}; % X-axis for each time window's FFT
ffts          = {}; % Data for each time window's FFT

% Default filtering options
am_enabled     = false;
zscore_enabled = true;

i = 1; % Audio time

% Prepare GUI
% -------------------------------------------------------------------------

% Load SPM
spm defaults eeg

% GUI states (stored in figure's Userdata)
GUI_RECORDING = 1;
GUI_STOPPED   = 0;
GUI_RESET_REQ = -1;

% Create Figure
fig = figure('Color','white','name','Brain-o-matic','MenuBar','none','NumberTitle','off');

% Top panel
p1 = uipanel(fig,'FontSize',12,...
             'BackgroundColor','white',...
             'Position',[0.05 .60 .9 .35]);

% Middle panel
p2 = uipanel(fig,'FontSize',12,...
             'BackgroundColor','white',...
             'Position',[0.05 .20 .9 .35]);

% Button panel
p3 = uipanel(fig,'FontSize',12,...
             'BackgroundColor','white',...
             'Position',[0.05 .05 .9 .10]);
         
% Plot axes         
h1 = axes('Parent',p1);
h2 = axes('Parent',p2);

% Buttons
uicontrol('Callback',@go,...
    'Style','pushbutton','String','Record / Play','Parent',p3,'Units','Normalized','Position',[0.2 0.1 0.18 0.8],'FontSize',12);
uicontrol('Callback',@stop,...
    'Style','pushbutton','String','Stop','Parent',p3,'Units','Normalized','Position',[0.4 0.1 0.18 0.8],'FontSize',12);
uicontrol('Callback',@reset,...
    'Style','pushbutton','String','Reset','Parent',p3,'Units','Normalized','Position',[0.6 0.1 0.18 0.8],'FontSize',12);

% Get list of audio device names
adi = audiodevinfo();
device_names = {adi.input.Name};
for d = 1:length(device_names)
    k = strfind(device_names{d},'(');
    if ~isempty(k)
        device_names{d} = strtrim(device_names{d}(1:k-1));
    end
end

% Input menu
m = uimenu('Label','Inputs');
uimenu(m,'Label','Load Wav','Callback',@set_input);
for d = 1:length(device_names)
    uimenu(m,'Label',device_names{d},'Callback',@set_input);
end

% Filtering menu
m = uimenu('Label','Filtering');
am_menu_item     = uimenu(m,'Label','AM demodulation','Callback',@set_am_filtering);
zscore_menu_item = uimenu(m,'Label','Z-score signal','Callback',@set_zscore,'Checked','On');

% Settings menu
m = uimenu('Label','Settings');
buff_menu_item = uimenu(m,'Label','Buffer length','Callback',@set_buffer_len);

% -------------------------------------------------------------------------
function set_input(varargin)
    % Callback for the set input dropdown selector
    
    selected_input = get(varargin{1},'Label');
    
    if strcmp(selected_input,'Load Wav')
        [FileName,PathName] = uigetfile('*.wav','Select the recording');
        
        if FileName == false
            return;
        end
        
        fn = fullfile(PathName,FileName);

        reader = dsp.AudioFileReader(fn,'SamplesPerFrame',frameLength);
        Fs = reader.SampleRate;
    else
        % A microphone was selected
        Fs = 10000;
        reader = dsp.AudioRecorder('SampleRate',Fs,'SamplesPerFrame',frameLength,'DeviceName',selected_input,'NumChannels',1);
    end

    % Reset audio timer
    i = 1;       

    % Initialize audio buffer
    buffer_length = (bufferLengthMins * 60) * Fs;
    signal_buffer = zeros(1,buffer_length);    
end

% -------------------------------------------------------------------------
function go(varargin)
    % Callback for Go button.
    
    % If already started then nothing to do
    if get(fig,'UserData') == GUI_RECORDING, return; end
    
    % If no input device selected, do nothing
    if isempty(reader)
        errordlg('Please select an input','No input');
        return;
    end

    % Set recording mode
    set(fig,'UserData', GUI_RECORDING);

    % Running sum of FFT
    fft_sum = [];
    
    while true;
        
        % Get data from input
        signal = reader;
        if isa(reader, 'dsp.AudioRecorder') || ...
                isa(reader,'dsp.AudioFileReader')
            signal = step(reader);
            signal = signal(:,1); % Left channel only
        end
        
        % AM demodulation if requested
        if am_enabled
            signal = amdemod(signal,am_freq,Fs);
        end
        
        % Normalize
        if zscore_enabled
            signal = (signal - mean(signal)) ./ std(signal);
        end
        
        % Append to signal buffer
        signal_buffer(i:i+frameLength-1) = signal;

        % Onset time of this frame
        secs = (i / 10000);

        % Plot timeseries
        axes(h1);
        x = (0:frameLength-1) ./ Fs;
        x = x + secs-1;
        plot(x,signal);
        xlabel('Secs','FontSize',12);
        ylabel('Signal','FontSize',12);    

        title(sprintf('Last 2 seconds of EEG signal (second %d)',round(secs)),'FontSize',16);
        set(gca,'YTickLabel',[]);

        % Update record of total samples
        i = i + frameLength;

        %if i < 100, continue; end

        % Plot any previous FFTs
        axes(h2);
        for j = 1:length(ffts)                
            plot(fttx{j},ffts{j},'Color',[0.7 0.7 0.7]);
            hold on;
        end

        % Compute new FFT (just this second)
        x = (0:(frameLength-1)) ./ Fs;
        y = signal';                
        [spectrum,ntaper,f] = ft_specest_mtmfft(y, x,...
            'taper','hanning','method','mtmfft','tapsmofrq',1);
        P1 = real(squeeze(spectrum));
        
        % Add to running sum
        if isempty(fft_sum)
            fft_sum = abs(P1);
        else
            fft_sum = fft_sum + abs(P1);
        end
        
        % Compute average FFT over windows
        fft_avg = fft_sum  ./ (round(i / frameLength));

        % Plot FFT avg
        plot(f,fft_avg,'Color',[0 0 0]);
        xlabel('Brain rhythm''s frequency (Hz)','FontSize',12);
        ylabel('Power','FontSize',12);
        set(gca,'YTickLabel',[]);
        title('Amount of each brain rhythm (whole recording)','FontSize',16);
        xlim([3 60]);
        hold off;

        % Refresh display
        drawnow();                 
        
        % If this is an audio file, flag to stop if finished
        if isa(reader,'dsp.AudioFileReader') && isDone(reader)
            set(fig,'UserData',GUI_STOPPED);
        end
        
        if get(fig,'UserData') == GUI_STOPPED
            % Save at end of recording
            L  = i - 1;        
            fttx{end+1} = f;
            ffts{end+1} = fft_avg;
            save_recording(signal_buffer);

            break;
        elseif get(fig,'UserData') == GUI_RESET_REQ
            % Reset has been requestsed
            break;
        end
        
    end

    if get(fig,'UserData') == GUI_RESET_REQ
        % Reset was requested - discard data
        perform_reset();
    else
        % Save at end of recording
        L  = i - 1;      
        fttx{end+1} = f;
        ffts{end+1} = fft_avg;
        save_recording(signal_buffer);
        set(fig,'UserData',GUI_STOPPED);
    end
end

% -------------------------------------------------------------------------
function save_recording(y)
    audiowrite(['recording_' char(datestr(now,'dd-mmmm-yyyy-HH-MM-SS')) '.wav'],y,Fs);
end

% -------------------------------------------------------------------------
function stop(varargin)
    set(fig,'UserData', 0);
end
% -------------------------------------------------------------------------
function reset(varargin)
    % Callback for reset button.
    
    % Confirm
    yesno = questdlg('Are you sure you wish to reset?',...
        'Confirm','Yes','No',struct('Default','No','Interpreter','none'));
    
    if ~strcmp(yesno,'Yes'), return; end
    
    if get(fig,'UserData') == GUI_RECORDING
        % Set a flag to stop the recording
        set(fig,'UserData', GUI_RESET_REQ);
    else
        % Stop it now
        perform_reset;
    end
    
end
% -------------------------------------------------------------------------
function perform_reset()
    %recordings = {};
    fttx = {};
    ffts = {};
    i = 1;
    axes(h1);
    cla;
    axes(h2);
    cla;
    set(fig,'UserData',GUI_STOPPED); 
end
% -------------------------------------------------------------------------
function set_am_filtering(varargin)
    % Callback for AM filtering menu item
    am_enabled = ~am_enabled;
    
    % Prompt for frequency
    if am_enabled
        str = inputdlg({'Enter frequency'},'AM demodulation',[1 50],{'1000'});
        if isempty(str)
            am_enabled = false;
            return;
        end
        am_freq = str2double(str);
    end
    
    % Update menu
    if am_enabled
        set(am_menu_item,'Checked','on');
    else
        set(am_menu_item,'Checked','off');
    end
end
% -------------------------------------------------------------------------
function set_zscore(varargin)
    zscore_enabled = ~zscore_enabled;
    
    if zscore_enabled
        set(zscore_menu_item,'Checked','on');
    else
        set(zscore_menu_item,'Checked','off');
    end
end

% -------------------------------------------------------------------------
function set_buffer_len(varargin)
    % Callback for the set buffer length menu item
    
    % Prompt for buffer length
    str = inputdlg('Set buffer length (mins):',...
        'Buffer Length',1,{num2str(bufferLengthMins)});
    
    if isempty(str), return, end
    
    % Prompt to confirm
    yesno = questdlg('All data will be cleared. Continue?',...
        'Yes','No',struct('Default','No','Interpreter','none')); 
    if ~strcmp(yesno,'Yes'), return; end
    
    perform_reset();
    bufferLengthMins = str2double(str);
end

end
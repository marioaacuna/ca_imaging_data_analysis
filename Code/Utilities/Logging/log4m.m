classdef log4m < handle 
    %LOG4M This is a simple logger based on the idea of the popular log4j.
    %
    % Description: Log4m is designed to be relatively fast and very easy to
    % use. It has been designed to work well in a matlab environment.
	% 
	% Description of logging levels:
	% TRACE: This level is mostly used to trace a code path.
	% DEBUG: This level is mostly used to log debugging output that may help
    %           identify an issue or verify correctness by inspection.
	% INFO: This level is intended for general user messages about the progress
    %           of the program.
	% WARNING (unrelated to Matlab's warning() function): This level is used to
    %           alert the user of a possible problem.
	% ERROR (unrelated to Matlab's error() function): This level is used for
    %           non-critical errors that can endanger correctness.
	% CRITICAL: This level is used for critical errors that definitely endanger
    %           correctness.
    
    properties (Constant)
        ALL = 0;
        TRACE = 1;
        DEBUG = 2;
        INFO = 3;
        WARN = 4;
        ERROR = 5;
        CRITICAL = 6;
        OFF = 7;
    end
        
    properties(Access = protected)
        logger;
        lFile;
    end
    
    properties(SetAccess = protected)
        fullpath = 'log4m.txt';  % Default file
        commandWindowLevel = log4m.ALL;
        logLevel = log4m.INFO;
        decorator = '#';
        last_message = '';
    end
    
    methods (Static)
        function obj = getLogger(logPath)
            %GETLOGGER Returns instance unique logger object.
            %   PARAMS:
            %       logPath - Relative or absolute path to desired logfile.
            %   OUTPUT:
            %       obj - Reference to signular logger object.
            %
            
            if(nargin == 0)
                logPath = 'log4m.log';
            elseif(nargin > 1)
                error('log4m:input_error', 'getLogger only accepts one parameter input');
            end
            
            localObj = log4m(logPath);
            obj = localObj;
        end
        
        function testSpeed(logPath)
            %TESTSPEED Gives a brief idea of the time required to log.
            %
            %   Description: One major concern with logging is the
            %   performance hit an application takes when heavy logging is
            %   introduced. This function does a quick speed test to give
            %   the user an idea of how various types of logging will
            %   perform on their system.
            %
            
            L = log4m.getLogger(logPath);
            
            disp('1e5 logs when logging only to command window');
            
            L.setCommandWindowLevel(L.TRACE);
            L.setLogLevel(L.OFF);
            tic;
            for i=1:1e5
                L.trace('log4mTest','test');
            end
            
            disp('1e5 logs when logging only to command window');
            toc;
            
            disp('1e6 logs when logging is off');
            
            L.setCommandWindowLevel(L.OFF);
            L.setLogLevel(L.OFF);
            tic;
            for i=1:1e6
                L.trace('log4mTest','test');
            end
            toc;
            
            disp('1e4 logs when logging to file');
            
            L.setCommandWindowLevel(L.OFF);
            L.setLogLevel(L.TRACE);
            tic;
            for i=1:1e4
                L.trace('log4mTest','test');
            end
            toc;
            
        end
    end
    
    
%% Public Methods Section %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods       
        function setFilename(self, logPath)
            %SETFILENAME Change the location of the text log file.
            %
            %   PARAMETERS:
            %       logPath - Name or full path of desired logfile
            %
            
            [fid,message] = fopen(logPath, 'a');
            
            if(fid < 0)
                error('log4m:file_error', ['Problem with supplied logfile path: ' message]);
            end
            fclose(fid);
            
            self.fullpath = logPath;
        end
          
     
        function setCommandWindowLevel(self,loggerIdentifier)
            self.commandWindowLevel = loggerIdentifier;
        end


        function setLogLevel(self,logLevel)
            self.logLevel = logLevel;
        end
        
        function setDecorator(self,decorator)
            self.decorator = decorator;
        end
        
        function setLastMessage(self,last_message)
            self.last_message = last_message;
        end
        

%% The public Logging methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function trace(self, message, varargin)
            self.writeLog(self.TRACE, message, varargin{:});
        end
        
        function debug(self, message, varargin)
            self.writeLog(self.DEBUG, message, varargin{:});
        end
        
        function info(self, message, varargin)
            self.writeLog(self.INFO, message, varargin{:});
        end

        function warn(self, message, varargin)
            self.writeLog(self.WARN, message, varargin{:});
        end
        
        function error(self, message, varargin)
            self.writeLog(self.ERROR, message, varargin{:});
        end

        function critical(self, message, varargin)
            self.writeLog(self.CRITICAL, message, varargin{:});
        end
        
    end

%% Private Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        
        function self = log4m(fullpath_passed)
            
            if(nargin > 0)
                path = fullpath_passed;
            end
			self.setFilename(path);
        end

%% WriteToFile        
        function writeLog(self, level, message, varargin)
            % Parse inputs
            if isempty(varargin)
                decorate               = false;
                contains_path          = false;
                write_to_file          = true;
                print_on_screen        = true;
                overwrite_last_message = false;
                append                 = false;
            else
                p = inputParser();
                addParameter(p, 'decorate',               false, @islogical)
                addParameter(p, 'contains_path',          false, @islogical)
                addParameter(p, 'write_to_file',          true,  @islogical)
                addParameter(p, 'print_on_screen',        true,  @islogical)
                addParameter(p, 'overwrite_last_message', false, @islogical)
                addParameter(p, 'append',                 false, @islogical)
                parse(p, varargin{:})
                decorate = p.Results.decorate;
                contains_path = p.Results.contains_path;
                write_to_file = p.Results.write_to_file;
                print_on_screen = p.Results.print_on_screen;
                overwrite_last_message = p.Results.overwrite_last_message;
                append = p.Results.append;
            end
            
            % Escape path, if any
            if contains_path
                message = strrep(message, '\', '\\');
            end
            % Escape percent symbol
            message = strrep(message, '%', '%%');
            
            % Get name of function caller
            scriptName = self.get_caller_name();

            % If necessary write to command window
            if self.commandWindowLevel <= level && print_on_screen
                % Make message to print on screen
                if decorate
                    message_length = length(message);
                    decoration = repmat(self.decorator, 1, message_length);
                    command_window_message = ['\n\n', decoration, '\n', message, '\n', decoration];
                else
                    command_window_message = message;
                end
                % Append new line
                command_window_message = [command_window_message, '\n'];
                % Delete last message, if necessary
                if overwrite_last_message
                    fprintf(repmat('\b', 1, length(self.last_message)-1))
                elseif append
                    fprintf('\b ')  % Backspace and space
                end
                % Print new message
                fprintf(command_window_message);
            else
                command_window_message = '';
            end
            % Keep last message in memory
            self.setLastMessage(command_window_message)
            
            %I f currently set log level is too high, just skip this log
            if self.logLevel > level || ~write_to_file
                return
            end
            
            % Set up our level string
            switch level
                case self.TRACE
                    levelStr = 'TRACE';
                case self.DEBUG
                    levelStr = 'DEBUG';
                case self.INFO
                    levelStr = 'INFO';
                case self.WARN
                    levelStr = 'WARN';
                case self.ERROR
                    levelStr = 'ERROR';
                case self.CRITICAL
                    levelStr = 'CRITICAL';
                otherwise
                    levelStr = 'UNKNOWN';
            end

            % Append new log to log file
            if decorate
                decoration = repmat(self.decorator, 1, 3);
                log_file_message = [decoration, ' ', message, ' ', decoration];
            else
                log_file_message = message;
            end

            try
                fid = fopen(self.fullpath,'a');
                fprintf(fid, '%-23s %-8s %-s %s\r\n', datestr(now,'yyyy-mm-dd HH:MM:SS.FFF'), levelStr, scriptName, log_file_message);
                fclose(fid);
            catch ME
                if strcmp(ME.identifier, 'MATLAB:FileIO:InvalidFid')
                    % Print on screen only
                    disp(log_file_message)
                else
                    % Rethrow exception
                    rethrow(ME)
                end
            end 
        end
    end
    
    %% Get name of caller function
    methods(Access = private, Static)
        function funcName = get_caller_name()
            % Get names of functions in the stack
            FunStack = dbstack('-completenames');
            % Caller function is the 4th one above this one. If not available,
            % return something else
            if length(FunStack) >= 4
                funcName = FunStack(4).name;
            else
                funcName = 'SCRIPT';
            end
        end
    end
end


classdef tbanz 
    %V1ANZ Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        path_raw_data = tbanz.get_path_raw_data() 
    end
    
    methods (Static)
        

        

        %% - - - Locate Data Storage Folder
        function path_raw_data = get_path_raw_data(obj) 
            if ismac && strcmp(getenv('LOGNAME'), 'luis')
                path_raw_data = '/Users/luis/Box/boxPHD/Toolbox/2pdata/';
            elseif ispc && strcmp(getenv('username'), 'Luis')
                path_raw_data = 'C:\Users\Luis\Box Sync\boxEXPER\Toolbox\2pdata\';
            elseif ispc && strcmp(getenv('username'), 'dario')
                path_raw_data = 'C:\2pdata\';
            else
                error('unknown computer')
            end   
            obj.path_raw_data = path_raw_data;
            
        end 
        
       

        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end


function data = loadConf( base_path, conf_name, frameIDs )
    fullfile( base_path, conf_name );
    fid = fopen( fullfile( base_path, conf_name ), 'r' );
    
    if fid == -1
        data=[];
        return;
    end
    
    tline = fgetl(fid);

    frameCount = 0;
    frameIDs = [];
    cnt = 0;
    while ischar(tline)
      % tokenize line
      tokens = strsplit( tline );

      if strcmp( tokens{1}, 'dataset' )
        data.dataset = tokens{2};
      end

      if strcmp( tokens{1}, 'color_resolution' )
        data.color_resolution = [str2num(tokens{2}), str2num(tokens{3})];
      end
      
      if strcmp( tokens{1}, 'depth_resolution' )
        data.depth_resolution = [str2num(tokens{2}), str2num(tokens{3})];
      end
      
      %decide what to do with the line based on first token
      if strcmp( tokens{1}, 'intrinsics' )
        data.intrinsics_name = tokens{2};
        data.K = reshape( readValuesFromTxt( fullfile( base_path, tokens{2}) ), 3, 3 )';
      end

      if strcmp( tokens{1}, 'depth_intrinsics' )
        data.intrinsics_name = tokens{2}
        data.K = reshape( readValuesFromTxt( fullfile( base_path, tokens{2}) ), 3, 3 )';
      end

      if strcmp( tokens{1}, 'n_images' )
        frameCount = str2num( tokens{2} );
        frameIDs = 1:frameCount;
        data.n_images = frameCount;
      end

      if strcmp(tokens{1}, 'image_directory')
        data.image_directory = tokens{2};
      end
      
       if strcmp(tokens{1}, 'color_directory')
        data.image_directory = tokens{2};
      end

      if strcmp(tokens{1}, 'depth_directory')
        data.depth_directory = tokens{2} ;
      end

      if strcmp(tokens{1}, 'scan')
        cnt = cnt + 1;
        
        depth_name = tokens{2};
        image_name = tokens{3};
        
        image_path = fullfile( base_path, data.image_directory, image_name );
        depth_path = fullfile( base_path, data.depth_directory, depth_name );
        
        if ~exist(image_path, 'file')
            continue;
        end;
        if ~exist(depth_path, 'file')
            continue;
        end;
        
        data.image{cnt} = fullfile( base_path, data.image_directory, image_name );
        data.depth{cnt} = fullfile( base_path, data.depth_directory, depth_name );
      end

      tline = fgetl(fid);
    end

    fclose(fid);
     
    if ~isfield(data, 'depth_resolution')
       data.depth_resolution = [640 480];
    end;
    if ~isfield(data, 'color_resolution')
       data.color_resolution = [640 480];
    end;
end

%% IO function
function values = readValuesFromTxt(filename)
    try
        values = textscan(urlread(filename),'%f');
    catch
        fid = fopen(filename,'r');
        values = textscan(fid,'%f');
        fclose(fid);
    end
    values = values{1};
end


function files = dirSmart(page, tag)
    [files, status] = urldir(page, tag);
    if status == 0
        files = dir(fullfile(page, ['*.' tag]));
    end
end


function fileStr = file2string(fname)
    fileStr = '';
    fid = fopen(fname,'r');
    tline = fgetl(fid);
    while ischar(tline)
        fileStr = [fileStr ' ' tline];
        tline = fgetl(fid);
    end
    fclose(fid);
end


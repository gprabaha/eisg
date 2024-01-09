function p = project_path()
p = fileparts( fileparts(fileparts(which('eisg.util.project_path'))) );
end
fileID = fopen('job_list.txt','w');
indices = floor( linspace( 1, 574, 100 ) );

for i = 1:numel( indices ) - 1
    file_path = "'/gpfs/milgram/pi/chang/pg496/repositories/eisg/scripts/run_sf_coherence_in_cluster.m'";
    if i==1
        matlab_command = sprintf('matlab -nodisplay -nosplash -r "to_process_inds=%d:%d; run(%s);"', indices(i), indices(i+1), file_path );
    else
        matlab_command = sprintf('matlab -nodisplay -nosplash -r "to_process_inds=%d:%d; run(%s);"', indices(i)+1, indices(i+1), file_path );
    end
    fcat_path = 'export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/gpfs/milgram/project/chang/pg496/repositories/categorical/lib/linux"';
    matlab_module = 'module load MATLAB/2020a';
    fprintf( fileID, "%s; %s; %s;\n", fcat_path, matlab_module, matlab_command );
end
fclose( fileID );
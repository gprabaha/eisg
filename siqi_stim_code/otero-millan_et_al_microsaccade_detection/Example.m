%%
addpath( '/Users/nick/Downloads/otero-millan_et_al_microsaccade_detection' );

% this script shows how to use the microsaccade detection method
% published in Otero-Millan et al. Journal of Vision 2014

% Set up variables --------------------------------------------------------
folder = '/Users/nick/Downloads/otero-millan_et_al_microsaccade_detection/data';


% =========================================================================
% EDIT THIS ===============================================================
% Make sure to fill the variables samples with the eye movement data in the
% proper colums. In blinks mark all the samples that are not good data. In
% most VOG systems it is best to remove 100 ms before and after the blink
% to remove all artifacts.
% =========================================================================
session = 'test';
samplerate = 1e3;
    
int_p = '/Volumes/external3/data/changlab/siqi/stim/intermediates';
pos_file = shared_utils.io.fload( fullfile(int_p, 'aligned_raw_samples/position/01222020_position_9.mat') );
m1_p = pos_file.m1;

fs = 1e3;
m1_p_deg = bfw.px2deg( m1_p );
samples = [];
samples(:, 1) = 1:size(m1_p_deg, 2);
samples(:, 2:3) = m1_p_deg';
samples(:, 4:5) = m1_p_deg';
blinks = false( size(m1_p_deg, 2), 1 );


% Loads the recording and prepares it por processing
recording = ClusterDetection.EyeMovRecording.Create(folder, session, samples, blinks, samplerate);

% Runs the saccade detection
[saccades stats] = recording.FindSaccades();

% Plots a main sequence
enum = ClusterDetection.SaccadeDetector.GetEnum;
figure
subplot(2,2,1)
plot(saccades(:,enum.amplitude),saccades(:,enum.peakVelocity),'o')
set(gca,'xlim',[0 1],'ylim',[0 100]);
xlabel('Saccade amplitude (deg)');
ylabel('Saccade peak velocity (deg/s)');


% Plots the traces with the labeled microsaccades
subplot(2,2,[3:4])
plot(samples(:,1), samples(:,2:end));
hold
yl = get(gca,'ylim');
u1= zeros(size(samples(:,1)))+yl(1);
u2= zeros(size(samples(:,1)))+yl(1);
u1((saccades(:,enum.startIndex))) = yl(2);
u2(saccades(:,enum.endIndex)) = yl(2);
u = cumsum(u1)-cumsum(u2);
plot(samples(:,1), u,'k')

xlabel('Time (ms)');
ylabel('Eye Position (deg)');

legend({'Left Horiz', 'Left Vert', 'Right Horiz' , 'Right Vert', 'Microsaccades'})

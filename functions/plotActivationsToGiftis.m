%% Plot activations on gifti and inflated gifti, at lateral, medial, and inferior angles
%
%   INPUTS
%       elecs = electrodes table of real electrodes, or give empty if you don't wish to plot
%       elecsAct = "electrodes table" format for the positions activations are at (i.e., bipolar-coordinate electrodes table)
%       giis = pial giftis, matching hemis
%       giiInfs = inflated pial giftis, matching hemis
%       sulcs = sulcus objects, matching hemis
%       hemis = cell array of hemispheres in subject
%       activations = vector of activations that matches electrodes
%       wm = weight max
%       outdir = output directory to save plots
%       activationStr = string of the type of activation being plotted (e.g., BBrsq)
%       visi =   Visibility of plots as they appear
%
function plotActivationsToGiftis(elecs, elecsAct, activations, giis, giiInfs, sulcs, hemis, wm, outdir, activationStr, visi)

    % View giftis and inflated gifti from 3 different angle. Configure theta and phi for each
    angleText = {'lateral', 'medial', 'inferior'};
    
    for ii = 1:length(hemis)
        
        % theta and phi, configured for each hemisphere
        ths = [choose(strcmpi(hemis{ii}, 'r'), 90, -90), choose(strcmpi(hemis{ii}, 'r'), -90, 90), choose(strcmpi(hemis{ii}, 'r'), 90, -90)];
        phis = [0, 0, -90];
    
        elecsActHemi = elecsAct(strcmpi(hemis{ii}, elecsAct.hemisphere), :);
        activationsHemi = activations(strcmpi(hemis{ii}, elecsAct.hemisphere));
    
        % normal gifti
        figure('Position', [200, 200, 600, 400], 'Visible', visi); ieeg_RenderGifti(giis{ii}); alpha 0.3
        hold on
        
        % isolate the electrodes and activation values in current hemi and plot
        if ~isempty(elecs)
            elecsHemi = elecs(strcmpi(hemis{ii}, elecs.hemisphere), :);
            plot3(elecsHemi.x, elecsHemi.y, elecsHemi.z, 'o', 'MarkerSize', 4, 'Color', 'k', 'MarkerFaceColor', 'w');
        end
    
        % plot weighted activations on gifti
        plotCortexWeightsAdjustable([elecsActHemi.x, elecsActHemi.y, elecsActHemi.z], activationsHemi, 'SigFraction', 0.05, 'MaxWt', wm, 'SizeJump', 2);
        hold off
    
        % rotate camera to each of 3 different views
        for jj = 1:3
            ieeg_viewLight(ths(jj), phis(jj));
            saveas(gcf, fullfile(outdir, sprintf('gifti%s_%s_%s.png', upper(hemis{ii}), activationStr, angleText{jj})));
        end
        close(gcf);
    
    
        % for inflated gifti, a new figure needs to be made for each angle, because opaque and els_popout
        for jj = 1:3
            xyzInf = els_popout([elecsActHemi.xInf, elecsActHemi.yInf, elecsActHemi.zInf], ths(jj), phis(jj), 0.02);
            figure('Position', [200, 200, 600, 400], 'Visible', visi); ieeg_RenderGifti(giiInfs{ii}, sulcs{ii});
            hold on
            
            % don't plot actual electrodes themselves on inflated brain (too busy)
    
            plotCortexWeightsAdjustable(xyzInf, activationsHemi, 'SigFraction', 0.05, 'MaxWt', wm, 'SizeJump', 2);
            hold off
            ieeg_viewLight(ths(jj), phis(jj));
            saveas(gcf, fullfile(outdir, sprintf('giftiInf%s_%s_%s.png', upper(hemis{ii}), activationStr, angleText{jj})));
            close(gcf);
        end
    
    end
end
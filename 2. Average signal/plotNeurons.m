function plotNeurons(av_neuron, startNeuron, endNeuron, figureTitle, xLimit)
    numNeurons = endNeuron - startNeuron + 1;

    fontSize = 13;  % Increase this value if you want bigger text

    for i = 1:numNeurons
        subplot(numNeurons, 1, i)
        plot(av_neuron(:, startNeuron + i - 1))

        % Set y-axis limits
        % ylim([-2, 2]) % for cond 1
        % ylim([-6, 6]) % for cond 2/3
        
        % Set x-axis limits
        xlim([0, xLimit])
        
        % Draw a red vertical line across the entire y-axis range
        % line([20, 20], [-0.06, 0.06], 'Color', 'r')
        vline(20, 'r');
        
        title(['SC #', num2str(startNeuron + i - 1)], 'FontSize', fontSize)
        
        % Turn off x-axis label for all but the last subplot
        if i ~= numNeurons
            set(gca, 'XTickLabel', []);
        end
        
        % Set font size for tick marks
        ax = gca;
        ax.FontSize = fontSize;
    end
    
    % Add a figure title
    sgtitle(figureTitle, 'FontSize', fontSize + 2);
    
    % Add centralized x and y labels
    h = axes('visible', 'off'); % Invisible axes
    h.XLabel.Visible = 'on';
    h.YLabel.Visible = 'on';
    xlabel('Frame #', 'FontSize', fontSize);
    ylabel('% ΔF/F', 'FontSize', fontSize);
end




%% Old function without centralized labels

% function plotNeurons(av_neuron, startNeuron, endNeuron, figureTitle, xLimit)
%     numNeurons = endNeuron - startNeuron + 1;
% 
%     for i = 1:numNeurons
%         subplot(numNeurons, 1, i)
%         plot(av_neuron(:, startNeuron + i - 1))
% 
%         % Set y-axis limits
%         % ylim([-0.02, 0.02])
%         % ylim([-0.06, 0.06])
% 
%         % Set x-axis limits
%         xlim([0, xLimit])
% 
%         % Draw a red vertical line across the entire y-axis range
%         % line([20, 20], [-0.06, 0.06], 'Color', 'r')
%         vline(20, 'r');
% 
%         title(['SC #', num2str(startNeuron + i - 1)])
%     end
%     % Add a figure title
%     sgtitle(figureTitle);
% end

% Without using the function:
% for i=1:10
%     subplot(10,1,i)
%     plot(av_neuron(:,i))
%     vline(20, 'r')
%     title(['Neuron #', num2str(i)])
% end

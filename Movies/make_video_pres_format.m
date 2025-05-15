function make_video_pres_format(title,X,Y1,T1,Y2,T2,params)

    % X is x axis
    % Y1 is no damage variable
    % T1 is time for no damage variable
    % Y2 is damage variable
    % T2 is time for damage variable

    % params for formatting:
    % - params.ylimit
    % - params.ytick
    
    % Open video writer
    v = VideoWriter(title,'MPEG-4');
    open(v);
    
    % Create time frames
    t = 0:T1(end)/2000:T1(end);

    % Assign variables
    if Y2 == 0
        x = X; y1 = interp1(T1,Y1,t);
    else
        x = X; y1 = interp1(T1,Y1,t); y2 = interp1(T2,Y2,t);
    end
    
    % Colour scheme
    % c1 = [0.4693    0.7545    0.5077]; % green
    c2 = [0.97 0.29 0.11]; % red
    % c3 = "#1B97E6"; % blue
    % c4 = "#FFCE3C"; % yellow
    
    % Create figure and set dimensions
    h = figure;
    % Set figure total dimension
    set(h,'Units','centimeters')
    % Absolute print dimensions of figure. 
    % [pos_from_left, pos_from_bottom, fig_width, fig_height]
    set(h,'Position',[0 0 10 13.5])
    ax = axes('Units','centimeters','InnerPosition',[1.5 2 8 11]);
    
    % Run through time
    for k = 1:length(t)/4
    
        % Plot
        hold(ax,'off')

        if Y2 == 0
            plot(x,y1(k,:),'LineWidth',1,'Color', "k");
        else
            plot(x,y1(k,:),'LineWidth',1,'Color', "k");
            hold(ax,'on')
            plot(x,y2(k,:),'LineWidth',1,'Color',c2);
        end
    
        % Formatting
        box(ax,"on")
        set(ax,'FontName','Times','FontSize',20);
        set(gcf,'Color','white')
    
        % Limits
        ylim(params.ylimit) % flux
        yticks(params.ytick)
        xlim([0 1])
        xticks([0 0.2 0.4 0.6 0.8 1])
    
        % Labels
        % ylabel('$\frac{\partial U}{\partial Z}$','Interpreter','latex','Rotation',0,'FontSize',30);
        xlabel('$Z$','Interpreter','latex','Rotation',0,'FontSize',20)
    
        % Legend
        legend('No damage','Damage','Interpreter','latex','Location','NorthEast','NumColumns',2,'Box','off','FontSize',15);
    
    
        % Save frame
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    
    close(v);

end
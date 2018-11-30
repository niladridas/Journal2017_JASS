% Plotting function
function plot_ot_enkf(samples, rep, x_est_OT,x_est_enkf,date_mee,time_pt_days,mode)
    stdy = 15;
     if mode == '1'
        % Pick a single run and plot
        ylabel_names = ['p','f','g','h','k','L'];
        measured_output = date_mee(stdy:end,7:12);
        x_est_OT = x_est_OT(stdy:end,:,1);
        x_est_enkf  = x_est_enkf(stdy:end,:,1);
        t_spot = time_pt_days(stdy:end);
        figure
        for i = 1:6
            subplot(3,2,i);
            a = plot(t_spot,x_est_OT(:,i));
            a.Color = 'red';
            a.LineStyle = '-';
            a.LineWidth = 1;
            a.Marker = 'o';
            a.MarkerSize = 1.5;
            a.MarkerEdgeColor = 'red';
            a.MarkerFaceColor = 'red';
            xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
            ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
            hold on
            c = plot(t_spot,measured_output(:,i));
            c.Color = [0 0.5 0];
            c.LineStyle = '-';
            c.LineWidth = 1;
            c.Marker = 'o';
            c.MarkerSize = 1.5;
            c.MarkerEdgeColor = [0 0.5 0];
            legend('OT','Obs')     
        end
        title_str = strcat('Plot of OT predictions for',{' '},num2str(samples),' samples');
        axes( 'Position', [0, 0.95, 1, 0.05] ) ;
        set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
        text( 0.55, 0, title_str, 'FontSize', 15', 'FontWeight', 'Bold', ...
          'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
        file_name = strcat('samplesOT',num2str(samples));
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(file_name,'-dpdf','-fillpage');
        figure
        for i = 1:6
            subplot(3,2,i);
            a = plot(t_spot,x_est_enkf(:,i));
            a.Color = 'red';
            a.LineStyle = '-';
            a.LineWidth = 1;
            a.Marker = 'o';
            a.MarkerSize = 2;
            a.MarkerEdgeColor = 'red';
            a.MarkerFaceColor = 'red';
            xlabel('Time (Days)','FontName','sans-serif','FontSize',10,'FontWeight','Bold');
            ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',10,'FontWeight','Bold');
        %     set(gca,'fontWeight','bold','Xtick',0:15);
            hold on
            c = plot(t_spot,measured_output(:,i));
            c.Color = [0 0.5 0];
            c.LineStyle = '-';
            c.LineWidth = 1;
            c.Marker = 'o';
            c.MarkerSize = 2;
            c.MarkerEdgeColor = [0 0.5 0];
            legend('EnKF','Obs')     
        %     grid on
        end
        title_str = strcat('Plot of EnKF predictions for',{' '},num2str(samples),' samples');
        axes( 'Position', [0, 0.95, 1, 0.05] ) ;
        set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
        text( 0.55, 0, title_str, 'FontSize', 15', 'FontWeight', 'Bold', ...
              'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
        file_name = strcat('samplesEnKF',num2str(samples));
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        print(file_name,'-dpdf','-fillpage');
    end
    
    
    
    if mode == '2'
        % Show the statistics of the prediction
        x_enkf_mean = mean(x_est_enkf(stdy:end,:,:),3);
        x_enkf_var =  sqrt(var(x_est_enkf(stdy:end,:,:),[],3));
        x_ot_mean = mean(x_est_OT(stdy:end,:,:),3);
        x_ot_var =  sqrt(var(x_est_OT(stdy:end,:,:),[],3));
        t_spot = time_pt_days(stdy:end);
        ylabel_names = ['p','f','g','h','k','L'];
        measured_output = date_mee(stdy:end,7:12);
        figure
        for i = 1:6
            subplot(3,2,i);
            %     a = plot(t_spot,x_est_OT(:,i));
            a =errorbar(t_spot,x_ot_mean(:,i),x_ot_var(:,i));
            a.Color = 'red';
            a.LineStyle = '-';
            a.LineWidth = 1;
            a.Marker = 'o';
            a.MarkerSize = 2;
            a.MarkerEdgeColor = 'red';
            a.MarkerFaceColor = 'red';
            xlabel('Time (Days)','FontName','sans-serif','FontSize',15,'FontWeight','Bold');
            ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',15,'FontWeight','Bold');
            %     set(gca,'fontWeight','bold','Xtick',0:15);
            hold on
            c = plot(t_spot,measured_output(:,i));
            c.Color = [0 0.5 0];
            c.LineStyle = '-';
            c.LineWidth = 1;
            c.Marker = 'o';
            c.MarkerSize = 2;
            c.MarkerEdgeColor = [0 0.5 0];
            legend('OT','True')     
        end
        title_str = strcat('Plot of variation of OT predictions for',{' '},num2str(samples),' samples');
        axes( 'Position', [0, 0.95, 1, 0.05] ) ;
        set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
        text( 0.55, 0, title_str, 'FontSize', 15, 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    
    
        figure
        for i = 1:6
            subplot(3,2,i);
            a =errorbar(t_spot,x_enkf_mean(:,i),x_enkf_var(:,i));
            a.Color = 'red';
            a.LineStyle = '-';
            a.LineWidth = 1;
            a.Marker = 'o';
            a.MarkerSize = 2;
            a.MarkerEdgeColor = 'red';
            a.MarkerFaceColor = 'red';
            xlabel('Time (Days)','FontName','sans-serif','FontSize',15,'FontWeight','Bold');
            ylabel(ylabel_names(i),'FontName','Times New Roman','FontSize',15,'FontWeight','Bold');
            %     set(gca,'fontWeight','bold','Xtick',0:15);
            hold on
            c = plot(t_spot,measured_output(:,i));
            c.Color = [0 0.5 0];
            c.LineStyle = '-';
            c.LineWidth = 1;
            c.Marker = 'o';
            c.MarkerSize = 2;
            c.MarkerEdgeColor = [0 0.5 0];
            legend('EnKF','True')     
        end
        title_str = strcat('Plot of variation of EnKF predictions for',{' '},num2str(samples),' samples');
        axes( 'Position', [0, 0.95, 1, 0.05] ) ;
        set( gca, 'Color', 'None', 'XColor', 'None', 'YColor', 'None' ) ;
        text( 0.55, 0, title_str, 'FontSize', 15', 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
    end
end
%% Pathnames and setup

% Path
BaseDir = '/home/BCBL/public/Lexiko/Lexiko2/Ave/Wrd/';

% subject-specific paths
SbjName = {'lx01',...
           'lx02',...
           'lx03',...
           'lx04',...
           'lx05',...
           'lx06',...
           'lx07',...
           'lx08',...
           'lx09',...
           'lx10',...
           'lx11',...
           'lx12',...
           'lx13',...
           'lx14',...
           'lx15',...
           'lx16' };

       
Cnd = {'ES',...
       'EU'};

   
% Read data into a cell array
sbj_idx = [1:10 12:16];
dip_idx = 1:2;



clear d y;


for k=1:15,
    for cnd=1:2,
            
        sbj = sbj_idx(k);
        
        % Set the path name
        PathName = strcat(BaseDir, SbjName{sbj}, '_word_', Cnd{cnd},'.fif');
     
        % Import the file
        d{sbj}{cnd} = fiff_read_dipwave(PathName);
        
        % Make array
        y(sbj,cnd,1:4,:) = d{sbj}{cnd}.evoked(:,:);

        % Get number of time points
        NPnts(sbj,cnd) = size(d{sbj}{cnd}.evoked.times,2);
        
        % Copy to simple array and scale to nA
        if NPnts(sbj,cnd) == 249,
            yy(sbj,cnd,1:4,:) = y(sbj,cnd).epochs(1:4,:).*1e9; 
        else
            for p=1:4,
                yy(sbj,cnd,p,:) = resample(y(sbj,cnd).epochs(p,:).*1e9, 249, NPnts(sbj,cnd)); 
            end
        end;
    end;
end;


y = yy;




% Error-band plot

N = 15;
FntSz = 18;

es.time = d{1}{1}.evoked.times;
eu.time = d{1}{1}.evoked.times;

figure;
subplot(2,2,1);

dipole = 1;
es.se = squeeze(mean(std(y(:,1,dipole,:),0,1)./sqrt(N)))';
eu.se = squeeze(mean(std(y(:,2,dipole,:),0,1)./sqrt(N)))';

es.mn = squeeze(mean(y(:,1,dipole,:),1));
eu.mn = squeeze(mean(y(:,2,dipole,:),1));

hold on;
ciplot(eu.mn+eu.se, eu.mn-eu.se, eu.time, [0.00 0.75 0.00]);
alpha(0.25);
ciplot(es.mn+es.se, es.mn-es.se, es.time, [0.75 0.00 0.00]);
alpha(0.5);

plot(es.time, es.mn, 'Color', [0.75 0.0  0.0], 'LineWidth', 3)
plot(eu.time, eu.mn, 'Color', [0.0  0.75 0.0], 'LineWidth', 3)
title('LH Temporal', 'FontSize', FntSz+2);

% plot(es.time(index), -2.5*(index(index)),'*','MarkerSize',12);

hold off;
axis tight; axis square; 
xlabel('s', 'FontSize', FntSz); ylabel('nA', 'FontSize', FntSz);

legend({'ES','EU'},'Location','NorthWest', 'FontSize', FntSz);
legend boxoff;

 
subplot(2,2,2)

dipole = 2;
es.se = squeeze(mean(std(y(:,1,dipole,:),0,1)./sqrt(N)))';
eu.se = squeeze(mean(std(y(:,2,dipole,:),0,1)./sqrt(N)))';

es.mn = squeeze(mean(y(:,1,dipole,:),1));
eu.mn = squeeze(mean(y(:,2,dipole,:),1));


hold on;
ciplot(eu.mn+eu.se, eu.mn-eu.se, eu.time, [0.00 0.75 0.00]);
alpha(0.25);
ciplot(es.mn+es.se, es.mn-es.se, es.time, [0.75 0.00 0.00]);
alpha(0.5);

plot(es.time, es.mn, 'Color', [0.75 0.0  0.0], 'LineWidth', 3)
plot(eu.time, eu.mn, 'Color', [0.0  0.75 0.0], 'LineWidth', 3)
title('RH Temporal', 'FontSize', FntSz+2);

% plot(es.time(index), -2*(index(index)),'*','MarkerSize',12);

hold off;
axis tight; axis square;
xlabel('s', 'FontSize', FntSz); ylabel('nA', 'FontSize', FntSz);


subplot(2,2,3)

dipole = 3;
es.se = squeeze(mean(std(y(:,1,dipole,:),0,1)./sqrt(N)))';
eu.se = squeeze(mean(std(y(:,2,dipole,:),0,1)./sqrt(N)))';

es.mn = squeeze(mean(y(:,1,dipole,:),1));
eu.mn = squeeze(mean(y(:,2,dipole,:),1));


hold on;
ciplot(eu.mn+eu.se, eu.mn-eu.se, eu.time, [0.00 0.75 0.00]);
alpha(0.25);
ciplot(es.mn+es.se, es.mn-es.se, es.time, [0.75 0.00 0.00]);
alpha(0.5);

plot(es.time, es.mn, 'Color', [0.75 0.0  0.0], 'LineWidth', 3)
plot(eu.time, eu.mn, 'Color', [0.0  0.75 0.0], 'LineWidth', 3)
title('LH Various', 'FontSize', FntSz+2);

hold off;
axis tight; axis square;
xlabel('s', 'FontSize', FntSz); ylabel('nA', 'FontSize', FntSz);

 
subplot(2,2,4)

dipole = 4;
es.se = squeeze(mean(std(y(:,1,dipole,:),0,1)./sqrt(N)))';
eu.se = squeeze(mean(std(y(:,2,dipole,:),0,1)./sqrt(N)))';

es.mn = squeeze(mean(y(:,1,dipole,:),1));
eu.mn = squeeze(mean(y(:,2,dipole,:),1));

hold on;
ciplot(eu.mn+eu.se, eu.mn-eu.se, eu.time, [0.00 0.75 0.00]);
alpha(0.25);
ciplot(es.mn+es.se, es.mn-es.se, es.time, [0.75 0.00 0.00]);
alpha(0.5);

plot(es.time, es.mn, 'Color', [0.75 0.0  0.0], 'LineWidth', 3)
plot(eu.time, eu.mn, 'Color', [0.0  0.75 0.0], 'LineWidth', 3)
title('RH Various', 'FontSize', FntSz+2);

% plot(es.time(index2), -5*(index2(index2)),'*','MarkerSize',12);

hold off;
axis tight; axis square;
xlabel('s', 'FontSize', FntSz); ylabel('nA', 'FontSize', FntSz);







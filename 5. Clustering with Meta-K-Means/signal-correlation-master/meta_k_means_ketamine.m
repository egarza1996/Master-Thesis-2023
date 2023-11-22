
function [allclusters, centroidcorr, dendmem, dunnsinitial]=meta_k_means_ketamine(eventsmat, distmetric);

if nargin==1;
    distmetric='correlation';
end

expr = eventsmat;
numden = size(expr,1);
numtimes = size(expr, 2);
counts = zeros(numden);
centroidmat = [];

dendevents=sum(eventsmat');
silentdends=find(dendevents<=1);
silentdendcnt=zeros(numden,1);
silentdendcnt(silentdends)=1;

if length(silentdends)>=numden-4;
    allclusters=[];
    dendmem=[];
    centroidcorr=[];
    dunnsinitial=[];
    
else
    % expr=zeros(numden-silentdends,numtimes); % original line
    expr=zeros(numden-length(silentdends),numtimes); % changed line
    cnt=1;
    for i=1:numden;
        
        if silentdendcnt(i)~=1;
            expr(cnt,:)=eventsmat(i,:);
            cnt=cnt+1;
        end
    end
    numden = size(expr,1);
    numtimes = size(expr, 2);
    counts = zeros(numden);
    centroidmat = [];

    
    % k clusters
    for k = 3
       
        % repeat t times
        for t = 1:500
            
            try
                
                [IDX, C, SUMD, D] = kmeans(expr, k, 'Distance', distmetric);

                % update counts
                % returns a matrix showing the number of times each pair of
                % dendrites were clustered together
                for x = 1:numden
                    for y = 1:numden
                        if IDX(x)==IDX(y)
                            counts(x,y) = counts(x,y) + 1;
                        end
                    end
                end

            catch
                IDX=[];
            end

        end
   
        if ~isempty(IDX);
            
            

            thr = 0.6; 



            clusters = [];
            for xx = 1:numden
                for yy = 1:xx-1
                    if counts(xx, yy) >= thr*t
                        clusters = [clusters xx];
                        clusters = [clusters yy];
                    end
                end
            end

            % sort list of dendrites
            clusters = sort(clusters);

            % delete the repeated dendrites from the list
            check = clusters(1);

            list = check;

            for c = 2:length(clusters)
                if clusters(c) ~= check
                    list = [list clusters(c)];
                    check = clusters(c);
                end
            end

            % vector of average point to centroid distances for each cluster
            numvect = [];

            cnt=1;
            dendmem=[];
            dendmem=dendmem';
            % find the clusters and print out graded memberships for the dendrites
            for n = 1:length(list)
                item = list(n);
                if item ~= 0
                    members = [item];
                    for nn = 1:numden
                        if nn ~= item
                            if counts(nn, item) >= thr*t
                                members = [members nn];
                                list(find(list == nn)) = 0;
                            end
                        end
                    end
                    eval(['cluster' num2str(cnt) '=members;']);
                    cnt=cnt+1;
                    % find centroid for cluster
                    centroid = zeros(1, numtimes);
                    for m = 1:length(members)
                        centroid = centroid + expr(members(m), :);
                    end
                    centroid = centroid/length(members);

                    centroidmat = [centroidmat; centroid];

                    numerator = 0;
                    % calculate the validity numerator
                    for m = 1:length(members)
                        numerator = numerator + (sum((expr(members(m),:)-centroid).^2))^0.5;
                    end
                    numerator = numerator/length(members);
                    numvect = [numvect numerator];

                    %graded membership using correlations
                    grmem = zeros(numden, 1);

                    for g = 1:numden
                        grmem(g, 1) = corr2(centroid, expr(g, :));
                    end
                    dendmem=[dendmem grmem];
                    grmem = [(1:numden); grmem'];

                end
            end

            allclusters=cell(cnt-1,2);
            for i=1:cnt-1;
                eval(['tempcluster=cluster' num2str(i) ';']);
                allclusters{i,1}=tempcluster;
                allclusters{i,2}=centroidmat(i,:);
            end


            centroidcorr = corrcoef(centroidmat');



            combclusters=0;

            [newclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, combclusters);

            endclustering=0;

            numclusters=10;

            while endclustering==0;
                numclusters=size(allclusters,1);
                combineclusters=[];
                for i=1:numclusters;
                    for j=(i+1):numclusters;
                        [newclusters, dunns, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, [i j]);
                        if dunns>=dunnsinitial;
                            combineclusters=[i j];
                            dunnsinitial=dunns;
                        end
                    end
                end
                if ~isempty(combineclusters);
                    [allclusters, dunnsinitial, centroidcorr, dendmem]=ioclustervalidity2(allclusters, expr, [combineclusters(1) combineclusters(2)]);
                else
                    endclustering=1;
                end
                if numclusters==3;
                    endclustering=1;
                end
            end

            numclusters=size(allclusters,1);

            for zz=1:numclusters;
                clusters=allclusters{zz,1};
                for zzz=1:length(clusters);
                    clusters(zzz)=clusters(zzz)+sum(silentdendcnt(1:clusters(zzz)));
                end
                allclusters{zz,1}=clusters;
            end
            
        else
            allclusters=[];
            centroidcorr=[];
            dendmem=[]; 
            dunnsinitial=0;
        end
        
    end
    
end



function [newclusters, dunns, centroidcorr, dendmem]=ioclustervalidity2(allclusters, eventsmat, combclusters);

numclusters=size(allclusters,1);
numtimes=size(eventsmat,2);

if combclusters~=0;
    a=allclusters{combclusters(1),1};
    b=allclusters{combclusters(2),1};
    c=[a b];
    c=sort(c);
    sizec=length(c);
    cnt=1;
    for i=1:sizec-1;
        d=0;
        for j=(i+1):sizec;
            if c(i)==c(j);
                d=1;
            end
        end
        if d==0;
            e(cnt)=c(i);
            cnt=cnt+1;
        end
    end
    e(cnt)=c(sizec);        
        
    eval(['cluster' num2str(combclusters(1)) '=e;']);
end



cnt=1;
for i=1:numclusters;
    if combclusters~=0;
        if cnt==combclusters(1);
            cnt=cnt+1;
        end
        if i~=combclusters(1) && i~=combclusters(2);
            eval(['cluster' num2str(cnt) '=allclusters{' num2str(i) ',1};']);
            cnt=cnt+1; 
        end
    else
        eval(['cluster' num2str(cnt) '=allclusters{' num2str(i) ',1};']);
        cnt=cnt+1; 
    end          
end

if combclusters~=0;
    numclusters=numclusters-1;
end

centroidmat=[];
dendmem=zeros(size(eventsmat,1),numclusters);

for i=1:numclusters;
    eval(['aa=cluster' num2str(i) ';']);
    nummembers=length(aa);
    centroid=zeros(1,numtimes);
    for j=1:nummembers;
        centroid(1,1:numtimes)=eventsmat(aa(j),1:numtimes)+centroid(1,1:numtimes);
    end
    eval(['centroid' num2str(i) '=centroid/nummembers;']);    
    centroidmat=[centroidmat; centroid/nummembers];
    for ii=1:size(eventsmat,1);
        dendmem(ii,i)=corr2(eventsmat(ii,:),centroid/nummembers);
    end
        
end

centroidcorr = corrcoef(centroidmat');

maxdiam=0;

for p=1:numclusters;
    eval(['clmembers=cluster' num2str(p) ';']);
    numpnts=length(clmembers);
    if numpnts~=1;
        for r=1:numpnts;
            for s=1:numpnts;
                %pairdist=(sum((eventsmat(clmembers(r),:)-eventsmat(clmembers(s),:)).^2))^0.5;
                pairdist=1-corr2(eventsmat(clmembers(r),:),eventsmat(clmembers(s),:));
                maxdiam=max(maxdiam,pairdist);
            end
        end
    end
end

mind=10000;

for p=1:numclusters;
    for r=1:numclusters;
        if p~=r;
             eval(['cl1members=cluster' num2str(p) ';']);
             eval(['cl2members=cluster' num2str(r) ';']);
             numpnts1=size(cl1members,2);
             numpnts2=size(cl2members,2);
             for s=1:numpnts1;
                for tt=1:numpnts2;
                    %pairdist=(sum((eventsmat(cl1members(s),:)-eventsmat(cl2members(tt),:)).^2))^0.5;
                    pairdist=1-corr2(eventsmat(cl1members(s),:),eventsmat(cl2members(tt),:));
                    if pairdist~=0;  
                        mind=min(mind,pairdist);
                    end
                end
             end
        end
     end
end

dunns = mind/maxdiam;

newclusters=cell(numclusters,2);

for i=1:numclusters;
    eval(['newclusters{' num2str(i) ',1}=cluster' num2str(i) ';']);
    eval(['newclusters{' num2str(i) ',2}=centroid' num2str(i) ';']);
end






allclusters=newclusters;
numdends=size(eventsmat,1);
numtimes=size(eventsmat,2);







                            
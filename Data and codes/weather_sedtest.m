
load weather_sed
%weather_sed=xlsread('C:\Users\杨依舟\Desktop\weather_indices.xlsx','weather_indices');
simitems={'Age';'CIA';'CIW';'PIA';'WIP';};

for i=1:length(simitems)
    agecol=strcmp(simitems,'Age');
    if ~any(agecol); error('Error: missing age variable'); end;
    if sum(agecol)>1; error('Error: multiple age variables'); end;
    data=zeros(length(weather_sed(:,1)),length(simitems));
    uncertainty=zeros(size(data));
    data(:,:)=weather_sed(:,:);
    uncertainty(:,:)=0.02;
    test=isnan(uncertainty(:,i))&~isnan(data(:,i));
    uncertainty(test,i)=mean(uncertainty(~isnan(uncertainty(:,i)),i));
end  

%分阶段筛选异常值，年龄值需要自己定
%     Age=data(:,i);
%     N=find(Age==551);
%     if length(N)>1 ;
%         N=N(1);
%     end 
%     M=find(Age==2501);
%     if length(M)>1;
%         M=M(1);
%     end

    %%Phanerozoic era
    % data1=weather_sed(1:N,2);
    % data1(find(data1<=0))=NaN;
    % T1=nanmean(log10(data1))+2*nanstd(log10(data1));
    % T2=nanmean(log10(data1))-2*nanstd(log10(data1));
    % data1(find(log10(data1)>T1 | log10(data1)<T2))=NaN;
    % 
    % data2=weather_sed(N+1:M,2);
    % data2(find(data2<=0))=NaN;
    % T3=nanmean(log10(data2))+2*nanstd(log10(data2));
    % T4=nanmean(log10(data2))-2*nanstd(log10(data2));
    % data2(find(log10(data2)>T3 | log10(data2)<T4))=NaN;
    %  
    % data3=weather_sed(M+1:end,2);
    % data3(find(data3<=0))=NaN;
    % T5=nanmean(log10(data3))+2*nanstd(log10(data3));
    % T6=nanmean(log10(data3))-2*nanstd(log10(data3));
    % data3(find(log10(data3)>T5 | log10(data3)<T6))=NaN;
    % data_f(:,1)=[data1;data2;data3];
    % data(:,2)=data_f(:,1);


MinAbsAgeUncert=200; 
test=uncertainty(:,agecol).*data(:,agecol) < MinAbsAgeUncert | isnan(uncertainty(:,agecol));
uncertainty(test,agecol)=MinAbsAgeUncert./data(test,agecol);

test=~isnan(weather_sed(:,1)); 
data=data(test,:);
uncertainty=uncertainty(test,:);

% 输出图形观察筛选效果
%     figure;
%     scatter(weather_sed(:,1),weather_sed(:,2),6,'MarkerEdgeColor','b');
%     hold on
%     scatter(data(:,1),data(:,2),6,'MarkerEdgeColor','r');

tic;
fprintf('Calculating sample weights: ')
    k=invweightAge(weather_sed(:,1));
prob=1./((k.*median(5./k))+1);
toc

 samplerows=100000;

 tic;
 fprintf('Allocating output matrix: ')
    % Generate output data matrix
 mczircon.data=NaN(samplerows,size(data,2));
 toc

 tic;
 fprintf('Resampling: ')
    % Fill the output data matrix with resampled data
i=1;
while i<samplerows
        % select weighted sample of data
    ru=rand(length(prob),1);
    sdata=data(prob>ru,:);
    suncertainty=uncertainty(prob>ru,:);
        % Randomize data over uncertainty interval
    rn=randn(size(sdata));
    sdata=sdata+rn.*suncertainty./2.*sdata;
    if i+size(sdata,1)-1<=samplerows
        mczircon.data(i:i+size(sdata,1)-1,:)=sdata;
    else
        mczircon.data(i:end,:)=sdata(1:samplerows-i+1,:);
    end
    i=i+size(sdata,1);
end
    
toc
mczircon.elements=simitems;
mczircon=elementify(mczircon);
save mczircon mczircon

agemin=0;
agemax=4000;
nbins =81;
binoverlap = 3;

i=length(simitems)-1;
figure;
while i>0
    Elem=simitems{i+1};
    test=~isnan(mczircon.(Elem));
    [c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(weather_sed(:,1)),nbins,binoverlap);
%     scatter(weather_sed(:,1),weather_sed(:,i+1),6,'MarkerEdgeColor','b');
%     hold on 
%     plot(c,m)
    subplot(2,2,i)
    errorbar(c,m,2*e,'.','LineWidth',1','MarkerSize',8)
    hold on 
    xlabel('Age (Ma)');
    ylabel(Elem); 
    xlim([agemin, agemax]);
    set(gca,'xdir','normal');
    set(gca,'fontsize',12,'fontweight','bold','FontName','Times New Roman');
    i=i-1;
end

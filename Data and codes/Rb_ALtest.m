
load Rb_AL
simitems={'Age';'Rb_AL';};

agecol=strcmp(simitems,'Age');
if ~any(agecol); error('Error: missing age variable'); end;
if sum(agecol)>1; error('Error: multiple age variables'); end;

data=zeros(length(Rb_AL),length(simitems));
uncertainty=zeros(size(data));
data(:,:)=Rb_AL(:,:);
uncertainty(:,:)=0.02;%误差为空值怎么样
test=isnan(uncertainty(:,1:2))&~isnan(data(:,1:2));
uncertainty(test,2)=mean(uncertainty(~isnan(uncertainty(:,2)),2));

Age=data(:,1);
N=find(Age==551);
if length(N)>1 ;
    N=N(1);
end 
M=find(Age==2501);
if length(M)>1;
    M=M(1);
end

%%Phanerozoic era
% data1=Rb_AL(1:N,2);
% data1(find(data1<=0))=NaN;
% T1=nanmean(log10(data1))+2*nanstd(log10(data1));%去除数据中的NAN值计算平均值，标准偏差
% T2=nanmean(log10(data1))-2*nanstd(log10(data1));
% data1(find(log10(data1)>T1 | log10(data1)<T2))=NaN;
% 
% data2=Rb_AL(N+1:M,2);
% data2(find(data2<=0))=NaN;
% T3=nanmean(log10(data2))+2*nanstd(log10(data2));
% T4=nanmean(log10(data2))-2*nanstd(log10(data2));
% data2(find(log10(data2)>T3 | log10(data2)<T4))=NaN;
%  
% data3=Rb_AL(M+1:end,2);
% data3(find(data3<=0))=NaN;
% T5=nanmean(log10(data3))+2*nanstd(log10(data3));
% T6=nanmean(log10(data3))-2*nanstd(log10(data3));
% data3(find(log10(data3)>T5 | log10(data3)<T6))=NaN;
% data_f(:,1)=[data1;data2;data3];
% data(:,2)=data_f(:,1);


MinAbsAgeUncert=0.2; 
test=uncertainty(:,agecol).*data(:,agecol) < MinAbsAgeUncert | isnan(uncertainty(:,agecol));
uncertainty(test,agecol)=MinAbsAgeUncert./data(test,agecol);

test=~isnan(Rb_AL(:,1)); % No point resampling by age if there isn't age data
data=data(test,:);
uncertainty=uncertainty(test,:);

figure;
scatter(Rb_AL(:,1),Rb_AL(:,2),6,'MarkerEdgeColor','b');
hold on
scatter(data(:,1),data(:,2),6,'MarkerEdgeColor','r');

tic;
fprintf('Calculating sample weights: ')
    k=invweightAge(Rb_AL(:,1));
    %U.k=k;
% end
prob=1./((k.*median(5./k))+1);
toc

samplerows=1E3;

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
mczircon=elementify(mczircon)
save mczircon mczircon

if ~exist('Rb_AL','var')
    load Rb_AL
end
%Elem='eHf_initial';
Elem='Rb_AL';
agemin=0;
agemax=4;
nbins =81;
binoverlap = 3;

test=~isnan(mczircon.(Elem));
%[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(zircon.Age),nbins,binoverlap);
[c,m,e]=bin(mczircon.Age(test),mczircon.(Elem)(test),agemin,agemax,length(mczircon.Age)./length(Rb_AL(:,1)),nbins,binoverlap);

figure; 
% scatter(Rb_AL(:,1),Rb_AL(:,2),6,'MarkerEdgeColor','b');
% hold on 
% plot(c,m,'LineWidth',1')
errorbar(c,m,2*e,'.','LineWidth',1','MarkerSize',8)
xlabel('Age (Ma)');
ylabel(Elem); 
xlim([0,4]);
%set(gca,'xdir','reverse ');
set(gca,'fontsize',12,'fontweight','bold','FontName','Times New Roman')
ylim([0 0.004]);
% %set(gca,'ytick',[-10 -5 0 5 10 15 20 25 30])
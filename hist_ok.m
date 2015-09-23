function hist_okt(X,n)

X=X(:);
X(isnan(X))=[];
X(isinf(X))=[];

%h=figure;
hist(X,n);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
set(gca,'FontSize',11);
%set(gca,'FontWeight','bold');
h1 = findobj(gca,'Type','patch');
set(h1,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1);

%M=nanmean(X(:));
%S=nanstd(X(:));
M=mean(X(:));
S=std(X(:));
RSD=S./M;
a=get(gca,'XLim');
x=max(a)-(max(a)-min(a))/4;
b=get(gca,'YLim');
y=max(b)-(max(b)-min(b))/10;

text(x,y,['\bf mean = ',num2str(round(M*100,3)/100),...
    '\newline','\bf std = ',num2str(round(S*100,3)/100),...
'\newline','\bf rsd = ',num2str(round(RSD*100,3)/100)]);
function Community=get_central_communities(weitrix)



degree=sum(weitrix,2);

[m,n]=size(weitrix);
average_degree=round(sum(degree)/m);
IndexMatrix=[1:m]';
AA=[IndexMatrix degree];
[B C]=sort(AA(:,2),'descend');
AA_descend=AA(C,:);    % sort the nodes according to the degree in descending order


l=length(find(AA_descend(:,2)<2));

AA_descend(end-l+1:end,:)=[]
;

AA_descend1=AA_descend;
for ii=1:m-l
    if ii>size(AA_descend1,1)
        break;
    end
    
    kk=AA_descend1(ii,1);     % node corresponding to the degree
    neighbor=IndexMatrix(logical(weitrix(:,kk)),:);   
    hola1=ismember(AA_descend1(:,1),neighbor);
    AA_descend1(hola1,:)=[];
end

BB=AA_descend1(:,1);
Community=num2cell(BB)';

end
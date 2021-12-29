function [eq]=Calculate_EQ(community,Adj_mat,degree,edges_num,vertex_num)
%    Calculate the value of EQ
%community    Is a cell structure

% global vertex_num degree edges_num EP;
% %clu_num = max(clu_assignment);
% degree=sum(Adj_mat,1);
% edges_num=sum(degree,2)/2;
% vertex_num=size(Adj_mat,1);

% for i = 1:clu_num
%     s_index = find(clu_assignment == i);
%     A=[];
%     for ii=1:length(s_index)
%         A=[A EP(s_index(ii),:)];
%     end
%     s_index1=unique(A);
%     community{i}=s_index1;
%     end

%Before calculating the EQ, first obtain a table, which represents 2*n, where n is the number of nodes that appear in the community (not including duplicates)
temp=unique([community{:,:}]);
% N=length(community);
% for i=1:N
%     a=community{i};
%     temp=[temp a];
%     temp=unique(temp);
% end;
temp2=[community{:,:}];
n=length(temp);

Table=zeros(2,n);

for j=1:n
    node=temp(j);
    num=sum(ismember(temp2,node));
    Table(1,j)=node;
    Table(2,j)=num;
end


N=length(community);

eq=0;
temp111=0;
temp222=0;
for i=1:N
    a=community{i};
    v_num=length(a);%   Is the number of nodes in each community
%     temp=0;
    temp11=0;
    temp22=0;
    for j=1:v_num
        v1=a(j);
        
        %  Count the number of v1 nodes in all communities
        contain_v1_num=Table(2,v1);
%         contain_v1_num=Table(2,find(Table(1,:)==v1));
%         contain_v1_num=contain(community,v1);
        for k=j+1:v_num
            v2=a(k);        
            
           
          contain_v2_num=Table(2,v2);
%             contain_v2_num=Table(2,find(Table(1,:)==v2));
            
            ttl=contain_v1_num*contain_v2_num;
%             contain_v2_num=contain(community,v2);   
            temp1=Adj_mat(v1,v2)/ttl;
            temp21=-(degree(v1)*degree(v2)/(2*edges_num))/ttl;
%             temp3=(temp1+temp2)/();
            temp11=temp11+temp1;
            temp22=temp22+temp21;
           
        end
    end
    temp111=temp111+temp11;
    temp222=temp222+temp22;
    
end
Q1=temp111/(2*edges_num);
Q2=temp222/(2*edges_num);
eq=Q1+Q2;


function num=contain_num(community,v)   
%  Function function: output the number of communities containing node v in the community
num=0;
n=length(community);
for i=1:n
    a=community{i};
    if any(a==v)==1
        num=num+1;
    end;
end;
end
end


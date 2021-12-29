    function Hab=extendNMI(A,B)
%NMI Normalized mutual information
% http://en.wikipedia.org/wiki/Mutual_information
% http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
% Author: http://www.cnblogs.com/ziqiao/   [2011/12/15] 
% A=[1 2 3 4;2 4 5 6];% Each row of the matrix represents the node number contained in a community, and A is divided into 2 groups
% B=[1 2 3 4;3 4 5 6];%  B is divided into 2 groups, and B is the real division!
numA=size(A,1); 
numB=size(B,1);
total1=max(max(B));
total2=max(max(A));
total=max(total1,total2);%When the algorithm cannot cover all nodes, special processing! Otherwise, an error will be reported!
%  Use a matrix idAOccur to store the partition of A, that is, the number of rows represents the number of communities, and the number of columns represents the number of nodes; the same for B
idAOccur=zeros(numA,total); 
for i=1:numA
    temp=find(A(i,:)~=0);
    a=A(i,temp);
    idAOccur(i,a)=1;
end
idBOccur=zeros(numB,total);
for i=1:numB
    temp=find(B(i,:)~=0);
    b=B(i,temp);
    idBOccur(i,b)=1;
end
Com11 = idAOccur * idBOccur';% Record the number of nodes shared by the two communities in the two divisions, which is a numA*numB matrix
sizeA=sum(idAOccur')';%  Record the size of each division in the A division, which is a vector of numA*1
sizeB=sum(idBOccur');%  Is a 1*numB vector

Com10=repmat(sizeA,1,numB)-Com11;
Com01=repmat(sizeB,numA,1)-Com11;
Com00=ones(size(Com11))*total-Com11-Com10-Com01;
Px = sum(idAOccur') / total;
Py = sum(idBOccur') / total;
P11=Com11/total;
P10=Com10/total;
P01=Com01/total;
P00=Com00/total;
min_A=[];% Store the degree of dissimilarity between each community in A and the corresponding community in B
for k=1:numA
    min_k=1000;
    Hk_l=0;
    for l=1:numB
        Hkl=0;h11=0;h10=0;h01=0;h00=0;
        if P11(k,l)~=0
            h11=-P11(k,l)*log2(P11(k,l)+eps);
            Hkl=Hkl+h11;
        end
        if P10(k,l)~=0
            h10=-P10(k,l)*log2(P10(k,l)+eps);
            Hkl=Hkl+h10;
        end
        if P01(k,l)~=0
            h01=-P01(k,l)*log2(P01(k,l)+eps);
            Hkl=Hkl+h01;
        end
        if P00(k,l)~=0
            h00=-P00(k,l)*log2(P00(k,l)+eps);
            Hkl=Hkl+h00;
        end       
        Py1=P11(k,l)+P01(k,l);
        Py0=P10(k,l)+P00(k,l);
        Hl=-Py1*log2(Py1+eps)+(-Py0*log2(Py0+eps));
        Hk_l=Hkl-Hl;
        if Hk_l<=min_k &&(h00+h11>h01+h10)
            min_k=Hk_l;
            templ=l;
        end
    end
    if min_k==1000  %      Restrictions are not met
        min_k=1;
    else
        Px1=P11(k,templ)+P10(k,templ);
        Px0=P01(k,templ)+P00(k,templ);
        Hk=-Px1*log2(Px1+eps)+(-Px0*log2(Px0+eps));
        min_k=min_k/Hk;
    end
    min_A=[min_A min_k];
end
Hab=sum(min_A)/numA;
end


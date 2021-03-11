function [ BS,Xbest,ft14 ] = HGCLFOA(fhd,sizepop,Dim,max_nfes,UB,LB,func_num)


fhd=str2func('cec17_func');
nfes=0;
max_nfes=Dim*10000;
maxgen=max_nfes/sizepop;
ft14=[];
g=1;
N_bi=0.35*sizepop;
N_p=sizepop-N_bi;
dmax=10^(-8);
A=[];
ci=0;
l=20;
SS={};



for i=1:sizepop  
    for j=1:Dim
        X(i,j)=LB+(UB-LB)*rand();
    end
    Smell(i,:)=feval(fhd,X(i,:)',func_num);
    nfes=nfes+1;
end

[X_s,X_i,NS, Smell_s, Smell_i]=Xguide(X,Smell,sizepop,N_bi,N_p);
BS=Smell_s(1,:);
Xbest=X_s(1,:);
gsmellbest(g,:)=Smell_s(1,:);
ft14= [ft14;ones(sizepop,1)*BS];

while nfes<max_nfes
    g=g+1;
    gnfes=0;
    
    
    weight(g)=1/(1+exp(-(0.5*(g-900)/20)));
    while weight(g)>1
        weight(g)=1;
    end
    %w=2^(1-((3*g)/maxgen));
    for i=2:N_bi
        D(i)=Euclidean(X_s(i,:),X_s(1,:));
        cost(i)=abs(Smell_s(i,:)-Smell_s(1,:));
    end
    
    d(2:N_bi)=D(2:N_bi)/max(D);
    c(2:N_bi)=cost(2:N_bi)/max(cost);
    P(1)=1;
    P(2:N_bi)=(d(2:N_bi)+exp(-c(2:N_bi)))/2;
    
   
    beta_f=1;
    beta_i=0.5;
    %            beta=((beta_f-beta_i)*(g/maxgen))+ beta_i;
    beta=1.5;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn*sigma;
    v=randn;
    step(g)=u/abs(v)^(1/beta);
    
    sa=2-2*(g/maxgen);
    
    count=0;
    
    for i=1:N_bi
        for j=1:NS(i)
            if rand()<=P(i)
                %                         X_i(count+1,:)=weight(g).*X_s(i,:)+(step*0.5).*(X_s(i,:)-X_i(count+1,:));
                X_i(count+1,:)=X_s(i,:)+(step(g)*0.5).*(X_s(i,:)-X_i(count+1,:));
                %
                
            else
                %                         if  rand()>sa/2
                q=randi([1,N_bi],1,1);
                X_mis=X_s(q,:)+rand*(X_s(q,:)-X_i(count+1,:));
                p=randi([1,N_bi],1,1);
                sita=(sa/2).*( X_i(count+1,:)-X_s(p,:));
                X_i(count+1,:)= X_mis+sita+randn(1,Dim);
                
                
            end
            
            X_i(count+1,:)=min(X_i(count+1,:),UB);
            X_i(count+1,:)=max(X_i(count+1,:),LB);
            
            Smell_i(count+1,:)=feval(fhd,X_i(count+1,:)',func_num);
            nfes=nfes+1;
            gnfes=gnfes+1;
            
            count=count+1;
        end
    end
    
    
    count=0;
    for i=1:N_bi
        for j=1:NS(i)
            if Smell_i(count+1,:)<Smell_s(i,:)
                tempX_s=X_s(i,:);
                X_s(i,:)=X_i(count+1,:);
                X_i(count+1,:)=tempX_s;
                
                tempSX_s=Smell_s(i,:);
                Smell_s(i,:)=Smell_i(count+1,:);
                Smell_i(count+1,:)=tempSX_s;
                
                if Smell_s(i,:)<BS
                    BS=Smell_i(count+1,:);
                    Xbest=X_s(i,:);
                end
                
            end
            
            
            count=count+1;
        end
    end
    
    
    [Smell_s,index_xs]=sort(Smell_s);
    for i=1:N_bi
        tempxs(i,:)=X_s(index_xs(i,:),:);
    end
    X_s=tempxs;
    
    
    
   
    
    avgSxs(g,:)=sum(Smell_s,1)/N_bi;
    if avgSxs(g,:)==avgSxs(g-1,:)
        
        for i=2:N_bi
            t=randi([1,i],1,1);
            while t == i
                t=randi([1,i],1,1);
            end
            
            %                             sita=abs((log10(g)/g).*(X_s(i,:)-X_s(1,:)));
            sita2=abs((sa/2).*( X_s(i,:)-X_s(t,:)));
            TX_s(i,:)=mvnrnd(X_s(t,:),sita2)+rand*X_s(t,:)-rand*X_s(i,:);
            %                             TX_s(i,:)=mvnrnd(X_s(1,:),sita)+rand*X_s(1,:)-rand*X_s(i,:);
            TSmell_s(i,:)=feval(fhd,TX_s(i,:)',func_num);
            nfes=nfes+1;
            gnfes=gnfes+1;
            
            if TSmell_s(i,:)<Smell_s(i,:)
                X_s(i,:)=TX_s(i,:);
                Smell_s(i,:)=TSmell_s(i,:);
                
                if Smell_s(i,:)<BS
                    BS=Smell_s(i,:);
                    Xbest=X_s(i,:);
                end
            end
        end
        
    else
        
        
        
        
        for i=1:N_bi
            %                       if i==1
            %                           w(i)=1/N_bi;
            %                       else
            sw=sum(log(N_bi+1)-log(1:i));
            w(i,:)=log(N_bi+1)/(sw);
            %                       end
        end
        %                   L0_mean=sum(w(1:N_bi,:).*X_s(1:N_bi,:));
        L0_mean1=mean(X_s(1:N_bi,:));
        for a=1:N_bi
            ws=w(a,:).*X_s(a,:);
        end
        Sw=sum(ws);
        Xs_mean2=Sw/N_bi;
        C1=cov(X_s);
        H=[A;X_s];
        C2=Cv(H,Xs_mean2);
        if ci<=l
            A=[cell2mat(SS(:,:));X_s];
            ci=ci+1;
        else
            A=[cell2mat(SS(ci-l+1:ci-1,:));X_s];
        end
        
        SS(g,:)={X_s};
        %                  C2=covarianceMatrix(X_s,Xs_mean2);
        for i=1:N_bi
            L_ms=(Xs_mean2+X_s(i,:))/2;
            mu=zeros(1,Dim);
            TX_s(i,:)=L_ms+mvnrnd(mu,C2);
        end
        TSmell_s(i,:)=feval(fhd,TX_s(i,:)',func_num);
        nfes=nfes+1;
        gnfes=gnfes+1;
        
        if TSmell_s(i,:)<Smell_s(i,:)
            X_s(i,:)=TX_s(i,:);
            Smell_s(i,:)=TSmell_s(i,:);
            
            if Smell_s(i,:)<BS
                BS=Smell_s(i,:);
                Xbest=X_s(i,:);
            end
        end
        
    end
    
    
    
    
    
    
    
    [Smell_s,index_xs ]=sort(Smell_s);
    for i=1:N_bi
        tempxs(i,:)=X_s(index_xs(i,:),:);
    end
    X_s=tempxs;
    
    count=1;
    [Smell_i,indexRS]=sort(Smell_i);
    for i=1:N_p
        temp(i,:)=X_i(indexRS(i,:),:);
    end
    X_i=temp;
    
    cs=[Smell_s;Smell_i(1,:)];
    CN=cs-max(cs);    % CN < 0
    
    NS=round(abs(CN/sum(CN))*N_p);
    NS(end)=[];
    
    i=N_bi;
    while sum(NS)>N_p
        if NS(i)>1
            NS(i)=NS(i)-1;
        else
            i=i-1;
        end
    end
    
    i=1;
    while sum(NS)<N_p
        NS(i)=NS(i)+1;
    end
    
    if find(NS==0)
        index=find(NS==0);
        for i=1:size(index,1)
            while NS(index(i))==0
                NS(index(i))=NS(index(i))+1;
                NS(i)=NS(i)-1;
            end
        end
    end
    NS=sort(NS,'descend');
    
    
    %              dmax=dmax-dmax/maxgen;
    gsmellbest(g,:)=Smell_s(1,:);
    ft14= [ft14;ones(gnfes,1)*Smell_s(1,:)];
    gbestX(g,:)=X_s(1,:);
    G(g,:)=g;
end

% plot(G,gsmellbest)
%            plot(G,step)
end





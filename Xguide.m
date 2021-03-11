function [L_0,L_p,NS, Smell_L0, Smell_p] = Xdesort(X,Smell,sizepop,N_bi,N_p)
            [Smellsort indexsort]=sort(Smell);
            for i=1:sizepop
                X_s(i,:)=X(indexsort(i,:),:);
            end

            %--------------------------划分精英,为不同个体分配不同引导对象--------------------------%
            BS=Smellsort(1,:);
            X_best=X_s(1,:);

            L_0=X_s(1:N_bi,:);
            Smell_L0=Smellsort(1:N_bi,:);  %精英种群

            L_p=X_s(N_bi+1:sizepop,:);
            Smell_p=Smellsort(N_bi+1:sizepop,:);  %差解种群

           

            cs=[Smell_L0;Smell_p(1,:)];
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
            NB=NS(2:end);
end
        
clear variables;
clc;
close all;


for i=1:1:3
    s=strcat('C:\Users\vladv\Desktop\Svertka\cat\cat', num2str(i), '.bmp');
    
    I=imread(s);
    
    %    INew=( I(:,:,1)+I(:,:,2)+I(:,:,3) )/3;
    
    h1=[1, 0, -1; 2, 0, -2; 1, 0, -1;];
%     I1=imfilter(I,h1);
    
    h2=[1, 2, 1; 0, 0, 0; -1, -2, -1;];
%     I2=imfilter(I,h2);
    
%     I=(I1+I2)/2;

    
    subplot(2,4,i)
    %     I=I3;
    I3=(I(:,:,1)*0.333+I(:,:,2)*0.333+I(:,:,3)*0.333);
    I3=abs(double(I3));
    I4{i}=(I3);
    
    imshow(I4{i});
    
    
end


%%%%% Входная картинка
inpu=I4;

target{1}=[1];
target{2}=[1];
target{3}=[0];

%%%%% Промежуточные по горизонтали [строк, столбцов;];
hid=[35,35;
     20,20;];  

 
sloi=[size(inpu{1});hid;size(target{1})];

N=size(sloi,1);


%%%%% Ядрa
for i=1:1:N-1
    core{i}=rand( sloi(i,1)-sloi(i+1,1)+1, sloi(i,2)-sloi(i+1,2)+1 )/10;
end
core{end}


a=0.3;
kolvo_epoh=100;
for epoh=1:1:kolvo_epoh
    
epoh

    for i=1:1:length(inpu)
    
        sloy{1}=inpu{i};
        sloy=forwRaspr(sloy,core,N);
    
        for H=N:-1:2
            if H==N
                errorHid{H-1}=target{i}-sloy{H};
            else
                errorHid{H-1}=back(errorHid{H},core{H});
            end
            delta{H-1}=errorHid{H-1};%.*(1-sloy{H}).*sloy{H};
            delW{H-1}=a*forw(sloy{H-1},delta{H-1});    %%%%%% ???    3*3'*(5*5) =  ??? А так ли? Проверить с перцептроном.
        end
        
        %%% Редактирование весов
        for H=1:1:N-1
            core{H}=core{H}+delW{H};
        end
    
    end
    
end


% inpu{3}
% inpu{2}

sloy{1}=inpu{3};
sloy=forwRaspr(sloy,core,N);
sloy{end}

core{end};
%%%%% Обратное распространение
% out=sloy{end};
% P{N}=out;
% for i=N-1:-1:1
%     P{i}=back(P{i+1},core{i});
% end



function sloy=forwRaspr(sloy,core,N)

%%%%% Прямое распространение
for i=1:1:N-1
   
    sloy{i+1}=forw(sloy{i},core{i});
    sloy{i+1}=activ(sloy{i+1});

end

end

function C=forw(A,B) %%% Находим С через А,   Картинка A, Ядро B, Выходная картинка С

C=zeros(size(A,1)-size(B,1)+1, size(A,2)-size(B,2)+1);

for i=1:1:size(A,1)-size(B,1)+1
    for j=1:1:size(A,2)-size(B,2)+1
        
        C(i,j)=sum(sum(  A(i:1:i+size(B,1)-1, j:1:j+size(B,2)-1).*B  )) ;
        
    end
end


end

function A=back(C,B) %%% Находим А через С,   Картинка A, Ядро B, Выходная картинка С

A=zeros(size(B,1)+size(C,1)-1, size(B,2)+size(C,2)-1);
%%% Перевернуть B на 180
B=B(end:-1:1,end:-1:1);

Map=zeros( size(C,1)+2*(size(B,1)-1), size(C,2)+2*(size(B,2)-1) );
Map( size(B,1):1:end-size(B,1)+1, size(B,2):1:end-size(B,2)+1 ) = C ;

for i=1:1:size(B,1)+size(C,1)-1
    for j=1:1:size(B,2)+size(C,2)-1
        
         A(i,j)=sum(sum(  Map(i:1:i+size(B,1)-1, j:1:j+size(B,2)-1).*B  )) ;
         
    end
end

end

function sloy=activ(sloy)

sloy(sloy<3000)=0; % sloy=1./( 1 + exp(-sloy./128) );
sloy(sloy>=3000)=1;

end




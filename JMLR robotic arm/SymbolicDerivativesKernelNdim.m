% We use symbolic calculus to compute the first-order derivatives (DmaxOrder=1) of the Gaussian
% kernel over inputs x and y both of dimension dim_in. The computations are
% then stored in a cell array GaussianFuncDeriv_mat_dim_*dim_in*_Dmax_1. This allows to
% avoid errors when performing computations involving these derivatives
% when computing kernel Gram matrices (typically, forgetting which variable
% is derived w.r.t.). The code writes scripts in text, to be passed to the
% symbolic differentiation diff() and to the function eval() to then become
% matlab functions.

dim_in=2;
a = sym('a',[1 dim_in]);
b = sym('b',[1 dim_in]);
param = sym('param',[1 dim_in]);

strCmd_create='syms ';% list of variables
strCmd_Gaussfunction='h=exp(-(';% Gaussian kernel to derive
strCmd_scaledSQdistfunction='';% Gaussian kernel to derive
for i=1:dim_in
    strCmd_create=[strCmd_create,['x',num2str(i),' ','y',num2str(i),' ']];%,'param',num2str(i),' '
    strCmd_Gaussfunction=[strCmd_Gaussfunction,['(x',num2str(i),'-','y',num2str(i),')^2','/(2*param(',num2str(i),')^2)+']];
    strCmd_scaledSQdistfunction=[strCmd_scaledSQdistfunction,['(x',num2str(i),'-','y',num2str(i),')^2','/(param(',num2str(i),')^2)+']];
end
strCmd_Gaussfunction=[strCmd_Gaussfunction(1:end-1),'));'];

strCmd_Matern52function=['h=(1+sqrt(5)*sqrt(', strCmd_scaledSQdistfunction(1:end-1),...
    ')+ 5/3*', strCmd_scaledSQdistfunction(1:end-1),')*exp(-(sqrt(5)*sqrt(', strCmd_scaledSQdistfunction(1:end-1),')));'];
eval(strCmd_create)
eval(strCmd_Gaussfunction)
% eval(strCmd_Matern52function)
% strCmd=[strCmd,'param'];
DmaxOrder=1;
% h=kFunc(x,y,param);
% =0.5;
% h=exp(-((x1-y1)^2+(x2-y2)^2)/(2*param^2));

%CODE WRITTEN FOR DmaxOrder=1, NEEDS CHANGES OTHERWISE
n_possibleDeriv=2*dim_in;
strCmd_DerivCellArr=repmat('DmaxOrder+1,',1,n_possibleDeriv);
strCmd_DerivCellArr=['kFuncDeriv_mat=cell(',strCmd_DerivCellArr(1:end-1),');'];
eval(strCmd_DerivCellArr);

bool_mat = false(n_possibleDeriv, 2^n_possibleDeriv);
count=1;
for i=1:n_possibleDeriv
    tempPos = nchoosek(1:n_possibleDeriv,i);
    for j=1:size(tempPos,1)
        bool_mat(tempPos(j,:),count)=true;
        count=count+1;
    end
end

for v=bool_mat
    strCmd_CompDeriv='diff(h,';
    for i=1:dim_in
        strCmd_CompDeriv=[strCmd_CompDeriv, repmat(['x',num2str(i),','], 1, v(i))];
    end
    for i=1:dim_in
        strCmd_CompDeriv=[strCmd_CompDeriv, repmat(['y',num2str(i),','], 1, v(i+dim_in))];
    end 
    strCmd_CompDeriv=[strCmd_CompDeriv(1:end-1),')'];
    strCmdX='[';
    strCmdY='[';
    for i=1:dim_in
    strCmdX=[strCmdX,['x',num2str(i),',']];
    strCmdY=[strCmdY,['y',num2str(i),',']];
    end
    strCmdX=[strCmdX(1:end-1),']'];
    strCmdY=[strCmdY(1:end-1),']'];
    v_str=reshape([num2str(v+1),repmat(',',n_possibleDeriv,1)]',[],1);
    v_str = ['kFuncDeriv_mat{',v_str(1:end-1)','}='];
	temp_eval=matlabFunction(simplify(subs(subs(eval(strCmd_CompDeriv),...   
        eval(strCmdX),a),eval(strCmdY),b)), 'vars', {a,b,param});%,'Steps',100,'IgnoreAnalyticConstraints',true
    eval([v_str,func2str(temp_eval),';']);
end
eval([v_str,func2str(matlabFunction(subs(subs(h,...
        eval(strCmdX),a),eval(strCmdY),b), 'vars', {a,b,param})),';']);

save(['GaussianFuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_',num2str(DmaxOrder)],'kFuncDeriv_mat');

% save(['Matern52FuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_',num2str(DmaxOrder)],'kFuncDeriv_mat');


%% VARIOUS FUNCTIONS TO CHECK AND VISUALIZE THE COMPUTATIONS ABOVE   
% figure    
% x_grid=linspace(-5,5,100);
% y_grid=linspace(-5,5,100);
% s=1;%linspace(-5,5,100)
% hold on
% plot(y_grid,kFuncDeriv_mat{1,1}(0,y_grid,1))
% plot(y_grid,kFuncDeriv_mat{1,2}(0,y_grid,1))
% plot(y_grid,kFuncDeriv_mat{2,1}(0,y_grid,1))
% plot(y_grid,kFuncDeriv_mat{2,2}(0,y_grid,1))

% matlabFunction(subs(subs(h,...
%         [x1,x2],a),[y1,y2],b), 'vars', {a,b,param});
    
    
% % kFuncDeriv_mat=cell(DmaxOrder+1,DmaxOrder+1,DmaxOrder+1,DmaxOrder+1);
% for r1=0:DmaxOrder
% for r2=0:DmaxOrder
% for s1=0:DmaxOrder
% for s2=0:DmaxOrder
%     stringCommand=['diff(h,', repmat('x1,', 1, r1), repmat('x2,', 1, r2),...
%         repmat('y1,', 1, s1),repmat('y2,', 1, s2)];
%     stringCommand=stringCommand(1:end-1);
%     stringCommand=[stringCommand,')'];
% %     stringCommand=[ num2str(j) '(t)']),r1,x2,,y1,s1,y2,s2
% 	kFuncDeriv_mat{r1+1,r2+1,s1+1,s2+1}=matlabFunction(subs(subs(eval(stringCommand),...
%         [x1,x2],a),[y1,y2],b), 'vars', {a,b,param});
% %     kFuncDeriv_mat_short{r1+1,r2+1,s1+1,s2+1}=@(x,y,param) tempFunc(param,x(1),x(2),y(1),y(2));
% end    
% end    
% end    
% end
% % kFuncDeriv_mat{1,1,1,1}=matlabFunction(subs(subs(h,...
% %         [x1,x2],a),[y1,y2],b), 'vars', {a,b,param});
%     
% save(['GaussianFuncDeriv_mat_dim_',num2str(dim_in),'_Dmax_',num2str(dim_in)],...
%     'kFuncDeriv_mat');
%%
% blob2= matlabFunction(subs(subs(h,...
%         [x1,x2],a),[y1,y2],b), 'vars', {a,b,param})
% blob2([0,0;1 0],[0,0; 1 0],1)
    %
% kFuncDeriv_mat{1,1,1,1}=@(x,y,param) hFunc(param,x(1),x(2),y(1),y(2));%matlabFunction(h);
% matlabFunction(kFuncDeriv_mat)
%%
% htest=kFuncDeriv_mat{1,1,1,2};
% hest2=@(param,x1,x2,y1,y2) hFunc(param,x1,x2,y1,y2);
% htest
% hest2(0.1,0,1,0,0)
% outputArg1 = inputArg1;
% outputArg2 = inputArg2;
% end

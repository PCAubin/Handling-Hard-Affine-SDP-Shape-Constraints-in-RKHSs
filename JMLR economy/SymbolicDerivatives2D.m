% We use symbolic calculus to compute the first-order derivatives (DmaxOrder=0) of the Gaussian
% kernel over 2D inputs x and y. The computations are
% then stored in a cell array kFuncDeriv_mat. This allows to
% avoid errors when performing computations involving these derivatives
% when computing kernel Gram matrices (typically, forgetting which variable
% is derived w.r.t.). The code writes scripts in text, to be passed to the
% symbolic differentiation diff() and to the function eval() to then become
% matlab functions.


clear all
a = sym('a',[1 2]);
b = sym('b',[1 2]);
syms x1 x2 y1 y2 param
DmaxOrder=0;
% h=kFunc(x,y,param);
% =0.5;
h=exp(-((x1-y1)^2+(x2-y2)^2)/(2*param^2));

kFuncDeriv_mat=cell(DmaxOrder+1,DmaxOrder+1,DmaxOrder+1,DmaxOrder+1);
% kFuncDeriv_mat_short=cell(DmaxOrder+1,DmaxOrder+1,DmaxOrder+1,DmaxOrder+1);
for r1=0:DmaxOrder
for r2=0:DmaxOrder
for s1=0:DmaxOrder
for s2=0:DmaxOrder
    stringCommand=['diff(h,', repmat('x1,', 1, r1), repmat('x2,', 1, r2),...
        repmat('y1,', 1, s1),repmat('y2,', 1, s2)];
    stringCommand=stringCommand(1:end-1);
    stringCommand=[stringCommand,')'];
%     stringCommand=[ num2str(j) '(t)']),r1,x2,,y1,s1,y2,s2
	kFuncDeriv_mat{r1+1,r2+1,s1+1,s2+1}=matlabFunction(subs(subs(eval(stringCommand),...
        [x1,x2],a),[y1,y2],b), 'vars', {a,b,param});
%     kFuncDeriv_mat_short{r1+1,r2+1,s1+1,s2+1}=@(x,y,param) tempFunc(param,x(1),x(2),y(1),y(2));
end    
end    
end    
end
kFuncDeriv_mat{1,1,1,1}=matlabFunction(subs(subs(h,...
        [x1,x2],a),[y1,y2],b), 'vars', {a,b,param});
    
% save('GaussianFuncDeriv_mat','kFuncDeriv_mat');
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


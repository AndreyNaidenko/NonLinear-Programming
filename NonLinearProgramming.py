import sympy as sp


def get_s(n,sym):
    '''
    n:int-количество переменных
    sym:char - символ 
    
    Пример:
        
        get_s(3,'x')
        
        return: [x1,x2,x3]
    '''
    s=[]
    for i in range(n):
        s.append(sp.symbols(sym+str(i+1)))
    return s
def get_f(f):
    '''
    f :string -функция
    Пример:f="x1+x2"
    '''
    F=sp.sympify(f)
    return F

def get_conditions(c,m):
    '''
    c :[],string -условия
    m :int -количество условий 
    Пример:
        x1+x2 <5
        x1+5*x2>10
        с[0]="5-x1-x2"
        c[1]="10-x1-5*x2"
    '''
    conditions=[]
    if len(c)!=m :
        print("Количество условий должно равнятся m.")
    else :
        for i in range(m):
            conditions.append(sp.sympify(c[i]))
     
        return conditions
    
def get_L(F,conditions):
    '''
    F:sympy - исходная функиция
    conditions :[] sympy - условия
    '''
    L=F
    y=get_s(len(conditions),'y')
    for i in range(len(conditions)):
       L=L+y[i]*conditions[i]
    return L    

def diff_L(L,n,m):
    '''
    L:sympy - функция L(x,y)
    n:int - количество переменных
    m:int -количество условий
    '''
    dL=[]
    x=get_s(n,'x')
    y=get_s(m,'y')
    for i in range(n):
        dL.append(sp.diff(L,x[i]))
    for i in range(m):
        dL.append(sp.diff(L,y[i]))
    return dL,x,y

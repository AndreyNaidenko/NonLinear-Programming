import sympy as sp
import numpy as np



def cond_to_tex(cond,sgn, add_variable_y=False, add_variable_v_w=False, add_variable_z=False):
    TeX=""
    TeX=TeX+"$\\begin{cases} \r\n"
    number_variable = 2
    for c in range(len(cond)):
        s="="
        if sgn[c][0]=="<":
            s="\\leq"
        if sgn[c][0]==">":
            s="\\geq"
        TeX=TeX+"\t"+sp.latex(cond[c])+s+" 0 \\\\ \r\n"
    if not add_variable_y and not add_variable_v_w and not add_variable_z:
        TeX = TeX + "x_i \geq 0, \\ i = \overline{1," + str(number_variable) + "}\\end{cases}$ \r\n"
    if add_variable_y and not add_variable_v_w and not add_variable_z:
        TeX = TeX + "x_i, y_i \geq 0, \\ i = \overline{1," + str(number_variable) + "}\\end{cases}$ \r\n"
    if add_variable_y and add_variable_v_w and not add_variable_z:
        TeX = TeX + "x_i, y_i, v_i, w_i \geq 0, \\ i = \overline{1," + str(number_variable) + "}\\end{cases}$ \r\n"
    if add_variable_y and add_variable_v_w and add_variable_z:
        TeX = TeX + "x_i, y_i, v_i, w_i, z_i \geq 0, \\ i = \overline{1," + str(number_variable) + "}\\end{cases}$ \r\n"
    return TeX
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

def get_conditions(c,cb,m):
    '''
    ai1*xi1+....+ain*xin = bi
    
    c :[],string - ai1*xi1+....+ain*xin
    cb:[],float -  bi
    m :int -количество условий 
    Пример:
        x1+x2   =5
        x1+5*x2 =10
        с[0]="x1-x2"
        c[1]="x1-5*x2"
        cb[0]=5
        cb[1]=10
    '''
    conditions=[]
    if len(c)!=m :
        print("Количество условий должно равнятся m.")
    else :
        for i in range(m):
            conditions.append(sp.sympify(c[i])-cb[i])
     
        return conditions
    
def get_L(F,conditions,variables):
    '''
    F:sympy - исходная функиция
    conditions :[] sympy - условия
    '''
    L=F
    y=get_s(len(conditions),'y')
    for i in range(len(conditions)):
       variables.append(y[i])
       L=L+y[i]*(-1*conditions[i])
    return L,variables   

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
def get_b(conditions):
    b=[]
    for i in range(len(conditions)):
        a=sp.Poly(conditions[i]).coeffs()
        b.append(-1*a[len(a)-1])
     
    return b
def canonical(conditions,sgn,variables):
    for i in range(int(len(conditions)/2)):
        s=sp.symbols('v'+str(i+1))
        variables.append(s)
        if sgn[i]=="<":
            conditions[i]=conditions[i]+s           
        else :
            conditions[i]=conditions[i]-s
        sgn[i]="="
    k=0
    for i in range(int(len(conditions)/2),len(conditions)):
        s=sp.symbols('w'+str(k+1))
        variables.append(s)
        k=k+1
        if sgn[i]=="<":
            conditions[i]=conditions[i]+s           
        else :
            conditions[i]=conditions[i]-s
        sgn[i]="="
    return conditions,sgn,variables
def simplex(F,F1,conditions,b,n,m,x):
    '''
    F  
    F1
    
    сделано для того что бы различить M и M
    в F хранится M 
    в F1 хранится 10000
    '''
    TeX="\\begin{tabular}{"
    for i in range(n+2):
        TeX=TeX+"|c"
    TeX=TeX+"|}\r\n"
        
        
    '''
    
    vb хранит базисные переменные в виде символов
    Пока вводится в ручную
    '''
    vb=[x[n-2],x[n-1],x[n-3],x[n-4]]
    
    f_coeffs=[]
    f_coeffs1=[]
    STOP=False
    f_coeffs=[]
    
    for i in range(n):
        f_coeffs.append(F.coeff(x[i]))
        f_coeffs1.append(F1.coeff(x[i]))
    ST = sp.Matrix.zeros(m,n+1)
    for i in range(m):
        ST[i,n]=sp.Rational(b[i])
    
   
    for i in range(m):
        c=np.zeros((n))
        for k in range(n):
            c[k]=conditions[i].coeff(x[k])
        for j in range(n):
            ST[i,j]=sp.Rational(c[j]) 
    basis=sp.Matrix.zeros(m,1)
    '''
    Тут нужно выбрать базис. 
    я его задаю в ручную.
    
    basis хранит M как 'M'
    basis1 хранит M как 100000
    '''
    M=sp.symbols('M')
    basis[0]=-M
    basis[1]=-M
    basis1=sp.Matrix.zeros(m,1)
    M=100000
    basis1[0]=-M
    basis1[1]=-M
    
    TeX=""
    NNNN=0
    TeX = TeX + "\\end{tabular}\r\n \\vspace{2mm}"
    while STOP != True and NNNN<6:
            
            f_coeffs=[] # коэффициенты функции F
            for i in range(n):
                f_coeffs.append(F.coeff(x[i]))
            M=sp.symbols('M') # M  ,которя в искуственном базисе M > > 1
            TeX=TeX+"\r\n \\vspace{2mm} \\begin{tabular}{"
            for i in range(n+3):
                 TeX=TeX+"|c"
            TeX=TeX+"|}\r\n  \hline \r\n"   
             
            TeX=TeX+"\tБазис&\tС"
            for i in range(n):
                TeX=TeX+"\t&$"+sp.latex(f_coeffs[i])+"$"
            TeX=TeX+"\t&\\\\ \r\n \hline \r\n \t& "
            
            for i in range (n):
                TeX=TeX+"\t&$"+sp.latex(x[i]) +"$"
            TeX=TeX+"\t&$b_{i}$ \\\\ \r\n \hline \r\n"
            
            for i in range (m):
                a="\t $"+sp.latex(vb[i])+"$&$"+sp.latex(basis[i])+"$"
                TeX=TeX+a
                 
                for j in range(n+1):
                    TeX=TeX+"\t&$"+sp.latex(ST[i,j]) +"$"
                TeX=TeX+" \\\\ \r\n \hline \r\n"
            
            
             
            TeX=TeX+"\\end{tabular}\r\n \\vspace{2mm}"
            for i in range(n):
                TeX=TeX+"$\\Delta_{"+str(i+1)+"}="+sp.latex(basis.transpose())+"\\cdot"+sp.latex(ST.col(i))+"+";
                if "-" in sp.latex(-1*f_coeffs[i]):
                    TeX=TeX+"("+sp.latex(-1*f_coeffs[i])+")="+sp.latex((basis.transpose()*ST.col(i))[0]-f_coeffs[i])+"$\r\n\\newline"
                else:
                    TeX = TeX + sp.latex(-1 * f_coeffs[i]) + "=" + sp.latex((basis.transpose() * ST.col(i))[0] - f_coeffs[i]) + "$\r\n\\newline \\vspace{2mm}"
            
            deltas=[] #массив хранящий дельты
            for i in range(n):
                deltas.append((basis1.transpose()*ST.col(i))[0]-f_coeffs1[i]) #вычисление дельт
           
             
            '''
            Цикл проходит по всем дельтам. Если дельта >= 0 , то записываем в нее 0 для того чтобы сохранить 
            размер массива
            '''
            for i in range(len(deltas)):
                if deltas[i]>=0:
                    deltas[i]=0
                deltas[i]=sp.Abs(deltas[i]) 
            k=deltas.index(max(deltas))# Находим номер максимальной дельты 
           
             
            if deltas[k]==0: # Если дельта равна нулю , то план оптимальный , т.к. в прошлом цикле все дельты больше нуля приравняли к 0
                STOP=True
                break
            TeX=TeX+"Вводим в базис\t$"+sp.latex(x[k])+"$\r\n"
            TeX=TeX+"\r\n \\vspace{2mm} Для определения вектора, подлежащего исключению из базиса, находим $min("
            mi=[]
            for i in range(m):
                if ST[i,k]==0:    
                    TeX=TeX+"\\infty,"
                    mi.append(1000000)
                else:                 
                    if ST[i,k]>0:
                        if ST[i,n]/ST[i,k]>=0: 
                            mi.append(ST[i,n]/ST[i,k])
                    else:
                        mi.append(1000000)
                    TeX=TeX+sp.latex(ST[i,n]/ST[i,k])+","
            r=mi.index(min(mi))
        
            
            TeX=TeX+")="+sp.latex(min(mi)) +"$\r\n \\vspace{2mm}"
            TeX=TeX+"Исключаем \t$"+sp.latex(vb[r])+"$\r\n"
            if(min(mi)>1000):
                break
            NST = sp.Matrix.zeros(m,n+1) 
            NST[r,n]=ST[r,n]/ST[r,k]
            for j in range(n):
                NST[r,j]=ST[r,j]/ST[r,k]
            
            for i in range(m):
                if i!=r:
                    NST[i,n]=ST[i,n]-ST[i,k]*(ST[r,n]/ST[r,k])
                for j in range(n):
                    if i!=r:
                        NST[i,j]=ST[i,j]-ST[i,k]*(ST[r,j]/ST[r,k])
            
                       
                       
                     
                        
                        
            ST=NST            
            basis1[r]=f_coeffs1[k]
            basis[r]=f_coeffs[k]
            vb[r]=x[k]
            
    TeX=TeX+"Все $\Delta_{i}$ не отрицательны , план является оптимальным \r\n"
         
    return TeX,vb,ST.col(n)
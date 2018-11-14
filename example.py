from NonLinearProgramming import *


'''
f="-2*x1-5*x2-1000*x6"

c=[]
c.append("2*x1-3*x2-x3-x6")
c.append("x1+2*x2+x4")
c.append("x1+x2+x5")
F=get_f(f)
cond=get_conditions(c,3)



simplex(F,cond,[180,100,95],6,3)
'''

f="2*x1+4*x2-x1**2-2*x2**2"
sgn=["<","<"]
c=["x1+2*x2","2*x1-x2"]
x=get_s(2,'x')
VAR=[x[0],x[1]]
cb=[8,12]
F=get_f(f)
cond=get_conditions(c,cb,2)
L,VAR=get_L(F,cond,VAR)
dL,x,y=diff_L(L,2,2)

TeX='\documentclass{article} \r\n'
packages=['[T2A]{fontenc}','[utf8]{inputenc}','[english,russian]{babel}','{amsmath}']
up="/usepackage"
for p in packages:
    TeX=TeX+"\\usepackage"+p+'\r\n'

TeX=TeX+"\\begin{document} \r\n"
TeX=TeX+"\\begin{center} \r\n"  
TeX=TeX+"$$F="+sp.latex(F) +"\\to max $$ \r\n"
TeX=TeX+cond_to_tex(cond,sgn);
TeX = TeX + "\\end{center} \r\n"
TeX=TeX+"\\newline \\vspace{5mm} Составляем функцию Лагранжа:"
TeX=TeX+"$$ L="+sp.latex(L)+"$$ \r\n"
TeX=TeX+"\\newline Запишем необходимые и достаточные условия существования седловой точки построенной функции:\r\n \\vspace{1mm}"
K=0

TeX = TeX + "$\\begin{cases}"
for i in range(len(x)):
    TeX=TeX+" \\vspace{1mm} \\dfrac{\\partial L}{\\partial "+sp.latex(x[i])+"}="+sp.latex(dL[K])+" \\leq 0\\\\  \\vspace{1mm} \r\n"
    K=K+1
for i in range(len(y)):
    TeX=TeX+"\\dfrac{\\partial L}{\\partial "+sp.latex(y[i])+"}="+sp.latex(dL[K])+" \\geq 0\\\\"
    K=K+1
TeX = TeX + "\\end{cases}$ \r\n \\vspace{1mm}"

K = 0

TeX = TeX + "$\\begin{cases}"
for i in range(len(x)):
    TeX=TeX+" \\vspace{1mm} " + sp.latex(x[i])+ "\\dfrac{\\partial L}{\\partial "+sp.latex(x[i])+"}="+sp.latex(x[i]*dL[K])+ "= 0\\\\  \\vspace{1mm} \r\n"
    K=K+1
for i in range(len(y)):
    TeX=TeX+" " + sp.latex(y[i]) +"\\dfrac{\\partial L}{\\partial "+sp.latex(y[i])+"}=" +sp.latex(y[i]*dL[K])+"= 0\\\\  \\vspace{1mm} \r\n"
    K=K+1
TeX = TeX + "\\end{cases}$\r\n"

sgn.append(">")
sgn.append(">")
new_cond=dL
#TeX=TeX+cond_to_tex(new_cond,sgn);
TeX=TeX+"\\newline \\vspace{2mm} Систему неравенств перепишем в следующим образом: \r\n \\vspace{2mm} "
for i in range(len(new_cond)):
    new_cond[i]=new_cond[i]*-1
    if sgn[i]==">":
        sgn[i]="<"
    else:
        sgn[i]=">"
TeX=TeX+cond_to_tex(new_cond,sgn, True)

TeX=TeX+"\\newline \\vspace{2mm} Вводя теперь дополнительные неотрицательные переменные $v_1, v_2, w_1$ и $w_2$, обратим неравенства в равенства: \r\n \\vspace{2mm} "

new_cond,sgn,VAR=canonical(new_cond,sgn,VAR)
TeX=TeX+cond_to_tex(new_cond,sgn, True, True);
b=get_b(new_cond)
z=get_s(2,'z')
m=sp.symbols('M')
VAR.append(z[0])
VAR.append(z[1])
G1=-1000*z[0]-1000*z[1]
G=-m*z[0]-m*z[1]
new_cond[0]=new_cond[0]+z[0]
new_cond[1]=new_cond[1]+z[1]
TeX = TeX + "\\vspace{2mm} Учитывая равенства, можно записать: $v_1 x_1 = 0, v_2 x_2 = 0, w_1 y_1 = 0, w_2 y_2 = 0.\r\n"
TeX=TeX+"\\newline \\vspace{2mm} Для нахождения системы линейных уравнений воспользуемся методом искусственного базиса. " \
        "В первое и второе уравнение системы соответственно добавим дополнительную неотрицательную переменную $z_1$ и $z_2$  " \
        "и рассмотрим задачу линейного программирования, сотостоящую в определении максимального значения функции \r\n \\vspace{2mm} "
TeX=TeX+"$$\overline{F}=-Mz_1-Mz_2$$ при условиях \r\n \\vspace{2mm}"
TeX=TeX+ "\\vspace{1mm} \r\n " + cond_to_tex(new_cond,sgn, True, True, True);
TeX=TeX+"\\newline \\vspace{8mm} Симплекс-метод: \r\n "
simpl,vb,bi=simplex(G,G1,new_cond,b,len(VAR),len(new_cond),VAR)
TeX=TeX+simpl

for i in range(len(vb)):
    if vb[i]==VAR[0]:
        a=bi[i]
    if vb[i]==VAR[1]:
        b=bi[i]
TeX=TeX+ "$F(X)="+sp.latex(F.evalf(subs={VAR[0]:a, VAR[1]:b})) + "$ \r\n"

TeX=TeX+"\\end{center}  \r\n" 
TeX=TeX+"\\end{document} \r\n"

my_file = open("some.txt", "w")
my_file.write(TeX)
my_file.close()

'''
c=['2*y2+7*y3-35*y1+y4','9*y2+4*y3-45*y1+y5','y2+y3+y6']

cb=[0,0,1]

f="5*y2+8*y3"

F=get_f(f)

F1=F
y=sp.symbols('y')
M=sp.symbols('M')

F=F-M*y

M=10000

F1=F1-M*y
v=get_s(6,'y')
c=get_conditions(c,cb,3)
for i in range(len(c)):
    c[i]=c[i]*-1
t=simplex(F,F1,c,cb,6,3,v)
print(t)
'''
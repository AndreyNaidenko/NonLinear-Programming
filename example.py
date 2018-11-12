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
TeX=TeX+"\\newline Составляем функцию L \r\n"
TeX=TeX+"$$ L(x,y)="+sp.latex(L)+"$$ \r\n"
TeX=TeX+"\\newline Находим частные производные от функции L \r\n"
K=0
for i in range(len(x)):
    TeX=TeX+"$$\\frac{\\partial L(x,y)}{\\partial "+sp.latex(x[i])+"}="+sp.latex(dL[K])+"$$\r\n"
    K=K+1
for i in range(len(y)):
    TeX=TeX+"$$\\frac{\\partial L(x,y)}{\\partial "+sp.latex(y[i])+"}="+sp.latex(dL[K])+"$$\r\n"
    K=K+1
    
sgn.append(">")
sgn.append(">")
new_cond=dL
TeX=TeX+cond_to_tex(new_cond,sgn);
for i in range(len(new_cond)):
    new_cond[i]=new_cond[i]*-1
    if sgn[i]==">":
        sgn[i]="<"
    else:
        sgn[i]=">"
TeX=TeX+cond_to_tex(new_cond,sgn)

new_cond,sgn,VAR=canonical(new_cond,sgn,VAR)
TeX=TeX+cond_to_tex(new_cond,sgn);
b=get_b(new_cond)
z=get_s(2,'z')
m=sp.symbols('M')
VAR.append(z[0])
VAR.append(z[1])
F1=F-1000*z[0]-1000*z[1]
F=F-m*z[0]-m*z[1]
new_cond[0]=new_cond[0]+z[0]
new_cond[1]=new_cond[1]+z[1]
TeX=TeX+cond_to_tex(new_cond,sgn);
TeX=TeX+simplex(F,F1,new_cond,b,len(VAR),len(new_cond),VAR)
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
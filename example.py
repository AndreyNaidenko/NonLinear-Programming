from NonLinearProgramming import *




f="2*x1+4*x2-x1**2-2*x2**2"

c=[]
c.append("8-x1-x2")
c.append("12-2*x1+x2")
F=get_f(f)
L=get_L(F,cond)
dL,x,y=diff_L(L,2,2)

TeX='\documentclass{article} \r\n'
packages=['[T2A]{fontenc}','[utf8]{inputenc}','[english,russian]{babel}','{amsmath}']
up="/usepackage"
for p in packages:
    TeX=TeX+"\\usepackage"+p+'\r\n'
TeX=TeX+"\\begin{document} \r\n"
TeX=TeX+"\\begin{center} \r\n"  
TeX=TeX+"$$F="+sp.latex(F) +"\\to max $$ \r\n"
TeX=TeX+"$\\begin{cases} \r\n"
cond=get_conditions(c,2)
for c in cond:
    TeX=TeX+"\t"+sp.latex(c)+"\\leq 0 \\\\ \r\n" 
TeX=TeX+"\\end{cases}$ \r\n" 
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
TeX=TeX+"\\end{center}  \r\n" 
TeX=TeX+"\\end{document} \r\n"
 
print(TeX)
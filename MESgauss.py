import math

def fun1(x):
    return 5*x**2+3*x+6
def fun2(x,y):
    return 5*(x**2)*(y**2)+3*x*y+6

class Gauss:
    a=1
    b=1
    def __init__(self,n):
        self.n=n
        if self.n==0:
            self.x=[0]
            self.w=[2]
        elif self.n==1:
            self.x=[-1/math.sqrt(3),1/math.sqrt(3)]
            self.w=[1,1]
        elif self.n==2:
            self.x=[-1*(math.sqrt(3/5)),0,math.sqrt(3/5)]
            self.w=[5/9,8/9,5/9]
        elif self.n==3:
            a1 = math.sqrt(3 / 7 - 2 / 7 * math.sqrt(6 / 5))
            a2 = math.sqrt(3 / 7 + 2 / 7 * math.sqrt(6 / 5))
            w1 = (18 + math.sqrt(30)) / 36
            w2 = (18 - math.sqrt(30)) / 36
            self.x = [-a2, -a1, a1, a2]
            self.w = [w2, w1, w1, w2]
        else:
            raise ValueError("Obslugiwane tylko n <0,3>")

    def mGaussa1d(self,funkcja):
        return sum(self.w[i] * funkcja(self.x[i]) for i in range(len(self.x)))

    def mGaussa2d(self,funkcja):
        wynik=0
        for i in range(self.n+1):
            for j in range(self.n+1):
                wynik+=self.w[i]*self.w[j]*funkcja(self.x[i],self.x[j])
        return wynik

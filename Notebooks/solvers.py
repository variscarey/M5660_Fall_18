import numpy as np

class Solver:
    
    def __init__(self,function,start=0,tolerance=1E-8,maxiter=20):
        self.fn=function
        self.guess=start
        self.tolerance=tolerance
        self.maxiter=maxiter
        self.iteration=0  #iteration counter
        self.root=None
    
    def state(self):
        print('Iterate=',self.guess)
        print('Iteration=',self.iteration)
        if self.root!=None:
            print('Root=',self.root)

class BracketSolver(Solver):
    def __init__(self,function,guess=0,tolerance=1E-8,maxiter=20,left=-1,right=1,verbose=False,switch=False,history=False):
            Solver.__init__(self,function,guess,tolerance,maxiter)
            self.left=left
            self.right=right
            self.lval=self.fn(left) #evaluate f at left bracket
            self.rval=self.fn(right) #evaluate f at right bracket
            self.switch=switch
            if self.lval*self.rval > 0: #error handling for invalide bracket
                print('Warning: No Root guarantee',self.fn(left),self.fn(right))
            if self.left >= self.right: #error handling for right endpoint before left endpoint
                print('Warning, Invalid Interval')
            self.verbose=verbose

    def state(self):
            print('Bracket=[',self.left,',',self.right,']')
            print('Values=(',self.lval,',',self.rval,')')
            Solver.state(self)
            
    def solve(self,type='bisection'):
        import matplotlib.pyplot  #for plotting
        from math import isnan #for error handling
        self.iteration+=1         
        old=None                  #for improved versions of regula falsi
        while self.iteration < self.maxiter:
            if self.verbose: #plot brackets, print system state
                self.state()
                matplotlib.pyplot.plot([self.left,self.right],[self.lval,self.rval],linestyle='dashed',marker='o')
            if type=='bisection': #perform bisection
                self.guess=.5*(self.left+self.right) #bisect bracket
                test=self.fn(self.guess)
                if isnan(test):
                    print('Current guess is Nan')
                    return
                if abs(test)<self.tolerance:  #we found a root
                    self.root=self.guess
                    return 
                else:   #iterate
                    #compute approximate fpi newton factor if self.switch=true.
                    if self.switch:
                        fpideriv=4*test*(2*test-self.rval-self.lval)/(self.rval-self.lval)**2 #f(x)f''(x)/f'(x)^2
                        print('Switch to fpi coeff',abs(fpideriv))
                        if abs(fpideriv) < .25:
                            return                    
                    if self.lval*test < 0:
                        self.rval=test
                        self.right=self.guess
                        self.iteration+=1
                    else:
                        self.lval=test
                        self.left=self.guess
                        self.iteration+=1
            elif type=='regula falsi':
                new=self.right - self.rval*(self.right-self.left)/(self.rval-self.lval)
#                    print self.left,new,self.right
#                    if self.left > new or self.right < new:
#                        new=(self.right + self.left)/2.0
                test=self.fn(new)
                matplotlib.pyplot.scatter(new,test)
                if abs(test)<self.tolerance:
                    self.guess=new
                    self.root=self.guess
                    return
                else:
                    if self.lval*test < 0:
                        self.rval=test
                        self.right=new
                        self.iteration+=1
                    else:
                        self.lval=test
                        self.left=new
                        self.iteration+=1
            elif type=='AB':
                new=self.right - self.rval*(self.right-self.left)/(self.rval-self.lval)
                print(self.left,self.right)
                print(self.lval,self.rval)
                test=self.fn(new)
                if old!=None:
                    if test*old < 0:                         # same bracket
                        if old==self.lval:
                            #m=1-self.rval*test
                            #if m<0:
                            m=.5
                            self.right=new
                            self.rval=test
                            self.lval*=m
                            old=self.lval
                        else:
                            #m=1-self.lval*test
                            #if m<0:
                            m=.5
                            self.left=new
                            self.lval=test
                            self.rval*=m
                            old=self.rval
                    else: 
                        if self.lval*test < 0:
                            self.rval=test
                            self.right=new
                            old=self.lval
                        else:
                            self.lval=test
                            self.left=new
                            old=self.rval
                else:
                    if self.lval*test < 0:
                            self.rval=test
                            self.right=new
                            old=self.lval
                    else:
                            self.lval=test
                            self.left=new
                            old=self.rval
                self.iteration+=1

        


class NewtonSolver(Solver):
    def __init__(self,function,deriv,guess=0,tolerance=1E-8,maxiter=20,verbose=False,history=False):
        Solver.__init__(self,function,guess,tolerance,maxiter)
        self.df=deriv
        self.residual=self.fn(guess)
        self.verbose=verbose
        self.history=history
        if history:
            self.sequence=np.array(guess)
            

        
    def state(self):
        print('Residual:',self.residual)
        Solver.state(self)
        
    def solve(self):
        
        while abs(self.residual) > self.tolerance and self.iteration < self.maxiter:  #if residual large or max iterations not reached, iterate.
            if self.verbose: #print state each iteration if verbose=True
                self.state()
            self.guess -= self.residual/self.df(self.guess) #Newton update (self.residual = f(x_n))
            if self.history:
                self.sequence=np.append(self.sequence,self.guess)
            self.residual=self.fn(self.guess) #compute residual
            self.iteration+=1  #increase iteration count
        return 

    
class SecantSolver(Solver):
    def __init__(self,function,guess=1,old=0,tolerance=1E-8,maxiter=20,verbose=False,history=False):
        Solver.__init__(self,function,guess,tolerance,maxiter)
        self.guess=guess
        self.residual=self.fn(guess)
        self.verbose=verbose
        self.history=history
        self.old=old
        if history:
            self.sequence=np.array((old,guess))
            

    def state(self):
        print('Residual:',self.residual)
        Solver.state(self)
        
    def solve(self):
        while abs(self.residual) > self.tolerance and self.iteration < self.maxiter:  #if residual large or max iterations not reached, iterate.
            if self.verbose: #print state each iteration if verbose=True
                self.state()
            temp = np.copy(self.guess)
            self.guess -= (self.guess-self.old)*self.residual/(self.residual-self.fn(self.old)) #Secant update (self.residual = f(x_n))
            if self.history:
                self.sequence=np.append(self.sequence,self.guess)
            self.residual=self.fn(self.guess) #compute residual
            self.iteration+=1  #increase iteration count
            #update old
            self.old=temp

        return 
                        

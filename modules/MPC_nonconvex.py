import numpy as np

from pyomo.environ import *
from pyomo.dae import *
import matplotlib.pyplot as plt

n_states = 3

class MPC:
    def __init__(self,rocket):
        self.rocket = rocket
        self.maxThrust = 20
        self.minThrust = 7

    def random(self):
        thrust = 1
        #random number between -1 and 1
        alpha = np.random.rand()*2-1
        beta = np.random.rand()*2-1
        self.rocket.input(thrust,alpha*0.1,beta*np.pi/2)
        
    def approxTf(self,initState):
        self.tf_guess = np.sqrt(2*initState[2]/9.81)

    def model(self):
        self.m = ConcreteModel()
        self.m.t = ContinuousSet(bounds=(0,self.tf))
        #self.m.tf = Var(domain=NonNegativeReals, initialize=self.tf_guess)
        self.m.x1 = Var(self.m.t, within=Reals)
        self.m.x2 = Var(self.m.t, within=Reals)
        self.m.x3 = Var(self.m.t, within=Reals)
        self.m.dx1 = DerivativeVar(self.m.x1, wrt=self.m.t)
        self.m.dx2 = DerivativeVar(self.m.x2, wrt=self.m.t)
        self.m.dx3 = DerivativeVar(self.m.x3, wrt=self.m.t)
        self.m.ddx1 = DerivativeVar(self.m.dx1, wrt=self.m.t)
        self.m.ddx2 = DerivativeVar(self.m.dx2, wrt=self.m.t)
        self.m.ddx3 = DerivativeVar(self.m.dx3, wrt=self.m.t)

        self.m.gamma = Var(self.m.t, within=NonNegativeReals, initialize=self.minThrust)

        self.m.u1 = Var(self.m.t, within=Reals)
        self.m.u2 = Var(self.m.t, within=Reals)
        self.m.u3 = Var(self.m.t, within=Reals)
    
    def constraints(self):
        tf = self.tf
        self.m.ode1 = Constraint(self.m.t, rule=lambda m, t: m.ddx1[t] == m.u1[t] - 0.05*m.dx1[t]*abs(m.dx1[t]))
        self.m.ode2 = Constraint(self.m.t, rule=lambda m, t: m.ddx2[t] == m.u2[t] - 0.05*m.dx2[t]*abs(m.dx2[t]))
        self.m.ode3 = Constraint(self.m.t, rule=lambda m, t: m.ddx3[t] == m.u3[t] - 9.81 - 0.05*m.dx3[t]*abs(m.dx3[t]))

        self.m.thrust_gamma1 = Constraint(self.m.t, rule=lambda m, t: m.u1[t]**2 + m.u2[t]**2 + m.u3[t]**2 <= self.maxThrust**2)
        self.m.thrust_gamma1 = Constraint(self.m.t, rule=lambda m, t: m.u1[t]**2 + m.u2[t]**2 + m.u3[t]**2 >= self.minThrust**2)
        
        self.m.thrust_angle = Constraint(self.m.t, rule=lambda m, t: m.u3[t] >= 0) #constraint on angle

        self.m.x_min = Constraint(self.m.t, rule=lambda m, t: m.x3[t] >= 0)

        #self.m.final1 = Constraint(rule=lambda m: m.x1[tf]**2 + m.x2[tf]**2 <= 5)
        self.m.final1 = Constraint(rule=lambda m: m.x1[tf]**2 + m.x2[tf]**2 <= 1)
        self.m.final3 = Constraint(rule=lambda m: m.x3[tf]**2 <= 1)

        self.m.final4 = Constraint(rule=lambda m: abs(m.dx1[tf]) <= 0.1)
        self.m.final5 = Constraint(rule=lambda m: abs(m.dx2[tf]) <= 0.1)
        self.m.final6 = Constraint(rule=lambda m: abs(m.dx3[tf]) <= 0.1)
        

    def inits(self,state):
        self.m.x1[0].fix(state[0])
        self.m.x2[0].fix(state[1])
        self.m.x3[0].fix(state[2])
        self.m.dx1[0].fix(state[3])
        self.m.dx2[0].fix(state[4])
        self.m.dx3[0].fix(state[5])
        
        

        
    def apply(self):
        J = lambda m,t: (m.u1[t]**2 + m.u2[t]**2 + m.u3[t]**2) +  0.5*(self.m.x1[t]**2 + self.m.x2[t]**2 + self.m.x3[t]**2) #change to absolute value
        self.m.J = Integral(self.m.t, wrt = self.m.t, rule = J) # + self.m.x1[1]**2 + self.m.x2[1]**2 + self.m.x3[1]**2

        self.m.obj = Objective(expr=self.m.J, sense=minimize)
        
        # transform and solve
        TransformationFactory('dae.finite_difference').apply_to(self.m, wrt=self.m.t, nfe=100)
        #discretizer = TransformationFactory('dae.collocation')
        #discretizer.apply_to(m2,wrt=m2.t,nfe=5000,ncp=3,scheme='LAGRANGE-RADAU')
        self.solver = SolverFactory('ipopt')
        #discretizer.reduce_collocation_points(m2,var=m.u,ncp=1,contset=m.t)
        self.solver.options['max_iter']= 10000 #number of iterations you wish
        self.solver.solve(self.m).write()
        

    def run(self,state):
        #self.tf = self.tf -1 
        self.model()
        self.constraints()
        self.inits(state)
        self.apply()
        return self.result()

    def result(self):
        tt = np.array([t for t in self.m.t])

        u1 = np.array([self.m.u1[t]() for t in self.m.t])
        u2 = np.array([self.m.u2[t]() for t in self.m.t])
        u3 = np.array([self.m.u3[t]() for t in self.m.t])

        x1 = np.array([self.m.x1[t]() for t in self.m.t])
        x2 = np.array([self.m.x2[t]() for t in self.m.t])
        x3 = np.array([self.m.x3[t]() for t in self.m.t])

        dx1 = np.array([self.m.dx1[t]() for t in self.m.t])
        dx2 = np.array([self.m.dx2[t]() for t in self.m.t])
        dx3 = np.array([self.m.dx3[t]() for t in self.m.t])

        return [u1,u2,u3,x1,x2,x3,dx1,dx2,dx3,tt]

""" m = ConcreteModel()
m.t = ContinuousSet(bounds=(0,1))
m.x = Var(m.t)
m.x_dot = DerivativeVar(m.x, wrt=m.t)
m.tf = Var(domain=NonNegativeReals,initialize=4)

J = lambda m,t: (t*m.x_dot[t] + (m.x_dot[t]/m.tf)**2)*m.tf
m.J2 = Integral(m.t, wrt = m.t, rule = J)

m.obj = Objective(expr=m.J2, sense=minimize)
m.x[0].fix(1)
m.x[1].fix(10)

# transform and solve
TransformationFactory('dae.finite_difference').apply_to(m, wrt=m.t, nfe=9000)
#discretizer = TransformationFactory('dae.collocation')
#discretizer.apply_to(m2,wrt=m2.t,nfe=5000,ncp=3,scheme='LAGRANGE-RADAU')
solver = SolverFactory('ipopt')
#discretizer.reduce_collocation_points(m2,var=m.u,ncp=1,contset=m.t)
solver.options['max_iter']= 10000 #number of iterations you wish
solver.solve(m).write() """
"""
discretizer = TransformationFactory('dae.collocation')
discretizer.apply_to(m,wrt=m.t,nfe=50,ncp=3,scheme='LAGRANGE-RADAU')
solver = SolverFactory('ipopt')
#solver.options['max_iter']= 10000 #number of iterations you wish
solver.solve(m).write()

"""
from sympy import symbols, Function, Derivative, Eq, Integer, Rational
from ..symbolic import Expectation, PDESystem, SymbolicPKF, FDModelBuilder, t

class Burgers(object):

    # Set the spatial coordinate system
    x = symbols('x')
    # Set the constants
    kappa = symbols('kappa')
    # Define the spatio-temporal scalar field
    u = Function('u')(t,x)

    dynamics = Eq(Derivative(u,t), -u*Derivative(u,x)+kappa*Derivative(u,x,2))

    pkf = SymbolicPKF(dynamics)

    g = pkf.fields[u].metric[0] # metric tensor
    s = pkf.fields[u].aspect[0] # aspect tensorc

    _p18_closure = Integer(3)*g**Integer(2)-Integer(2)*Derivative(g,x,2)
    p18_closure = {
        'metric': _p18_closure,
        'aspect': _p18_closure.subs(g,1/s).doit().expand(),
    }

    unclosed_term = list(pkf.unclosed_terms)[0]
    closed_pkf = SymbolicPKF(dynamics)
    closed_pkf.set_closure({unclosed_term:p18_closure['aspect']})

    def __init__(self, n=241):
        def factory():
            exec(FDModelBuilder(Burgers.closed_pkf.in_aspect, class_name='ClosedPKFBurgers').code)
            return ClosedPKFBurgers
        self.model = factory()(shape=(n,))

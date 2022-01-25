from unittest import TestCase
from sympkf.symbolic import Eq, omega, israndom
from sympy import symbols, Function

class TestIsRandom(TestCase):
    """
    Test of `israndom`

    """

    def test_Eq(self):
        # Test the case randomness of an equation depending of the 
        # randomness of its lhs, rhs
        t, x = symbols('t x')
        e = Function('e')(t,x,omega)
        e1 = Function('e_1')(t,x,omega)
        f = Function('f')(t,x)

        # Test heterogeneous situations        
        eq = Eq(e,0)        
        self.assertTrue(israndom(eq))
        eq = Eq(e,f)        
        self.assertTrue(israndom(eq))

        # Test homogeneous situations        

        ## 1. Deterministic situations
        eq = Eq(f,0)        
        self.assertFalse(israndom(eq))
        eq = Eq(f,f)        
        self.assertFalse(israndom(eq))
        ## 2. Stochastic situations
        eq = Eq(e,e1)        
        self.assertTrue(israndom(eq))






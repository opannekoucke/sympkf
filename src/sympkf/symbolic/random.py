import sympy
from sympy import Symbol, Function, Wild, Add, Mul, Pow, Derivative, Integral, S
from .util import Eq, Matrix

__all__ = [ 'Expectation', 'omega', 'Omega', 'israndom', 'Eq' ]



class Omega(Symbol):
    ''' Random set

    Only one instance of Omega is possible: omega
    '''

    _instance = 0

    def __new__(cls, name):
        cls._instance += 1
        if cls._instance > 1:
            raise Exception("Two instance of Omega not allowed")
        else:
            return super(Omega, cls).__new__(cls, name)

omega = Omega('omega')

debug = True

class Expectation(Function):
    """ Expectation operator.

    There is no automatic evaluation, offering the possibility
    do manage more correctly the computation.

    Parameters
    ==========

    arg : sympy object

    Examples
    ========

    >>> from sympy import Function, symbols
    >>> from sympkf.symbolic import Expectation, omega
    >>> x, t = symbols('x t')
    >>> f = Function('f')(x,t)
    >>> X = Function('X')(x,t,omega)
    >>> Expectation(X+f*X, evaluate=False)
    E(X+f*X)
    >>> Expectation(X+f*X)
    E(X) + f*E(X)
    >>> Expectation(X*Derivative(X*f,t))
    Derivative(f,t)*E(X**2)+f*E(X*Derivative(X,t))

    """

    nargs = 1
    _latex_name = '{\mathbb E}'
    is_linear = True

    def _eval_derivative(self, v):
        return Derivative(self,v, evaluate=False)


    @classmethod
    def eval(cls, arg, **hints):


        # Get hints keys value
        deep = hints['deep'] if 'deep' in hints else False
        inside = hints['inside'] if 'inside' in hints else False
        # - DEBUG -
        # doit = hints['doit'] if 'doit' in hints else True
        doit = hints['doit'] if 'doit' in hints else False

        """ #####  Nor more used
        if doit:
            #arg = arg.doit().expand()
            arg = arg.expand()
        """
        arg = arg.expand()

        # Idendité des valeurs certaines
        # BUG Expectation -1-
        #   should replace `arg.has(omega)` by `israndom(arg)̀
        if israndom(arg)==False:
            """
            if str(arg)[0]=='_':
                # Case where this is a dummy variable (in `Derivative` with ̀evaluate=True`)
                return cls(arg, evaluate=False)
            else:
                # Case of classical deterministic variable
                return arg
            """
            return arg
        # Linéarité
        # 1 - addition classique
        elif isinstance(arg, Add):
            # arg = first + others
            args = arg.args
            #Old:
            # first, others = args[0], args[1:]
            # return Add(cls(first), cls(Add(*others)))
            return Add(*[cls(arg) for arg in args])
        # 2 - multiplication par une quantité certaine
        elif isinstance(arg, Mul):
            # Objectif mettre arg.args sous la forme
            # arg.args =  scalar * random(omega)
            # pour pouvoir retourner a * E(b)
            # Note:
            #   Idempotence E o E = E is include in this case ! due to 'has' method
            scalar = Wild('scalar', properties=[lambda k: israndom(k) == False])
            random = Wild('random', properties=[lambda k: israndom(k)])

            out = arg.match(scalar * random)
            if out[scalar] != S.One:
                return Mul(out[scalar], cls(out[random]))
            else:
                return cls(arg, evaluate=False)
        # Commutativité par rapport aux opérateurs de dérivations time/space
        elif isinstance(arg, Derivative):
            return Derivative(cls(arg.args[0]), arg.args[1])
        # Commutativité par rapport aux opérateurs de intégraux time/space
        elif isinstance(arg, Integral):
            return Integral(cls(arg.args[0]), arg.args[1])
        #
        # Apply to particular structure
        #
        # -- Eq
        elif isinstance(arg, Eq):
            lhs, rhs = arg.args
            #return Eq(cls(lhs), cls(rhs))
            return arg.__class__(cls(lhs), cls(rhs))
        # -- Matrix
        elif isinstance(arg, (sympy.Matrix,sympy.ImmutableMatrix)):
            return Matrix(*arg.shape, [cls(expr) for expr in arg])
        else:
            return cls(arg, evaluate=False)

    def _has(self, pattern):
        # This function is important to check if an expression is random or not.
        #print(f" ----   _has({pattern})  -----")
        if isinstance(pattern, Omega):
            return False
        else:
            return self.args[0]._has(pattern)

    def _latex(self, printer):
        """
        Make the latex printing with standard \\mathbb E notation.
        :param printer:
        :return:

        :Example:

            >>> from sympkf.util.symbolic import Expectation, omega
            >>> from sympy import Function
            >>> X = Function('X')(omega)
            >>> type(Expectation(X))
            Expectation
            >>> latex(Expectation(X))
            '\\mathbb E\\left(X{\\left(\\omega \\right)}\\right)'

        """
        return self._latex_name+ r"\left(%s\right)" % printer._print(self.args[0])


class IncoherentEquality(Exception): pass

def israndom(arg, verbose=False):
    """ Test if `arg` is a random variable or not """
    ''' List of rules
    1.  If `arg` doesn't contains `omega` it should not be random
    2.  Expectation of an expression is non-random
    '''
    def lprint(expr):
        if verbose:
            print(expr)
    if arg.has(omega) is False:
        # Rule 2.
        lprint('israndom: check omega')
        return False
    elif isinstance(arg, Expectation):
        # Rule 1.
        lprint('israndom: Expectation')
        return False
    elif isinstance(arg, sympy.Eq):
        lhs = israndom(arg.args[0], verbose=verbose)
        rhs = israndom(arg.args[1], verbose=verbose)
        # Warning : it is possible that lhs!=rhs e.g. 
        # if rhs is random while lhs is not
        # example: 
        #     e(t,x,omega) = 0
        #   in this case, 'lhs' is random (there is omega) while 'rhs' 
        #   is deterministic (there is no omega)
        return any([lhs,rhs])
    elif isinstance(arg, Function):
        lprint('israndom: Function')
        return any( israndom(elm, verbose=verbose) for elm in arg.args )
    elif isinstance(arg, Symbol):
        lprint('israndom: Symbol')
        return isinstance(arg, Omega)
    elif isinstance(arg, (Derivative, Integral)):
        lprint(f'israndom: {type(arg)}')
        return israndom(arg.args[0], verbose=verbose)
    elif isinstance(arg, (Mul, Add, Pow)):
        lprint(f'israndom: {type(arg)}')
        return any(israndom(elm, verbose=verbose) for elm in arg.args)
    elif isinstance(arg, (sympy.Matrix, sympy.ImmutableMatrix)):
        return any(israndom(elm, verbose=verbose) for elm in arg)
    else:
        raise NotImplementedError

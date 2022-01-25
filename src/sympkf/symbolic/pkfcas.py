from .random import Expectation, omega
from .util import PDESystem, Eq, remove_eval_derivative, upper_triangle
from .tool import clean_latex_name
from .constants import t as time_symbol
import collections


from sympy import Derivative, symbols, Function, sqrt, Integer, Rational, Mul, Matrix
import sympy



#__all__ =['SymbolicPKF','Field','PartialOrderDerivative']

# Add Abs operator for positive function in order to evaluate absolute values of positive function
# Example:
# --------
#  >>> func = Function('f')(t,x)
#  >>> Abs(f)
#   |f(t,x)|
#  >>> pos_func = Function('f',positive=True)(t,x)
#  >>> Abs(f)
#   f(t,x)
def _eval_Abs(self):
    if self.is_positive:
        return self

setattr(Function, '_eval_Abs', _eval_Abs)

class SymbolicPKF(object):
    """ Parametric dynamics associated to a given system of evolution equations

    The workflow is the follow one:
      1. Calculation of the linearized dynamics using the Reynolds decomposition of random fields
      2. Calculation of the Reynolds averaged dynamics
      3. Calculation of the error dynamics
      4. Calculation of the variance dynamics
      5. Calculation of the normalized error dynamics
      6. Calculation of the metric dynamics
    """
    def __init__(self, pde_system, closure=None):

        if not isinstance(pde_system, PDESystem):
            pde_system = PDESystem(pde_system)

        self.pde_system = pde_system

        self.fields = collections.OrderedDict([(field, Field(field))
                                               for field in self.pde_system.prognostic_functions]
                                              )
        self.time_coordinate = self.pde_system.time_coordinate
        self.tl_system = None
        self._first_order_reynolds_system = None
        self._second_order_reynolds_system = None
        self._expectation_system = None
        self._tl_expectation_system = None
        self._error_system = None
        self._epsilon_system = None
        self._variance_system = None
        self._cross_covariance_system = None
        self._std_system = None
        self._metric_system = None

        self._system_in_metric = None
        self._system_in_diffusion = None
        self._system_in_aspect = None

        self._closed_system_in_metric = None
        self._closed_system_in_diffusion = None
        self._closed_system_in_aspect = None        

        self._internal_closure = None

        # external closure
        self._closure = {} if closure is None else closure

        # Introduce or modify labels for error/normalized error/variance/cross-covariance
        self._label_covariance = "V_"

    def _compute_reynolds_system(self):
        """ Computation of Reynolds system at the second and the first orders
        :return:
        """

        # Compute the linear reynolds system
        first_order_reynolds = []
        second_order_reynolds = []

        epsilon = symbols('\\varepsilon')
        subs_reynolds = {            
            field: Expectation(self.fields[field].random) + epsilon * self.fields[field].error            
            for field in self.fields
        }

        # bug-Expectation-hotfix 5595950431936087884
        #
        tmp_func = Function('tmp-func')
        # /bug

        for equation in self.pde_system.equations:
            # 1) Reynolds substitution
            equation = equation.subs(subs_reynolds).doit().expand()
            lhs, rhs = equation.args

            # 2) Calculation of the Second order
            # bug-Expectation-hotfix 5595950431936087884
            #
            rhs = rhs.subs(Expectation, tmp_func)
            #
            # / bug

            taylor = rhs.series(epsilon, 0, 3).removeO()

            # bug-Expectation-hotfix 5595950431936087884
            #            
            taylor = taylor.subs(tmp_func, Expectation)            
            # /bug

            equation = Eq(lhs, taylor).expand()

            # 3) Computation of the second-order Taylor's expansion
            second_order_reynolds.append(equation.subs(epsilon, Integer(1)))

            # 4) Computation of the first-order Taylor's expansion
            equation = equation.subs(epsilon*epsilon, Integer(0))
            first_order_reynolds.append(equation.subs(epsilon, Integer(1)))

        self._first_order_reynolds_system = first_order_reynolds
        self._second_order_reynolds_system = second_order_reynolds

    @property
    def first_order_reynolds_system(self):
        """ Create the linear dynamics by including the Reynolds decomposition of random fields.
        :return:
        """

        if self._first_order_reynolds_system is None:
            self._compute_reynolds_system()

        return self._first_order_reynolds_system

    @property
    def second_order_reynolds_system(self):
        """ Calculation of the second-order Reynolds equations """

        if self._second_order_reynolds_system is None:
            self._compute_reynolds_system()

        return self._second_order_reynolds_system

    @property
    def expectation_system(self):
        """ Calculate the expectation dynamics of averaged fields """

        if self._expectation_system is None:

            # dic for substitution error -> variance**(1/2)*epsilon
            subs_error = {
                mfield.error: mfield.variance ** Rational(1 , 2) * mfield.epsilon for mfield in self.fields.values()
            }

            expectation_system = []
            for equation in self.second_order_reynolds_system:
                # 1) Computation of the expectation
                equation = Expectation(equation)
                # 2) Substitute error in Variance/Corrélation
                equation = equation.subs(subs_error).doit().expand()
                # 3) Apply internal closure
                equation = equation.subs(self.internal_closure).doit()

                expectation_system.append(equation)

            self._expectation_system = expectation_system

        return self._expectation_system

    @property
    def tl_expectation_system(self):
        """ Calculate the TL expectation system for later use in error_sytem"""
        if self._tl_expectation_system is None:
            tl_expectation_system = []
            for equation in self.first_order_reynolds_system:
                # 1) Computation of the expectation
                equation = Expectation(equation)
                # 2) Apply internal closure
                equation = equation.subs(self.internal_closure).doit()
                tl_expectation_system.append(equation)

            self._tl_expectation_system = tl_expectation_system

        return self._tl_expectation_system

    @property
    def error_system(self):
        """ Error equations are difference between reynolds equation and expectation equation """
        if self._error_system is not None:
            return self._error_system

        else:
            # Compute the error system
            error_system = []

            for linear_reynolds_equation, expectation_equation in zip(self.first_order_reynolds_system,
                                                                      self.tl_expectation_system):
                error_equation = linear_reynolds_equation - expectation_equation
                error_system.append(error_equation.expand())

            self._error_system = error_system

            return self._error_system

    @property
    def variance_system(self):

        if self._variance_system is not None:
            return self._variance_system

        else:
            # Compute the variance system
            subs_error_trend = {equation.args[0]: equation.args[1] for equation in self.error_system}
            subs_error_to_epsilon = {
                self.fields[field].error: self.fields[field].epsilon * sqrt(self.fields[field].variance)
                for field in self.fields
            }
            # Built the system of equation
            variance_system = []
            for field in self.fields:
                variance = self.fields[field].variance
                error = self.fields[field].error

                definition = Eq(
                    Derivative(variance, time_symbol),
                    Expectation(Derivative(error ** Rational(2), time_symbol), evaluate=False)
                )

                lhs, rhs = definition.args

                rhs = rhs.doit()
                rhs = rhs.subs(subs_error_trend)
                rhs = rhs.expand()
                rhs = rhs.subs(subs_error_to_epsilon).doit()
                rhs = rhs.expand()

                # Substitution with internal_closure and higher terms
                rhs = self._apply_internal_closure(rhs)

                equation = Eq(lhs, rhs)

                variance_system.append(equation)

            self._variance_system = variance_system

            return self._variance_system

    @property
    def epsilon_system(self):
        """
        Computation of the trend of epsilon (dynamics of the normalized error). This is used in the computation
        of the dynamics of the metric tensor: the substitution with epsilon is easier than with the error.
        :return:
        """


        if self._epsilon_system is not None:
            return self._epsilon_system

        else:            

            subs_error_trend = {
                eq.args[0]: eq.args[1] for eq in self.error_system
            }
            subs_variance_trend = {
                eq.args[0]: eq.args[1] for eq in self.variance_system
            }
            subs_error = {
                self.fields[field].error: self.fields[field].epsilon * sqrt(self.fields[field].variance)
                for field in self.fields
            }

            epsilon_system = []
            for field in self.fields:
                epsilon = self.fields[field].epsilon
                error = self.fields[field].error
                sqrt_variance = sympy.sqrt(self.fields[field].variance)

                lhs = Derivative(epsilon, time_symbol)

                rhs = Derivative(error / sqrt_variance, time_symbol).doit()
                rhs = rhs.subs(subs_error_trend).doit()
                rhs = rhs.subs(subs_variance_trend).doit()
                rhs = rhs.subs(subs_error).doit()

                equation = Eq(lhs, rhs.expand())

                epsilon_system.append(equation)

            self._epsilon_system = epsilon_system

            return self._epsilon_system

    @property
    def std_system(self):
        """ Calculation of the dynamics of the standard deviations

        Comment:
            This system is not include in the calculation workflow of the parametric dynamics
        """


        if self._std_system is not None:
            return self._std_system

        else:
            # Compute the std system            

            subs_variance = {
                self.fields[field].variance: self.fields[field].std ** Integer(2)
                for field in self.fields
            }

            std_system = []
            for equation in self.variance_system:
                # Replace variance by std
                equation = equation.subs(subs_variance).doit()
                # Extract trends
                trends = equation.args[0].atoms(Derivative)
                for trend in trends:
                    if time_symbol == trend.args[1]:
                        break
                equation = equation.isolate(trend)
                std_system.append(equation)

            self._std_system = std_system

            return self._std_system

    @property
    def multivariate_couples(self):
        """ Return couples of multivariate fields """
        fields = list(self.fields.keys())
        couples = ( (f1, f2) for i,f1 in enumerate(fields) for j,f2 in enumerate(fields) if i<j )
        return couples


    @property
    def cross_covariance_system(self):
        """ Return the dynamics of the cross-covariance in multivariate situations """
        if self._cross_covariance_system is not None:
            return self._cross_covariance_system

        elif len(list(self.fields.keys()))==1:
            # Univariate situation            
            return None

        else:
            # Compute the cross-covariance system

            # Set substution dictionnary
            subs_error_trends = {equation.args[0]: equation.args[1] for equation in self.error_system}
            subs_error_to_epsilon = {
                self.fields[field].error: self.fields[field].epsilon * sqrt(self.fields[field].variance)
                for field in self.fields
            }
        

            # Set the cross-covariance meta-data and compute the cros-covariance dynamics 
            self._cross_covariance_system = []
            
            for couple in self.multivariate_couples:
                # 1) Extract fields and meta-data
                f1, f2 = couple        
                mf1, mf2 = self.fields[f1], self.fields[f2] 

                # 2) extract error fields
                e1, e2 = mf1.error, mf2.error                                                
                V12 = self.internal_closure[Expectation(e1*e2)]

                # 3) Definition and computation of the dynamics
                lhs = Derivative(V12, time_symbol)
                rhs = Expectation(Derivative(e1*e2, time_symbol).doit()).subs(subs_error_trends)
                rhs = rhs.subs(subs_error_to_epsilon).doit()
                #.. todo ??:
                #  should include substitution from 'epsilon' in place of 'error'

                rhs = self._apply_internal_closure(rhs)                

                # 4) update of the cross-covariance system
                self._cross_covariance_system.append( Eq(lhs,rhs.expand()) )
                
            return self._cross_covariance_system

    @property
    def metric_system(self):
        """ Return the dynamics of the metric using the variances (not the stds) """
        if self._metric_system is not None:
            return self._metric_system

        else:

            #t = self.time_coordinate

            subs_epsilon_trends = {
                eq.args[0]: eq.args[1] for eq in self.epsilon_system
            }

            # Compute the pre-system
            metric_system = []
            for meta_field in self.fields.values():
                # Create the system for each component of the metric tensors
                for i,xi in enumerate(meta_field.spatial_coordinates):
                    for j, xj in enumerate(meta_field.spatial_coordinates):
                        if j<i:
                            continue

                        # Set the lhs: D_t g_ij
                        lhs = Derivative(meta_field.metric_func(i,j), time_symbol)

                        # Compute the rhs: E[ D_t(D_i eps D_j eps)]
                        #  - Definition of the rhs
                        rhs = Expectation(
                                    Derivative(
                                        Derivative(meta_field.epsilon,xi)*Derivative(meta_field.epsilon,xj)
                                    ,time_symbol).doit()
                                    )
                        # Substitutes the trends of error
                        rhs = rhs.subs(subs_epsilon_trends).doit()

                        # Subs usual statistic properties of epsilon
                        rhs = self._apply_internal_closure(rhs)

                        # Forms the equation
                        equation = Eq(lhs, rhs.expand())

                        metric_system.append(equation)

            self._metric_system = metric_system
            return self._metric_system

    @property
    def in_metric(self):
        """ Return de the pkf dynamics is term of the variance/metric tensor """
        if self._system_in_metric is None:

            # Compute the expectation system
            full_system_in_metric = []

            systems = [
                        self.expectation_system, 
                        self.variance_system, 
                        self.cross_covariance_system, 
                        self.metric_system
                    ] 
            
            for system in systems:

                if system is None: # to handle the multivariate situation: None is for univariate
                    continue
                
                # 1. Closes the system (by default closure is the empty dictionary {})
                if self.closure != {}:
                    closed_system = []
                    for equation in system:
                        equation = equation.subs(self.closure).expand()
                        closed_system.append(equation)
                    system = closed_system

                # 2. Clean from expectation of random fields
                full_system_in_metric += self._clean_system_from_expectation(system)

            self._system_in_metric = full_system_in_metric

        return self._system_in_metric

    @property
    def in_diffusion(self):
        """ Return de the pkf dynamics is term of the variance/diffusion tensor

        Description
        -----------

            The derivation relies on the identity $\nu g = I/2$ whose trend is $(\dot\nu)g+\nu(\dot g)=0$
            so that $\dot \nu = -\nu(\dot g)g^{-1}$. Again, with $\nu g =I/2$ leading to $g^{-1} = 2\nu$
            it results that
            $$\dot \nu = -2\nu(\dot g)\nu.$$

        """

        if self._system_in_diffusion is None:

            t = self.time_coordinate

            # 1. Set dictionary for substitution

            #  1.1 Create the substitution dictionary for migration : metric -> diffusion
            metric_to_diffusion = collections.OrderedDict()
            for mfield in self.fields.values():
                metric = upper_triangle(mfield.metric)
                metric_in_diffusion = mfield.diffusion.inv() * Rational(1 , 2)
                metric_in_diffusion = upper_triangle(metric_in_diffusion)
                for metric_ij, diffusion_ij in zip(metric, metric_in_diffusion):
                    metric_to_diffusion[metric_ij] = diffusion_ij

            #  1.2 Dictionary for the metric trends
            subs_metric_trend = {}
            for equation in self._apply_closure(self.metric_system):
                lhs, rhs = equation.args
                subs_metric_trend[lhs] = rhs

            #  2. Migration of expectation and variance systems
            diffusion_system = []
            for system in [self.expectation_system, self.variance_system]:
                # -1- apply external closure
                system = self._apply_closure(system)
                # -2- switch from metric to diffusion
                for equation in system:
                    equation = equation.subs(metric_to_diffusion)
                    diffusion_system.append(equation)

            # 3. Computation of the system at a symbolic level
            #     forms the equation $$ \pdt \nu = - 2\nu \pdt g \nu $$
            #     The computation of the system is made as a loop over univariate fields
            for mfield in self.fields.values():
                # Extract tensors
                diffusion = mfield.diffusion
                metric = mfield.metric

                # Computation of the rhs: $- 2\nu \pdt g \nu$
                trend_metric = Derivative(metric, t).doit()
                rhs = -Integer(2)*diffusion*trend_metric*diffusion
                rhs = rhs.doit()

                # Computation of the lhs: $\pdt \nu$
                lhs = Derivative(diffusion, t).doit()

                # Set the system by substituting terms
                for lhs_term, rhs_term in zip(upper_triangle(lhs), upper_triangle(rhs)):
                    # Replace metric trend by its value
                    rhs_term = rhs_term.subs(subs_metric_trend)
                    rhs_term = rhs_term.doit()
                    rhs_term = rhs_term.simplify()
                    rhs_term = rhs_term.expand()

                    # Replace metric terms by their values
                    rhs_term = rhs_term.subs(metric_to_diffusion)
                    rhs_term = rhs_term.doit()
                    rhs_term = rhs_term.simplify()
                    rhs_term = rhs_term.expand()

                    # Set the equation
                    equation = Eq(lhs_term, rhs_term)
                    diffusion_system.append(equation)

            # 3. Clean Expectation of fields
            diffusion_system = self._clean_system_from_expectation(diffusion_system)

            self._system_in_diffusion = diffusion_system

        return self._system_in_diffusion

    @property
    def in_aspect(self):
        """ Return de the pkf dynamics is term of the variance/aspect tensor

        Description
        -----------

            The derivation relies on the identity $\nu \bs = I$ whose trend is 
            $(\dot\bs)g+\bs(\dot g)=0$
            so that $\dot \bs = -\bs(\dot g)g^{-1}$. Again, with $\bs g =I$ 
            leading to $g^{-1} = \bs$
            it results that
            $$\dot \bs = -\bs(\dot g)\bs.$$

        """

        if self._system_in_aspect is None:            

            # 1. Set dictionary for substitution

            #  1.1 Create the substitution dictionary for migration : metric -> aspect
            metric_to_aspect = collections.OrderedDict()
            for mfield in self.fields.values():
                metric = upper_triangle(mfield.metric)
                # add 'aspect' in mfield
                metric_in_aspect = mfield.aspect.inv()
                metric_in_aspect = upper_triangle(metric_in_aspect)
                for metric_ij, aspect_ij in zip(metric, metric_in_aspect):
                    metric_to_aspect[metric_ij] = aspect_ij

            #  1.2 Dictionary for the metric trends
            subs_metric_trend = {}
            for equation in self._apply_closure(self.metric_system):
                lhs, rhs = equation.args
                subs_metric_trend[lhs] = rhs

            #  2. Migration of expectation and variance systems
            aspect_system = []

            systems = [ self.expectation_system, 
                        self.variance_system, 
                        self.cross_covariance_system
                    ]

            for system in systems:
                
                if system is None: # to handle the multivariate situation: None is for univariate
                    continue

                # -1- apply external closure
                system = self._apply_closure(system)
                # -2- switch from metric to diffusion
                for equation in system:
                    equation = equation.subs(metric_to_aspect)
                    aspect_system.append(equation)

            # 3. Computation of the system at a symbolic level
            #     forms the equation $$ \pdt \bs = - \bs \pdt g \bs $$
            #     The computation of the system is made as a loop over univariate fields

            #t = self.time_coordinate

            for mfield in self.fields.values():
                # Extract tensors
                aspect = mfield.aspect
                metric = mfield.metric

                # Computation of the rhs: $- \bs \pdt g \bs$
                trend_metric = Derivative(metric, time_symbol).doit()
                rhs = - aspect * trend_metric * aspect
                rhs = rhs.doit()

                # Computation of the lhs: $\pdt \bs$
                lhs = Derivative(aspect, time_symbol).doit()

                # Set the system by substituting terms
                for lhs_term, rhs_term in zip(upper_triangle(lhs), upper_triangle(rhs)):
                    # Replace metric trend by its value
                    rhs_term = rhs_term.subs(subs_metric_trend)
                    rhs_term = rhs_term.doit()
                    rhs_term = rhs_term.simplify()
                    rhs_term = rhs_term.expand()

                    # Replace metric terms by their values
                    rhs_term = rhs_term.subs(metric_to_aspect)
                    rhs_term = rhs_term.doit()
                    rhs_term = rhs_term.simplify()
                    rhs_term = rhs_term.expand()

                    # Set the equation
                    equation = Eq(lhs_term, rhs_term)
                    aspect_system.append(equation)

            # 3. Clean Expectation of fields
            aspect_system = self._clean_system_from_expectation(aspect_system)

            self._system_in_aspect = aspect_system

        return self._system_in_aspect


    @property
    def closure(self):
        return self._closure

    def set_closure(self, closure):
        # 1. Update the closure
        self._closure.update(closure)

        # 2. Reset systems in metric/diffusion/aspect
        self._system_in_metric = None
        self._system_in_diffusion = None
        self._system_in_aspect = None

    def _clean_system_from_expectation(self, system):
        """ Eliminate expectation of random fields from equation to simplify the representation and to prepare
        the translation in computational codes """

        clean_expectation = {}
        for mfield in self.fields.values():
            clean_expectation[Expectation(mfield.random)] = mfield.value

        new_system = []
        for equation in system:
            new_system.append( equation.subs(clean_expectation))

        return new_system

    @property
    def internal_closure(self):
        if self._internal_closure is None:

            self._internal_closure = {}

            # 1. Set univariate closure
            for meta_field in self.fields.values():
                self._internal_closure.update(meta_field.internal_closure)
            
            # 2. Set multivariate closure (only in multivariate situations)           
            for couple in self.multivariate_couples:
                # 1) Extract fields and meta-data
                f1, f2 = couple        
                mf1, mf2 = self.fields[f1], self.fields[f2] 

                # 2) extract error fields
                e1, e2 = mf1.error, mf2.error
                V1, V2 = mf1.variance, mf2.variance
                std1, std2 = sqrt(V1), sqrt(V2)
                eps1, eps2 = mf1.epsilon, mf2.epsilon

                # 3) Definition of the cross_variance
                
                # 3.a) Extract he cross covariance label
                V12 = self.get_covariance(f1,f2)

                # 3.c) Update internal closure
                self._internal_closure[Expectation(e1*e2)] = V12
                self._internal_closure[Expectation(eps1*eps2)] = V12/(std1*std2)

        return self._internal_closure
    
    def get_covariance(self, f1, f2):
        if all([field in self.fields for field in [f1,f2]]):
            # 1. Get associated metafields
            mf1 = self.fields[f1]
            mf2 = self.fields[f2]

            # 2. Selection of the coordinates
            # .. todo: 
            #   Modify the selection of the coordinates to account of two-point covariances between surface / volumique fields
            # this could be made from the cup product of the coordinates mf1.coordinates and mf2.coordinates
            # e.g. f1(t,x) f2(t,x,y) => V12(t,x,y) ??            
            cf1 = mf1.coordinates
            cf2 = mf2.coordinates
            assert cf1==cf2, ValueError("f1 and f2 have different coordinate system")
            coordinates = cf1

            return Function(self._label_covariance+f1.name+f2.name)(*coordinates)
        else:
            raise ValueError("f1 or f2 are not prognostic fields")

    @property
    def subs_tree(self):
        """
        :return: substitution tree where only the internal closure is used
        and which corresponds to the dictionnary of univariate terms E[D^alpha eps D^beta eps] 
        given as function of terms in E[eps D^gamma eps], for orders larger then 1.
        """
        subs_tree = {}
        for meta_field in self.fields.values():
            # Extract the tree
            meta_subs_tree = meta_field.subs_tree()
            # Close the tree from the defined internal_closure
            closed_subs_tree = {key:value.subs(self.internal_closure).doit() 
                                for key, value in meta_subs_tree.items()}

            # Update the tree
            subs_tree.update(closed_subs_tree)

        return subs_tree

    @staticmethod
    def check_univariate(expr):
        """ Check the univariate terms from an expression """

        expectations = expr.atoms(Expectation)

        univariates = {}

        for term in expectations:
            univariate = UnivariateTerm.is_univariate(term)
            if univariate is not None:
                function = univariate.function
                if function in univariates:
                    univariates[function].add(univariate)
                else:
                    univariates[function] = {univariate}

        return univariates

    def _apply_closure(self, system):
        """ Apply external closure """
        if self.closure == {}:
            return system
        else:
            closed_system = []
            for equation in system:
                equation = equation.subs(self.closure).expand()
                closed_system.append(equation)
            return closed_system

    def _apply_internal_closure(self, rhs):
        """ Apply the internal_closure on an expression (generally the rhs of an equation) """

        epsilon_to_mfields = {}
        for field,mfield in self.fields.items():
            epsilon_to_mfields[mfield.epsilon] = mfield

        # -1- Get univariate terms
        univariates = self.check_univariate(rhs)

        # -2- Retain epsilon's that are in self.fields !!
        for epsilon in univariates:
            if epsilon not in epsilon_to_mfields:
                univariates.pop(epsilon)

        # -3- Compute max degree for each univariate terms, by functions
        max_degrees = {}
        for epsilon in univariates:
            max_degrees[epsilon] = max([univariate.degree for univariate in univariates[epsilon]])

        # -4- Subs univariate terms E[D^alpha eps D^beta eps] by terms in E[eps D^gamma eps]
        #  ---- Replace only present terms ..
        for epsilon in univariates:
            max_degree = max_degrees[epsilon]
            subs_tree = epsilon_to_mfields[epsilon].subs_tree(max_degree)
            closed_subs_tree = {key:value.subs(self.internal_closure).doit() for key, value in subs_tree.items()}
            # can extract only required terms..
            #  -- but has to handle terms like E(eps D^beta eps) which should not be replace and in tre..
            rhs = rhs.subs(closed_subs_tree)

        # -3- Closure
        rhs = rhs.subs(self.internal_closure)
        return rhs

    @property
    def unclosed_terms(self):
        unclosed_terms = set()

        # -1- Search in variance/metric system
        systems = [self.variance_system, self.metric_system]
        for system in systems:
            for equation in system:
                unclosed_terms.update(equation.args[1].atoms(Expectation))

        # -2- Eliminates Expectation(field)
        for field,mfield in self.fields.items():
            unclosed_terms.discard(Expectation(mfield.random))

        return unclosed_terms

class UnivariateTerm(object):
    """
    Handle terms in E[ D^alpha \eps D^beta \eps]
    """

    """ ..todo ::  details what is univariate terms
         Univariate terms are of the form E(Dx^alpha eps Dx^beta eps)
    """
    def __init__(self, term, function, alpha, beta, degree):
        self.term = term
        self.function = function
        self.alpha = alpha
        self.beta = beta
        self.degree = degree

    def __repr__(self):
        return f"Univariate term: {self.function},{self.alpha}, {self.beta}, {self.degree}"

    @classmethod
    def is_univariate(cls, term):

        # Get derivatives
        derivatives = term.atoms(Derivative)

        # Get functions
        functions = set()
        for derivative in derivatives:
            functions.update(derivative.atoms(Function))

        if len(functions) == 1:  # Univariate term
            function = functions.pop()

            if len(derivatives) == 1: # E[eps D^beta eps] or E[(D^alpha eps)**2] or E[(D^alpha eps)**k] k>2
                derivative = derivatives.pop()
                if term is Expectation(function * derivative): # E[eps D^beta eps]
                    alpha = 0
                    beta = derivative.args[1:]
                    degree = derivative.derivative_count
                    return cls(term, function, alpha, beta, degree)
                elif term is Expectation(derivative*derivative): #E[(D^alpha eps)**2]
                    alpha = derivative.args[1:]
                    beta = derivative.args[1:]
                    degree = 2*derivative.derivative_count
                    return cls(term, function, alpha, beta, degree)
                else:
                    # E[(D^alpha eps)**k] k>2
                    return None

            elif len(derivatives) == 2:
                if term is Expectation(Mul(*derivatives)):
                    # -1- Compute the total degree
                    degree = 0
                    for derivative in derivatives:
                        degree += derivative.derivative_count
                    # -2- Extract the two derivatives
                    alpha, beta = derivatives
                    alpha = alpha.args[1:]
                    beta = beta.args[1:]
                    return cls(term, function, alpha, beta, degree)
                else:
                    return None
            else:
                return None
        else:
            return None


def gamma_def(epsilon, k, i, j):    
    return Expectation(Derivative(epsilon, k) * Derivative(epsilon, i, j))

def gamma_subs(metric, k, i, j):
    return Rational(1 , 2) * (Derivative(metric(k, j), i) + Derivative(metric(i, k), j) \
                              - Derivative(metric(i, j), k))

def skewness_def(epsilon, k, i, j):
    return Expectation(epsilon * Derivative(epsilon, k, i, j))

def skewness_subs(metric, k, i, j):
    return -Rational(1 , 2) * (Derivative(metric(k, j), i) + Derivative(metric(i, k), j) \
                               + Derivative(metric(i, j), k))

class Field(object):
    def __init__(self, field):
        self.value = field
        self.code = clean_latex_name(field.func)

        self.coordinates = field.args
        self.spatial_coordinates = tuple([coord for coord in self.coordinates if coord is not time_symbol])

        self.coords_code = tuple(clean_latex_name(coord) for coord in self.coordinates)
        self.spatial_coords_code = tuple(clean_latex_name(coord) for coord in self.spatial_coordinates)

        # Associated random fields
        self.random = Function(str(field.func))(*field.args, omega)
        self.epsilon = Function('{\\varepsilon_{' + self.code + '}}')(*self.coordinates, omega)
        self.error = Function('{e_{' + self.code + '}}')(*self.coordinates, omega)

        # Associated statistics
        # -- Variance field
        self.variance = Function('{V_{' + self.code + '}}', positive=True)(*self.coordinates)

        # -- Standard deviation field
        self.std = Function('{\\sigma_{' + self.code + '}}', positive=True)(*self.coordinates)

        # -- Tensor fields
        shape_tensor = 2 * (len(self.spatial_coordinates),)
        self.metric = Matrix(*shape_tensor, self.metric_func)
        self.diffusion = Matrix(*shape_tensor, self.nu_func)
        self.aspect = Matrix(*shape_tensor, self.aspect_func)

        # --trends
        self.trends = {
                    'field':Derivative(self.value,time_symbol),
                    'variance':Derivative(self.variance,time_symbol),
                    'error': Derivative(self.error, time_symbol),
                    'epsilon': Derivative(self.epsilon, time_symbol),
                    'metric': Derivative(self.metric, time_symbol),
                    'diffusion': Derivative(self.diffusion, time_symbol),
                    'aspect': Derivative(self.aspect, time_symbol),
        }

        self.subs_tree = UnivariateTree(self.epsilon, self.spatial_coordinates)

        self._internal_closure = None

    @property
    def internal_closure(self):

        if self._internal_closure is None:
            # Computation of the default internal_closure
            closure = {}
            # -0- error is centered
            closure.update({Expectation(self.error): Integer(0)})
            closure.update({Expectation(self.error ** Integer(2)): self.variance})
            # -1- epsilon is centered
            closure.update({Expectation(self.epsilon): Integer(0)})
            # -2- epsilon is normalized
            closure.update({Expectation(self.epsilon ** Integer(2)): Integer(1)})
            # -3- correlation is flat
            closure.update({Expectation(self.epsilon * Derivative(self.epsilon, coord)): Integer(0)
                            for coord in self.spatial_coordinates})
            # -4- metric is labeled
            closure.update({
                -self.metric_definition(i, j): -self.metric_func(i, j)
                for i in range(len(self.spatial_coordinates)) for j in range(i, len(self.spatial_coordinates))}
            )

            # -5- skewness is function of the metric
            skewness_closure = {}
            metric_func = lambda xi, xj: self.metric_func(self.spatial_coordinates.index(xi),
                                                          self.spatial_coordinates.index(xj))
            for partial_order in PartialOrderDerivative.all_of_degree(self.spatial_coordinates, 3):
                skewness_closure[skewness_def(self.epsilon, *partial_order.as_sequence)] = skewness_subs(metric_func,
                                                                                                         *partial_order.as_sequence)

            closure.update(skewness_closure)

            self._internal_closure = closure

        return self._internal_closure

    def index_code(self,i,j):
        if j<i:
            i,j = j,i
        return self.spatial_coords_code[i] + self.spatial_coords_code[j]

    def metric_func(self,i,j):
        return Function('{g_{' + self.code + ',' + self.index_code(i, j) + '}}', real=True)(*self.coordinates)

    def metric_definition(self,i,j):
        return -Expectation(self.epsilon * Derivative(self.epsilon, self.spatial_coordinates[i], self.spatial_coordinates[j]))

    def nu_func(self,i,j):
        return Function('{\\nu_{' + self.code + ',' + self.index_code(i, j) + '}}', real=True)(
            *self.coordinates)

    def aspect_func(self,i,j):
        return Function('{s_{' + self.code + ',' + self.index_code(i, j) + '}}', real=True)(
            *self.coordinates)

class PartialOrderDerivative(object):
    """ Handler for partial order derivatives """

    def __init__(self, coordinates, partial_orders):
        self._coordinates = coordinates
        self._as_tuple = self._builder(partial_orders)
        self._degree = sum(self._as_tuple)

    @property
    def coordinates(self):
        return self._coordinates

    @property
    def as_tuple(self):
        return self._as_tuple

    @property
    def as_couples(self):
        return [(coord, order) for coord, order in zip(self._coordinates, self._as_tuple)]

    @property
    def as_sequence(self):
        # format du type [x0,x1,x1]
        sequence = []
        for x_i, alpha_i in zip(self._coordinates, self._as_tuple):
            if alpha_i > 0:
                sequence += alpha_i * [x_i]
        return sequence

    def _builder(self, partial_orders):

        builder_names = ['tuple', 'couples', 'sequence']

        test = [getattr(self, 'is_as_' + name)(partial_orders) for name in builder_names]

        if any(test):
            name = builder_names[test.index(True)]
            return getattr(self, '_set_from_' + name)(partial_orders)
        else:
            raise TypeError(f"{partial_orders} is not a valid derivative order for coordinates {self._coordinates}")

    def is_as_tuple(self, partial_orders):

        if len(partial_orders) != len(self._coordinates):
            return False

        if isinstance(partial_orders, tuple):
            test = [isinstance(alpha_i, int) for alpha_i in partial_orders]
            return False not in test

        return False

    def is_as_couples(self, partial_orders):
        # Can not be empty (ie can not be [])
        if partial_orders == []:
            return False

        try:
            # can be a simple couple (x_i,alpha_i)
            if isinstance(partial_orders, tuple):
                xi, alpha_i = partial_orders
                return xi in self._coordinates and isinstance(alpha_i, int)
            # can be a list of as_couples [.. (x_i,alpha_i) ..]
            if isinstance(partial_orders, list):
                test = [xi in self._coordinates and isinstance(alpha_i, int) for xi, alpha_i in partial_orders]
                return False not in test
        except:
            return False

    def is_as_sequence(self, partial_orders):
        if partial_orders == []:
            return True
        test = [xi in self._coordinates for xi in partial_orders]
        return False not in test

    def _set_from_tuple(self, partial_orders_tuple):
        return partial_orders_tuple

    def _set_from_couples(self, couples):
        orders = len(self._coordinates) * [0]
        if isinstance(couples, tuple):
            xi, alpha_i = couples
            orders[self._coordinates.index(xi)] += alpha_i

        else:
            for xi, alpha_i in couples:
                orders[self._coordinates.index(xi)] += alpha_i

        return tuple(orders)

    def _set_from_sequence(self, sequence):
        return tuple([sequence.count(xi) for xi in self._coordinates])

    def __eq__(self, rhs):
        return self._as_tuple == rhs._as_tuple and self._coordinates == rhs._coordinates

    def new(self, partial_orders):
        return PartialOrderDerivative(self._coordinates, partial_orders)

    def copy(self):
        return PartialOrderDerivative(self._coordinates, self._as_tuple)

    @classmethod
    def all_of_degree(cls, coordinates, derivative_order):
        """ Return partial order derivative of all derivative at a given degree

        Description
        -----------

         The algorithm employs dynamical programming based on sets to avoid duplicate outputs.

         Each generation is computed from the previous by moving the partial order of the first coordinate toward the others.

        Example
        -------

          >>> coords = symbols(' '.join(['x'+str(i) for i in range(3)]))
          >>> for index in PartialOrderDerivative.all_of_degree(coords,4):
          >>>       print(index.as_couples)
            [(x0, 4.0), (x1, 0.0), (x2, 0.0)]
            [(x0, 3.0), (x1, 0.0), (x2, 1.0)]
            [(x0, 3.0), (x1, 1.0), (x2, 0.0)]
            [(x0, 2.0), (x1, 0.0), (x2, 2.0)]
            [(x0, 2.0), (x1, 1.0), (x2, 1.0)]
            [(x0, 2.0), (x1, 2.0), (x2, 0.0)]
            [(x0, 1.0), (x1, 2.0), (x2, 1.0)]
            [(x0, 1.0), (x1, 0.0), (x2, 3.0)]
            [(x0, 1.0), (x1, 3.0), (x2, 0.0)]
            [(x0, 1.0), (x1, 1.0), (x2, 2.0)]
            [(x0, 0.0), (x1, 1.0), (x2, 3.0)]
            [(x0, 0.0), (x1, 4.0), (x2, 0.0)]
            [(x0, 0.0), (x1, 2.0), (x2, 2.0)]
            [(x0, 0.0), (x1, 3.0), (x2, 1.0)]
            [(x0, 0.0), (x1, 0.0), (x2, 4.0)]

        """
        import numpy as np

        nb_coordinates = len(coordinates)
        start = np.zeros(nb_coordinates, dtype=int)
        start[0] = derivative_order
        start = tuple(start)

        fathers = {start}
        generation = 0

        while True:

            # Exit if all generation has been done.
            if generation > derivative_order:
                break

            # Yield partial order for derivative [ .. (xi,partial_order_i) .. ] for all i's
            for partial_orders in fathers:
                partial_orders = tuple([int(order) for order in partial_orders])
                yield PartialOrderDerivative(coordinates, partial_orders)

            # Compute new generation
            generation += 1
            sons = set()
            for father in fathers:
                # takes one ball in 0 and distribute it to others
                father = np.asarray(father)
                flux = np.zeros(nb_coordinates, dtype=int)
                flux[0] = -1
                for move in range(1, nb_coordinates):
                    flux[move] = 1
                    son = father + flux
                    sons.add(tuple(son))
                    flux[move] = 0
            fathers = sons


class UnivariateTree(object):
    """ Compute and handle subtitution tree for terms E[D^alpha epsilon D^beta epsilon] 
    of degree |alpha|+|beta|<= max_degree


    Description
    -----------

        Dynamical structure that return the tree at a given degree from stored levels
        It replaces terms in E[D^alpha epsilon D^beta epsilon]  by terms in E[epsilon D^gamma epsilon] 


    """
    def __init__(self, epsilon, spatial_coordinates, self_substitute = True):
        self.epsilon = epsilon
        self.spatial_coordinates = spatial_coordinates
        self._self_substitute = self_substitute
        self._degree = 0
        self._levels = {}
        self._tree = {}

    def __call__(self, degree=None):
        """ Return a tree of degree 'degree'  """

        if degree is None:
            degree = self._degree

        if degree > self._degree:
            # Compute the tree until degree
            for add_degree in range(self._degree, degree):
                self._add_tree_level()

        tree = {}
        for level, level_dico in self._levels.items():
            if level>degree:
                continue
            tree.update(level_dico)
        return tree

    def _add_tree_level(self):

        self._degree += 1
        degree = self._degree

        level = {}

        # Compute all termes E[D^alpha eps D^beta eps] such as |alpha|+|beta| = degree, with 1<=|alpha|<= degree//2
        for alpha_degree in range(1, degree // 2 + 1):
            beta_degree = degree - alpha_degree

            for alpha in PartialOrderDerivative.all_of_degree(self.spatial_coordinates, alpha_degree):

                sub_alpha = alpha.new(alpha.as_sequence[1:])

                for beta in PartialOrderDerivative.all_of_degree(self.spatial_coordinates, beta_degree):

                    # term = Expectation(Derivative(self.epsilon, *alpha) * Derivative(self.epsilon, *beta))
                    term = Derivative(self.epsilon, *alpha.as_couples) * Derivative(self.epsilon, *beta.as_couples)
                    term = Expectation(term)

                    lhs = Derivative(Derivative(self.epsilon, *sub_alpha.as_couples) * Derivative(self.epsilon, *beta.as_couples), alpha.as_sequence[0])
                    rhs = lhs.doit()

                    equation = Eq(Expectation(lhs), Expectation(rhs))

                    equation = equation.isolate(term)

                    # Add into tree
                    subs_term = equation.args[1]
                    if self._self_substitute:
                        subs_term = subs_term.subs(self._tree)
                    self._tree[term] = subs_term
                    level[term] = subs_term

        self._levels[degree] = level

    @property
    def degree(self):
        return self._degree


class Closure(object):

    """
    Exploration of closure from the Taylor expansion of heterogeneous correlation functions,
    written in the aspect tensor approximation (aspect tensor of the correlation 
    approximately equals the local parameter tensor `s` of the correlation)

    \begin{equation}
    \rho(x,y) = \frac{|s_x|^{1/4}|s_y|^{1/4}}{|\frac{1}{2}(s_x+s_y)|^{1/2}}
        \exp\left(-||x-y||^2_{(s_x+s_y)^{-1}}\right)
    \end{equation}
    """

    def __init__(self, mfield):
        if not isinstance(mfield, Field):
            mfield = Field(mfield)
        self.mfield = mfield
        self.spatial_coordinates = mfield.spatial_coordinates

        self._x = self.spatial_coordinates
        self._dx = symbols(' '.join([f"d{xi}" for xi in self.spatial_coordinates]))
        if not isinstance(self._dx, tuple):
            self._dx = (self._dx,)

        self.epsilon = mfield.epsilon
        self.metric = mfield.metric
        self.diffusion = mfield.diffusion
        self.aspect = mfield.aspect
        self._order = 0
        self._taylor = None
        self._in_metric = None
        self._in_diffusion = None
        self._in_aspect = None
        self._closure_by_order = {}

    def correlation_in_diffusion(self, px, py):
        """
        Prototype of correlation function with parameters defined in the so-called 
        diffusion tensor.
        """
        from sympy import exp
        cpx = {xi: pxi for xi, pxi in zip(self.spatial_coordinates, px)}
        cpy = {yi: pyi for yi, pyi in zip(self.spatial_coordinates, py)}
        nux = self.diffusion.subs(cpx)
        nuy = self.diffusion.subs(cpy)
        normalization = (nux.det() ** Rational(1 , 4) * nuy.det() ** Rational(1 , 4)) / (
                    Rational(1 , 2) * (nux + nuy)).det() ** Rational(1 , 2)

        vx = Matrix(px)
        vy = Matrix(py)
        dxy = vx - vy
        dot_prod = dxy.T @ (nux + nuy).inv() @ dxy
        h_correlation = exp(-Rational(1 , 2) * dot_prod[0])
        return normalization * h_correlation

    def correlation(self, px, py):
        """ Heterogeneous Gaussian correlation function in aspect tensor 

        \begin{equation}
        \rho(x,y) = \frac{|s_x|^{1/4}|s_y|^{1/4}}{|\frac{1}{2}(s_x+s_y)|^{1/2}}
            \exp\left(-||x-y||^2_{(s_x+s_y)^{-1}}\right)
        \end{equation}
        
        (see eg Pannekoucke 2021)
        """
        from sympy import exp
        cpx = {xi: pxi for xi, pxi in zip(self.spatial_coordinates, px)}
        cpy = {yi: pyi for yi, pyi in zip(self.spatial_coordinates, py)}
        sx = self.aspect.subs(cpx)
        sy = self.aspect.subs(cpy)
        normalization = (sx.det() ** Rational(1 , 4) * sy.det() ** Rational(1 , 4)) / (
                    Rational(1 , 2) * (sx + sy)).det() ** Rational(1 , 2)

        vx = Matrix(px)
        vy = Matrix(py)
        dxy = vx - vy
        dot_prod = dxy.T @ (sx + sy).inv() @ dxy
        h_correlation = exp(-dot_prod[0])
        return normalization * h_correlation

    def taylor(self, order):
        """ Compute the Taylor expansion until degree 'order' """
        if self._order < order:
            px = self._x
            py = tuple([xi + dxi for xi, dxi in zip(self._x, self._dx)])
            taylor = self.correlation(px, py)
            for dxi in self._dx:
                taylor = taylor.series(dxi, 0, order).removeO()
                taylor = remove_eval_derivative(taylor)

            self._taylor = taylor
            self._order = order

        return self._taylor

    def closure_by_order(self, order):
        """
        Compute the closure for termes E[eps * D^alpha eps]
        :param order:
        :return: a
        """
        order = order + 1
        if self._order < order:
            # 1. Compute terms until order
            factorial = 1
            for partial_order in range(order):

                partial_order = int(partial_order)
                closure = {}

                if partial_order != 0:
                    factorial *= partial_order

                for alpha in PartialOrderDerivative.all_of_degree(self.spatial_coordinates, partial_order):
                    # -1- Extract taylor coefficient of dx^alpha
                    coeff_alpha = self.taylor(order).expand()
                    for xi, dxi in zip(self._x,self._dx):
                        alpha_i = alpha.as_sequence.count(xi)
                        coeff_alpha = coeff_alpha.coeff(dxi, alpha_i)

                    # -2- Extract coefficient of dx^alpha in taylor expansion
                    term = Expectation(self.epsilon * Derivative(self.epsilon, *alpha.as_couples))
                    closure[term] = Integer(factorial) * coeff_alpha
                self._closure_by_order[partial_order] = closure

            # 2. Update the order
            self._order = order

        return self._closure_by_order

    def closure(self, order):
        closure = {}
        for elm in self.closure_by_order(order).values():
            closure.update(elm)

        return closure

    def in_metric_by_order(self, order):
        subs_closure_by_order = self.closure_by_order(order)
        g_dico = {} # denotes the metric tensor
        s_dico = {} # denotes the aspect tensor (here the parameter of the formulation)
        # -1- Supress partial derivative in metric terms
        second_order = 2
        system_s_to_g = []
        for alpha in PartialOrderDerivative.all_of_degree(self.spatial_coordinates, second_order):
            # -1- get term
            term = Expectation(self.epsilon * Derivative(self.epsilon, *alpha.as_couples))
            # -2- Eliminate derivatives
            subs_term = subs_closure_by_order[second_order][term]
            subs_term = subs_term.subs({derivative: 0 for derivative in subs_term.atoms(Derivative)})
            subs_closure_by_order[second_order][term] = subs_term
            # -3-
            if len(self.spatial_coordinates) > 1:
                i, j = alpha.as_tuple[0], alpha.as_tuple[1]
            else:
                i, j = 0, 0

            g_dico[term] = -self.mfield.metric_func(i, j)
            s_dico[term] = self.mfield.aspect_func(i, j)
            system_s_to_g.append(Eq(g_dico[term], subs_term))

        # -2- expression en fonction de la notation de la métrique
        s_to_g = sympy.solve(system_s_to_g, *s_dico.values())

        # -3- Replace in all subs_closure
        for order in subs_closure_by_order:
            for term in subs_closure_by_order[order]:
                subs_closure_by_order[order][term] = subs_closure_by_order[order][term].subs(s_to_g).doit()

        return subs_closure_by_order

    def in_diffusion(self, order):
        raise NotImplementedError

    def in_aspect(self, order):
        raise NotImplementedError


class LocGaussianClosure(Closure):
    """
    Closure associated with the local Gaussian written in aspect tensor.
    """

    def correlation(self, px, py, metric=True):
        from sympy import exp
        cpx = {xi: pxi for xi, pxi in zip(self.spatial_coordinates, px)}
        cpy = {yi: pyi for yi, pyi in zip(self.spatial_coordinates, py)}

        vx = Matrix(px)
        vy = Matrix(py)
        dxy = vx - vy

        if metric:
            gux = self.metric.subs(cpx)
            dot_prod = dxy.T @ gux @ dxy
        else:
            sx = self.aspect.subs(cpx)
            dot_prod = dxy.T @ (Integer(2)*sx).inv() @ dxy

        h_correlation = exp(-dot_prod[0])
        return h_correlation

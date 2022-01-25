from .model_builder import ModelBuilder
from .model_builder_manifold import MetaDerivative as AbstractMetaDerivative

class MetaDerivative(AbstractMetaDerivative):

    @staticmethod
    def slice_coding_rule(coord, step: int):
        return f'self.index({coord},{step})'

    @staticmethod
    def dx_coding_rule(coord):
        return f"self.dx[self.coordinates.index('{coord}')]"

    def _set_as_code(self):
        """ Transform derivative into finite differences """

        from .finite_difference import finite_difference, code_finite_difference

        # -1- Transform derivative into regular finite difference
        finite = finite_difference(self.as_symbolic)

        # -2- Replace functions and dx's
        code = code_finite_difference(finite.simplify(), self.slice_coding_rule, self.dx_coding_rule)

        return code


class FDModelBuilder(ModelBuilder):
    """ Finite Difference Model Builder from symbolic system of partial differential equations

    Create a simple mesh structure where indexes are contains in the dictionary self.index (see template)

    """

    from .templates.body_FDModelBuilder import template as template_FDModelBuilder
    template_body = template_FDModelBuilder

    def _create_meta_derivative(self, derivative):
        return MetaDerivative(derivative)


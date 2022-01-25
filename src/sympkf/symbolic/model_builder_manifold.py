from .model_builder import ModelBuilder
from .model_builder import MetaDerivative as AbstractMetaDerivative

class MetaDerivative(AbstractMetaDerivative):

    def _set_as_code(self):
        return f"self.base_space.diff({self.function.as_keyvar},{self.as_partial_orders})"

class ManifoldModelBuilder(ModelBuilder):
    """ Model Builder from symbolic system of partial differential equations

    This builder is limited to simple functions do not consider mix of 2D/3D, but can handle 2D functions as well as 3D functions

    """

    from .templates.body_ManifoldModelBuilder import template as template_ManifoldModelBuilder
    template_body = template_ManifoldModelBuilder

    def _create_meta_derivative(self, derivative):
        return MetaDerivative(derivative)



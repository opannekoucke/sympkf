template='''
    def __init__(self, base_space, **kwargs):

        super().__init__() # Time scheme is set from Model.__init__()

        # Set base space (manifold)
        #------------------------------
        assert base_space.dimension == len(self.coordinates), \
                                                "Dimension of base_space is different from the number of coordinates"
        self.base_space = base_space
        
        {% if constant_functions_flag %}
        #-----------------------
        # Set constant functions
        #-----------------------
          
        # Set a default nan value for constants
        {% for function in constant_functions %}self.{{function.as_keyvar}} = np.nan # @@ set constant value @@
        {% endfor %}
                
        # Set constant function values from external **kwargs (when provided)
        for key in kwargs:
            if key in self.constant_functions:
                setattr(self, key, kwargs[key])
        
        # Alert when a constant is np.nan
        for function in self.constant_functions:
            if getattr(self, function) is np.nan:
                print(f"Warning: function `{function}` has to be set")
        {% endif %}

        {% if constants_flag %}
        #---------------------------
        # Set constants of the model
        #---------------------------
          
        # Set a default nan value for constants
        {% for constant in constants %}self.{{constant.as_keyvar}} = np.nan # @@ set constant value @@
        {% endfor %}
                
        # Set constant values from external **kwargs (when provided)
        for key in kwargs:
            if key in self.constants:
                setattr(self, key, kwargs[key])
        
        # Alert when a constant is np.nan
        for constant in self.constants:
            if getattr(self, constant) is np.nan:
                print(f"Warning: constant `{constant}` has to be set")
        {% endif %}

'''

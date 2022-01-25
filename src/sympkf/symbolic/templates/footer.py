template = '''
    def trend(self, t, state):
        """ Trend of the dynamics """

        # Init output state with pointer on data
        #-------------------------------------------

        #   a) Set the output array
        dstate = np.zeros(state.shape)

        #   b) Set pointers on output array `dstate` for the computation of the physical trend (alias only).
        {% for function in prognostic_functions -%}
        {{function.trend}} = dstate[{{function.index}}]
        {% endfor %}

        # Load physical functions from state
        #------------------------------------
        {% for function in prognostic_functions -%}
        {{function.as_keyvar}} = state[{{function.index}}]
        {% endfor %}  
        
        {% if constant_functions_flag %}
        # Alias for constant functions
        #-----------------------------
        {% for function in constant_functions -%}
        {{function.as_keyvar}} = self.{{function.as_keyvar}}
        {% endfor %}         
        {% endif %}        
        
        {% if constants_flag -%}
        # Alias for constants
        #--------------------
        {% for constant in constants -%}
        {{constant.as_keyvar}} = self.{{constant.as_keyvar}}
        {% endfor %}         
        {% endif %}        
         
        {% if exogenous_functions_flag %}
        # Compute exogenous functions
        #-----------------------------
        exogenous = self.compute_exogenous(t, state) # None if no exogenous function exists
        {% for function in exogenous_functions -%}
        {{function.as_keyvar}} = exogenous[{{function.index}}]
        {% endfor %} 
        {% endif -%}
        
        # Compute derivative
        #-----------------------
        #
        #  Warning: might be modified to fit appropriate boundary conditions. 
        #
        {% for derivative in spatial_derivatives -%}
        {{derivative.as_keyvar}} = {{derivative.as_code}}
        {% endfor %}
        
        # Implementation of the trend
        #--------------------------------
        {% for line in trend_code -%}
        {{line}}
        {% endfor %}
        
        return dstate        

    {% if exogenous_functions_flag %}
    def compute_exogenous(self, t, state):
        """ Computation of exogenous functions """
        return None # if no exogenous functions are needed
    {% endif %}

'''

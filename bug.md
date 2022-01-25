Known bugs
==========


bug-Expectation
 - tag: 5595950431936087884
 - description: Expectation operator disappear when using Taylor expansion with `.series()` 
 - hotfix: replace the Expectation operator by a temporary function that should be replaced before giving the taylor expansions of the dynamics.
 - fix: refactor the Expectation operator, so it better fit to `sympy` way of code
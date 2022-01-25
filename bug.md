Known bugs
==========


bug-Expectation
 - tag: 5595950431936087884
 - description: Expectation operator disappear when using Taylor expansion with `.series()` 
 - hotfix: replace the Expectation operator by a temporary function that should be replaced before giving the taylor expansions of the dynamics.
 - fix: refactor the Expectation operator, so it better fit to `sympy` way of code

 bug-Expectation-latex
  - tag: 2379311772708329900
  - description: diplay of square of Expectation operator abort under Jupyter
  - hotfix: ?? todo
  - fix: ?? todo
# - initialize to steady state
# - apply perturbation to inputs
# - fix all algebraic variables (requires preciese knowledge of the 
#   algebraic variables) to their steady state values
# - deactivate all algebraic equations
# - solve single finite element (should be a dynamic problem, because inputs reach 
#   differential variables through boundary conditions)
# - algebraic update (?) (based on algebraic equations, which are deactivated) 
# - (re-fix algebraic variables... does setting value unfix variables?)
# - unfix/activate "first wave" of algebraic variables/equations
#   (could be as little as a single constraint/variable)
#   (may or may not change trajectory of differential variables)


# what happens when you have a problem, solved at a point, then
# add a constraint/variable (add a dimension), and try to solve for that 
# varaible?
# (what properties of the constraint set, projected into the dimension of the new 
# variable, guarantee convergence ??? (if any... - convexity, likely) )

# implicit euler step: algebraic variable at t0 is free (alg_update)
# algebraic variable at t1 is the only degree of freedom that has been added,
# but its initialization is incorrect, and when it gets recalculated to 
# satisfy its algebraic equation, the differential equations are no longer 
# satisfied 

# why do I expect this problem to be any easier to solve than just solving the
# the problem all at once ?????
# ^ no clear reason, yet, but it's at least something to try

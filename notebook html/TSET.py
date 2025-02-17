import sympy as sp

# Define symbolic variables
Tp, Tn, Tnp, Tpp, sigma_n, sigma_p, P = sp.symbols('Tp Tn Tnp Tpp sigma_n sigma_p P')

# Coefficient matrix (10x6)
A = sp.Matrix([
    [(Tp**48 / Tp**40) * (1 + (sigma_n / sigma_p) * P), 
     2 * (Tp**48 / Tp**40), 
     2 * (sigma_n / sigma_p) * (Tp**48 / Tp**40) * P, 
     -((sigma_n / sigma_p) + P), 
     -2, 
     -2 * (sigma_n / sigma_p) * P],
    
    [(Tn**48 / Tn**40) * ((sigma_n / sigma_p) + P), 
     2 * (Tn**48 / Tn**40) * P, 
     2 * (sigma_n / sigma_p) * (Tn**48 / Tn**40), 
     -((sigma_n / sigma_p) + P), 
     -2 * P, 
     -2 * (sigma_n / sigma_p)],

    [(sigma_n / sigma_p) + 1, 
     2, 
     2 * (sigma_n / sigma_p), 
     -((sigma_n / sigma_p) + 1), 
     -2, 
     -2 * (sigma_n / sigma_p)],

    [(sigma_n / sigma_p) * (Tnp**48 / Tnp**40), 
     2 * P * (Tnp**48 / Tnp**40), 
     2 * (sigma_n / sigma_p) * (Tnp**48 / Tnp**40) * P, 
     -(sigma_n / sigma_p), 
     -2 * P, 
     -2 * (sigma_n / sigma_p) * P],

    [(sigma_n / sigma_p) * (Tpp**48 / Tpp**40) * P, 
     (Tpp**48 / Tpp**40), 
     0, 
     -(sigma_n / sigma_p) * P, 
     -1, 
     0],

    [(sigma_n / sigma_p) * (Tpp**48 / Tp**48) * P - (sigma_n / sigma_p) * P - 1, 
     (Tpp**48 / Tp**48) - 2, 
     -2 * (sigma_n / sigma_p) * P, 
     0, 
     0, 
     0],

    [(sigma_n / sigma_p) * (Tnp**48 / Tn**48) - (sigma_n / sigma_p) - P, 
     2 * P * ((Tnp**48 / Tn**48) - 1), 
     2 * (sigma_n / sigma_p) * ((Tnp**48 / Tn**48) * P - 1), 
     0, 
     0, 
     0],

    [0, 0, 0, 
     (sigma_n / sigma_p) * (Tpp**40 / Tp**40) * P - (sigma_n / sigma_p) * P - 1, 
     (Tpp**40 / Tp**40) - 2, 
     -2 * (sigma_n / sigma_p) * P],

    [0, 0, 0, 
     (sigma_n / sigma_p) * ((Tnp**40 / Tn**40) - P - 1), 
     2 * P * ((Tnp**40 / Tn**40) - 1), 
     2 * (sigma_n / sigma_p) * ((Tnp**40 / Tn**40) * P - 1)],

    [0, 0, 0, 0, 1, -1]
])

# Define variables x1, x2, ..., x6 (6 variables to fit the 10x6 system)
x1, x2, x3, x4, x5, x6 = sp.symbols('x1 x2 x3 x4 x5 x6')

# Construct the right-hand side vector (zero vector for homogenous system)
b = sp.Matrix([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

# Solve the system A * x = b for variables [x1, x2, x3, x4, x5, x6]
solution = sp.linsolve((A, b), x1, x2, x3, x4, x5, x6)
print(solution)

#!/usr/bin/env python3

import sys
import re
from typing import Dict, Tuple

def parse_term(term: str) -> Tuple[float, int]:
    """Parse a term like '5 * X^2' and return (coefficient, power)"""
    term = term.strip()
    
    # Match pattern: coefficient * X^power or just coefficient * X^0
    match = re.match(r'([+-]?\s*[\d.]+)\s*\*\s*X\s*\^\s*(\d+)', term)
    if match:
        coef = float(match.group(1).replace(' ', ''))
        power = int(match.group(2))
        return (coef, power)
    
    raise ValueError(f"Invalid term format: {term}")

def parse_equation(equation: str) -> Dict[int, float]:
    """Parse the equation and return coefficients as {power: coefficient}"""
    # Split by '='
    if '=' not in equation:
        raise ValueError("Equation must contain '='")
    
    left, right = equation.split('=')
    
    # Split terms by + or - (keeping the sign)
    def split_terms(expr: str):
        # Add '+' at the beginning if the expression doesn't start with a sign
        expr = expr.strip()
        if expr and expr[0] not in ['+', '-']:
            expr = '+' + expr
        
        terms = re.split(r'(?=[+-])', expr)
        return [t.strip() for t in terms if t.strip()]
    
    coefficients = {}
    
    # Process left side
    for term in split_terms(left):
        if term:
            coef, power = parse_term(term)
            coefficients[power] = coefficients.get(power, 0) + coef
    
    # Process right side (subtract from left)
    for term in split_terms(right):
        if term:
            coef, power = parse_term(term)
            coefficients[power] = coefficients.get(power, 0) - coef
    
    return coefficients

def get_degree(coefficients: Dict[int, float]) -> int:
    """Get the degree of the polynomial (highest non-zero power)"""
    degree = 0
    for power, coef in coefficients.items():
        if abs(coef) > 1e-10 and power > degree:
            degree = power
    return degree

def print_reduced_form(coefficients: Dict[int, float]):
    """Print the equation in reduced form"""
    max_power = max(coefficients.keys()) if coefficients else 0
    terms = []
    
    for power in range(max_power + 1):
        coef = coefficients.get(power, 0)
        if power == 0 or abs(coef) > 1e-10 or not terms:
            sign = '+' if coef >= 0 else '-'
            abs_coef = abs(coef)
            
            if terms:  # Not the first term
                terms.append(f"{sign} {abs_coef} * X^{power}")
            else:  # First term
                if coef >= 0:
                    terms.append(f"{abs_coef} * X^{power}")
                else:
                    terms.append(f"-{abs_coef} * X^{power}")
    
    if not terms:
        print("Reduced form: 0 * X^0 = 0")
    else:
        print(f"Reduced form: {' '.join(terms)} = 0")

def sqrt(n: float) -> float:
    """Calculate square root using Newton's method"""
    if n < 0:
        raise ValueError("Cannot compute square root of negative number")
    if n == 0:
        return 0
    
    x = n
    while True:
        root = 0.5 * (x + n / x)
        if abs(root - x) < 1e-10:
            return root
        x = root

def solve_degree_0(coefficients: Dict[int, float]):
    """Solve equation of degree 0"""
    c = coefficients.get(0, 0)
    
    if abs(c) < 1e-10:
        print("Any real number is a solution.")
    else:
        print("No solution.")

def solve_degree_1(coefficients: Dict[int, float]):
    """Solve equation of degree 1: ax + b = 0"""
    a = coefficients.get(1, 0)
    b = coefficients.get(0, 0)
    
    solution = -b / a
    print("The solution is:")
    print(solution)

def solve_degree_2(coefficients: Dict[int, float]):
    """Solve equation of degree 2: ax^2 + bx + c = 0"""
    a = coefficients.get(2, 0)
    b = coefficients.get(1, 0)
    c = coefficients.get(0, 0)
    
    # Calculate discriminant
    discriminant = b * b - 4 * a * c
    
    if discriminant > 1e-10:
        print("Discriminant is strictly positive, the two solutions are:")
        sqrt_disc = sqrt(discriminant)
        sol1 = (-b + sqrt_disc) / (2 * a)
        sol2 = (-b - sqrt_disc) / (2 * a)
        print(sol1)
        print(sol2)
    elif abs(discriminant) <= 1e-10:
        print("Discriminant is zero, the solution is:")
        sol = -b / (2 * a)
        print(sol)
    else:
        print("Discriminant is strictly negative, the two complex solutions are:")
        real_part = -b / (2 * a)
        imag_part = sqrt(-discriminant) / (2 * a)
        
        # Format as fractions if they simplify nicely
        print(f"{real_part} + {imag_part}i")
        print(f"{real_part} - {imag_part}i")

def solve_equation(equation: str):
    """Main function to solve the equation"""
    try:
        coefficients = parse_equation(equation)
        print_reduced_form(coefficients)
        
        degree = get_degree(coefficients)
        print(f"Polynomial degree: {degree}")
        
        if degree > 2:
            print("The polynomial degree is strictly greater than 2, I can't solve.")
        elif degree == 2:
            solve_degree_2(coefficients)
        elif degree == 1:
            solve_degree_1(coefficients)
        else:  # degree 0
            solve_degree_0(coefficients)
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    if len(sys.argv) > 1:
        equation = sys.argv[1]
    else:
        print("Enter equation:")
        equation = input()
    
    solve_equation(equation)

if __name__ == "__main__":
    main()
#!/usr/bin/env python3

def parse_equation(equation_string):
    """
    Parse the equation into left and right terms.
    Returns two lists of terms.
    """
    left_side, right_side = equation_string.split("=")
    return left_side.strip().split(), right_side.strip().split()


def extract_coefficients(terms):
    """
    Convert a list of terms into a dictionary of coefficients {exponent: coefficient}.
    """
    coefficients = {}
    sign = 1
    i = 0

    while i < len(terms):
        token = terms[i]

        if token == "+":
            sign = 1
            i += 1
            continue
        elif token == "-":
            sign = -1
            i += 1
            continue

        # Example of a term: "5", "*", "X^2"
        coefficient = float(token) * sign
        exponent = int(terms[i + 2].split("^")[1])
        coefficients[exponent] = coefficients.get(exponent, 0.0) + coefficient

        i += 3
    return coefficients


def reduce_equation(left_terms, right_terms):
    """
    Move all terms to the left-hand side and return reduced coefficients.
    """
    left_coeffs = extract_coefficients(left_terms)
    right_coeffs = extract_coefficients(right_terms)

    for exp, coeff in right_coeffs.items():
        left_coeffs[exp] = left_coeffs.get(exp, 0.0) - coeff

    return left_coeffs


def get_polynomial_degree(coefficients):
    """
    Return the degree of the polynomial.
    """
    non_zero_exps = [exp for exp, coeff in coefficients.items() if coeff != 0]
    if not non_zero_exps:
        return 0
    return max(non_zero_exps)


def format_reduced_form(coefficients):
    """
    Return a string for the reduced form.
    """
    if not coefficients:
        return "0 * X^0 = 0"

    terms = []
    for exp in sorted(coefficients.keys()):
        coeff = coefficients[exp]
        terms.append(f"{coeff} * X^{exp}")
    return " + ".join(terms).replace("+ -", "- ") + " = 0"


def solve_polynomial(coefficients):
    """
    Solve polynomial equations of degree 0, 1, or 2.
    """
    degree = get_polynomial_degree(coefficients)
    print(f"Polynomial degree: {degree}")

    if degree > 2:
        print("The polynomial degree is strictly greater than 2, I can't solve.")
        return

    a = coefficients.get(2, 0.0)
    b = coefficients.get(1, 0.0)
    c = coefficients.get(0, 0.0)

    if degree == 0:
        if c == 0:
            print("Any real number is a solution.")
        else:
            print("No solution.")
    elif degree == 1:
        solution = -c / b
        print("The solution is:")
        print(solution)
    elif degree == 2:
        discriminant = b ** 2 - 4 * a * c
        if discriminant > 0:
            root1 = (-b + discriminant ** 0.5) / (2 * a)
            root2 = (-b - discriminant ** 0.5) / (2 * a)
            print("Discriminant is strictly positive, the two solutions are:")
            print(root1)
            print(root2)
        elif discriminant == 0:
            root = -b / (2 * a)
            print("Discriminant is zero, the solution is:")
            print(root)
        else:
            real_part = -b / (2 * a)
            imaginary_part = (abs(discriminant) ** 0.5) / (2 * a)
            print("Discriminant is strictly negative, the two complex solutions are:")
            print(f"{real_part} + {imaginary_part}i")
            print(f"{real_part} - {imaginary_part}i")


def main():
    equation_string = input("Enter the polynomial equation: ").strip()
    left_terms, right_terms = parse_equation(equation_string)
    coefficients = reduce_equation(left_terms, right_terms)
    reduced_form = format_reduced_form(coefficients)
    print("Reduced form:", reduced_form)
    solve_polynomial(coefficients)


if __name__ == "__main__":
    main()

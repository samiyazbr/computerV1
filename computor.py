#!/usr/bin/env python3

"""
Computor v1: Solve polynomial equations up to degree 2.

Features
- Parse input equation strings like: "5 * X^0 + 4 * X^1 - 9.3 * X^2 = 1 * X^0"
- Reduce to canonical form on the left side: sum(coeff * X^degree) = 0
- Display polynomial degree
- Solve for degree 0, 1, and 2; classify discriminant and print solutions
- Handle zero, negative, and non-integer coefficients

Optional robustness
- Accept free-form inputs such as "5 + 4 * X + X^2 = X^2" (implied coefficients and exponents)
- Basic malformed-input handling with helpful error messages

This program avoids external libraries and aims for clear, readable code for extension.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple, List
import sys

# Ensure Python 3.10 compatibility
if sys.version_info < (3, 10):
    print("This program requires Python 3.10 or higher.")
    sys.exit(1)


@dataclass
class Polynomial:
    """Represents a polynomial as a mapping degree -> coefficient (floats)."""

    coefficients: Dict[int, float]

    def degree(self) -> int:
        non_zero_degrees = [deg for deg, coef in self.coefficients.items() if coef != 0.0]
        if not non_zero_degrees:
            return 0
        return max(non_zero_degrees)

    def normalized(self) -> "Polynomial":
        """Remove zero coefficients to keep the representation tidy."""
        return Polynomial({deg: coef for deg, coef in self.coefficients.items() if coef != 0.0})

    def __add__(self, other: "Polynomial") -> "Polynomial":
        result: Dict[int, float] = dict(self.coefficients)
        for deg, coef in other.coefficients.items():
            result[deg] = result.get(deg, 0.0) + coef
        return Polynomial(result).normalized()

    def __sub__(self, other: "Polynomial") -> "Polynomial":
        result: Dict[int, float] = dict(self.coefficients)
        for deg, coef in other.coefficients.items():
            result[deg] = result.get(deg, 0.0) - coef
        return Polynomial(result).normalized()


def parse_number(text: str) -> float:
    """Parse a numeric string into a float. Accept integers and decimals."""
    text = text.strip()
    if not text:
        raise ValueError("Empty number")
    try:
        return float(text)
    except Exception as exc:
        raise ValueError(f"Invalid number: {text}") from exc


def tokenize_side(side: str) -> List[str]:
    """Split a polynomial side into additive term strings preserving signs.

    Example: "- 3 * X^2 + 4X - X + 5" -> ["- 3 * X^2", "+ 4X", "- X", "+ 5"]
    """
    # Normalize spaces around operators, keep '-' as a sign. Convert all minus to '+ -' then split on '+'.
    s = side.replace("-", "+ -")
    parts = s.split("+")
    terms: List[str] = []
    for part in parts:
        p = part.strip()
        if not p:
            continue
        # Prepend '+' for positive terms to normalize downstream formatting if missing
        if not (p.startswith("+") or p.startswith("-")):
            p = "+ " + p
        terms.append(p)
    return terms


def parse_term(term: str) -> Tuple[int, float]:
    """Parse a single term into (degree, coefficient).

    Supported forms (case-insensitive X):
    - "+ 3 * X^2", "- 3 * X^2", "+ X^2", "- X^2"
    - "+ 4 * X", "+ 4X", "- X"
    - "+ 5" (constant)
    The '*' between coefficient and X is optional. The '^exp' is optional (default 1 for X present).
    """
    t = term.strip().replace(" ", "")
    if not t:
        raise ValueError("Empty term")

    # Normalize sign
    sign = 1
    if t[0] == '+':
        t = t[1:]
    elif t[0] == '-':
        sign = -1
        t = t[1:]

    if not t:
        raise ValueError("Dangling sign without a term")

    # Case-insensitive X
    t = t.replace('x', 'X')

    # If no X, it's a constant number
    if 'X' not in t:
        coef = parse_number(t) * sign
        return (0, coef)

    # Strict required format for variable terms: "<number>*X^<int>"
    # Optional whitespace already removed; '*' and '^' are required for non-constant terms
    # Find '*'
    if '*' not in t:
        raise ValueError(f"Invalid term (missing '*'): {term}")
    coef_str, after_star = t.split('*', 1)
    if not after_star.startswith('X'):
        raise ValueError(f"Invalid term (missing 'X'): {term}")
    after_x = after_star[1:]
    if not after_x.startswith('^'):
        raise ValueError(f"Invalid term (missing '^'): {term}")
    exp_str = after_x[1:]
    if exp_str == "":
        raise ValueError(f"Missing exponent after '^' in term: {term}")
    try:
        degree = int(exp_str)
    except Exception as exc:
        raise ValueError(f"Non-integer exponent in term: {term}") from exc
    if degree < 0:
        raise ValueError(f"Negative exponents are not supported: {term}")

    coef = parse_number(coef_str) * sign

    return (degree, coef)


def parse_side_to_poly(side: str) -> Polynomial:
    """Parse a side of an equation into a Polynomial."""
    terms = tokenize_side(side)
    coefficients: Dict[int, float] = {}
    for raw_term in terms:
        degree, coef = parse_term(raw_term)
        coefficients[degree] = coefficients.get(degree, 0.0) + coef
    return Polynomial(coefficients).normalized()


def parse_equation(equation: str) -> Tuple[Polynomial, Polynomial]:
    """Parse a full equation string into (left_poly, right_poly)."""
    if '=' not in equation:
        raise ValueError("Equation must contain '=' sign")
    left_str, right_str = equation.split('=', 1)
    if left_str.strip() == "" or right_str.strip() == "":
        raise ValueError("Both sides of the equation must be non-empty")
    left_poly = parse_side_to_poly(left_str)
    right_poly = parse_side_to_poly(right_str)
    return left_poly, right_poly


def reduce_equation(left: Polynomial, right: Polynomial) -> Polynomial:
    """Move all terms to the left-hand side and return the reduced polynomial."""
    return (left - right).normalized()


def format_coefficient(value: float) -> str:
    """Format a coefficient for display in reduced form.

    - Print integers without trailing .0
    - Otherwise print as decimal or fraction depending on simplicity
    """
    # Print integers without trailing .0
    if value.is_integer():
        return str(int(value))
    # Keep up to 10 decimals and trim
    text = f"{value:.10f}".rstrip('0').rstrip('.')
    return text if text else "0"


def format_reduced_form(poly: Polynomial) -> str:
    """Return the reduced form string like: "4 * X^0 + 4 * X^1 - 9.3 * X^2 = 0"."""
    if not poly.coefficients:
        return "0 = 0"
    parts: List[str] = []
    for degree in sorted(poly.coefficients.keys()):
        coef = poly.coefficients[degree]
        if coef == 0:
            continue
        sign = "-" if coef < 0 else "+"
        abs_coef = -coef if coef < 0 else coef
        coef_str = format_coefficient(abs_coef)
        term = f"{coef_str} * X^{degree}"
        parts.append((sign, term))

    if not parts:
        return "0 = 0"

    # Build string with correct spacing and no leading '+'
    first_sign, first_term = parts[0]
    output = first_term if first_sign == "+" else f"- {first_term}"
    for sign, term in parts[1:]:
        if sign == "+":
            output += f" + {term}"
        else:
            output += f" - {term}"

    return f"{output} = 0"


def solve_degree_0(c: float) -> str:
    if c == 0:
        return "All real numbers are solutions."
    return "No solution."


def solve_degree_1(a: float, b: float) -> Tuple[str, List[str]]:
    # a * X + b = 0  => X = -b / a
    if a == 0:
        return ("degenerate", [solve_degree_0(b)])
    x = -b / a
    return ("linear", [format_float_solution(x)])


def sqrt_float(x: float) -> float:
    # Basic square root using exponentiation; input assumed non-negative
    return x ** 0.5


def format_float_solution(x: float) -> str:
    return f"{x:.6f}"


def solve_degree_2(a: float, b: float, c: float) -> Tuple[str, List[str]]:
    if a == 0:
        # Fallback to linear
        kind, sols = solve_degree_1(b, c)
        return ("linear_fallback", sols)

    # Compute discriminant
    D = b ** 2 - 4.0 * a * c
    if D > 0:
        sqrtD = sqrt_float(D)
        x1 = (-b + sqrtD) / (2.0 * a)
        x2 = (-b - sqrtD) / (2.0 * a)
        return ("positive", [format_float_solution(x1), format_float_solution(x2)])
    elif D == 0:
        x = (-b) / (2.0 * a)
        return ("zero", [format_float_solution(x)])
    else:
        # Complex solutions
        sqrt_abs = sqrt_float(-D)
        real = (-b) / (2.0 * a)
        imag = sqrt_abs / (2.0 * a)
        # Format as a + bi and a - bi with 6 decimals
        real_s = format_float_solution(real)
        imag_s = format_float_solution(abs(imag))
        s1 = f"{real_s} + {imag_s}i"
        s2 = f"{real_s} - {imag_s}i"
        return ("negative", [s1, s2])


def solve_equation(poly: Polynomial) -> None:
    deg = poly.degree()
    print(f"Polynomial degree: {deg}")
    # Extract coefficients a2, a1, a0
    a2 = poly.coefficients.get(2, 0.0)
    a1 = poly.coefficients.get(1, 0.0)
    a0 = poly.coefficients.get(0, 0.0)

    if deg == 0:
        print(solve_degree_0(a0))
        return
    if deg == 1:
        kind, sols = solve_degree_1(a1, a0)
        if kind == "degenerate":
            print(sols[0])
        else:
            print("The solution is:")
            print(sols[0])
        return
    if deg == 2:
        kind, sols = solve_degree_2(a2, a1, a0)
        if kind == "positive":
            print("Discriminant is strictly positive, the two solutions are:")
            print(sols[0])
            print(sols[1])
        elif kind == "zero":
            print("Discriminant is zero, the solution is:")
            print(sols[0])
        elif kind == "negative":
            print("Discriminant is strictly negative, the two complex solutions are:")
            print(sols[0])
            print(sols[1])
        else:
            # linear fallback
            print("Degenerated to a linear equation. The solution is:")
            print(sols[0])
        return

    print("The polynomial degree is strictly greater than 2, I can't solve.")


def display_results(reduced: Polynomial) -> None:
    reduced_str = format_reduced_form(reduced)
    print(f"Reduced form: {reduced_str}")
    solve_equation(reduced)


def main(argv: List[str]) -> int:
    """Entry point. Accept equation string either as a single argument or via stdin."""
    if len(argv) >= 2:
        # Accept the entire expression as a single argument, or join all arguments
        equation = " ".join(argv[1:])
    else:
        try:
            equation = input("Enter equation: ").strip()
        except EOFError:
            print("No input provided.")
            return 1

    if not equation:
        print("No input provided.")
        return 1

    try:
        left, right = parse_equation(equation)
        reduced = reduce_equation(left, right)
        display_results(reduced)
    except ValueError as err:
        print(f"Error: {err}")
        return 2
    except Exception as err:
        print(f"Unexpected error: {err}")
        return 3
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))



import sympy as sp
from copy import deepcopy
from numbers import Number


number_types = (Number, sp.Pow, sp.Symbol)


class State:

    def __init__(self, numbers, amplitudes=1):

        if not isinstance(numbers, list):
            self.numbers = [[numbers]]
            self.amplitudes = [amplitudes]

        elif not isinstance(numbers[0], list):
            self.numbers = [numbers]
            self.amplitudes = [amplitudes]

        else:
            self.numbers = numbers
            self.amplitudes = amplitudes

    def __add__(self, other):

        if not isinstance(other, State):
            raise TypeError('Invalid argument for addition.')

        new_numbers = deepcopy(self.numbers)
        new_amplitudes = deepcopy(self.amplitudes)

        for number, amplitude in zip(other.numbers, other.amplitudes):
            if number in new_numbers:
                i = new_numbers.index(number)
                new_amplitudes[i] += amplitude
            else:
                new_numbers.append(number)
                new_amplitudes.append(amplitude)

        return State(new_numbers, new_amplitudes)

    def __sub__(self, other):

        if not isinstance(other, State):
            raise TypeError('Invalid argument for addition.')

        return self.__add__((-1)*other)

    def __mul__(self, other):

        if isinstance(other, number_types):
            return self.__rmul__(other)
        else:
            raise TypeError('Invalid argument for multiplication.')

    def __rmul__(self, other):

        if isinstance(other, number_types):
            new_amplitudes = [other*a for a in self.amplitudes]
            return State(self.numbers, new_amplitudes)

        elif isinstance(other, Operator):
            return other.operate(self)

        else:
            raise TypeError('Invalid argument for multiplication.')

    def create(self, i=0):

        new_numbers = deepcopy(self.numbers)
        new_amplitudes = deepcopy(self.amplitudes)

        for k, n in enumerate(new_numbers):
            new_amplitudes[k] *= sp.sqrt(n[i]+1)
            new_numbers[k][i] += 1

        return State(new_numbers, new_amplitudes)

    def annihilate(self, i=0):

        new_numbers = deepcopy(self.numbers)
        new_amplitudes = deepcopy(self.amplitudes)

        for k, n in enumerate(new_numbers):
            new_amplitudes[k] *= sp.sqrt(n[i])
            new_numbers[k][i] -= 1

        return State(new_numbers, new_amplitudes)

    def norm(self):
        return sp.sqrt(scalar_product(self, self))

    def show(self):

        string = ''

        for num, amp in zip(self.numbers, self.amplitudes):

            if amp == 0:
                pass

            else:
                if amp == 1:
                    str_a = ''
                else:
                    str_a = str(sp.simplify(amp))

                str_n = ''
                for n in num:
                    str_n += str(n) + ', '

                string += str_a + '|' + str_n[:-2] + '> + '

        print(string[:-2])


class Operator:

    def __init__(self, terms, coeffs=1):

        if not isinstance(terms, list):
            self.terms = [[terms]]
            self.coeffs = [coeffs]

        elif not isinstance(terms[0], list):
            self.terms = [terms]
            self.coeffs = [coeffs]

        else:
            self.terms = terms
            self.coeffs = coeffs

    def __add__(self, other):

        if not isinstance(other, Operator):
            raise TypeError('Invalid argument for addition.')

        new_terms = deepcopy(self.terms)
        new_coeffs = deepcopy(self.coeffs)

        for term, coeff in zip(other.terms, other.coeffs):
            if term in new_terms:
                i = new_terms.index(term)
                new_coeffs[i] += coeff
            else:
                new_terms.append(term)
                new_coeffs.append(coeff)

        return Operator(new_terms, new_coeffs)

    def __sub__(self, other):

        if not isinstance(other, Operator):
            raise TypeError('Invalid argument for addition.')

        new_coeffs = [-c for c in other.coeffs]
        return self + Operator(other.terms, new_coeffs)

    def __mul__(self, other):

        if isinstance(other, number_types):

            new_coeffs = [other*c for c in self.coeffs]
            return Operator(self.terms, new_coeffs)

        elif isinstance(other, Operator):

            new_terms, new_coeffs = [], []

            for t1, c1 in zip(self.terms, self.coeffs):
                for t2, c2 in zip(other.terms, other.coeffs):
                    new_terms.append(t1 + t2)
                    new_coeffs.append(c1*c2)

            return Operator(new_terms, new_coeffs)

        elif isinstance(other, State):
            return other.__rmul__(self)

        else:
            raise TypeError('Invalid argument for multiplication.')

    def __rmul__(self, other):

        if isinstance(other, number_types):
            return self.__mul__(other)

        elif isinstance(other, Operator):
            return other.__mul__(self)

        else:
            raise TypeError('Invalid argument for multiplication.')

    def __pow__(self, exp):

        if not isinstance(exp, int) or exp < 1:
            raise TypeError('Exponent must be a positive integer.')

        result = deepcopy(self)
        for i in range(exp-1):
            result *= deepcopy(self)

        return result

    def operate(self, other):

        state = self.coeffs[0]*deepcopy(other)

        for operator in reversed(self.terms[0]):

            if operator[0] == '+':
                state = state.create(operator[1])
            elif operator[0] == '-':
                state = state.annihilate(operator[1])
            else:
                pass

        for i, term in enumerate(self.terms[1:]):

            state_add = deepcopy(other)

            for operator in reversed(term):

                if operator[0] == '+':
                    state_add = state_add.create(operator[1])
                elif operator[0] == '-':
                    state_add = state_add.annihilate(operator[1])
                else:
                    pass

            state += self.coeffs[i+1]*state_add

        return state

    def show(self):

        string = ''

        for coeff, term in zip(self.coeffs, self.terms):

            if coeff == 0:
                pass

            else:
                if coeff == 1:
                    str_c = ''
                else:
                    str_c = str(sp.simplify(coeff))

                str_t = ''
                for factor in term:
                    str_t += str(factor[0]) + str(factor[1]) + ' '

                string += str_c + '(' + str_t[:-1] + ') + '

        print(string[:-2])


def scalar_product(state_1, state_2):

    def basis_scalar_product(numbers_1, numbers_2):

        if len(numbers_1) != len(numbers_2):
            return 0

        else:
            s = 1
            for n_1, n_2 in zip(numbers_1, numbers_2):
                s *= sp.KroneckerDelta(n_1, n_2)
            return s

    s = 0
    for i, n1 in enumerate(state_1.numbers):
        for k, n2 in enumerate(state_2.numbers):
            s += sp.conjugate(state_1.amplitudes[i])*state_2.amplitudes[k] \
                    * basis_scalar_product(n1, n2)

    return sp.simplify(s)


if __name__ == '__main__':

    m = sp.Symbol('m', real=True, positive=True)
    n = sp.Symbol('n', real=True, positive=True)
    p = sp.Symbol('p', real=True, positive=True)
    q = sp.Symbol('q', real=True, positive=True)

    psi = State([m, n])
    phi = State([p, q])
    psi.show()      # |psi> = |m, n>
    phi.show()      # |phi> = |p, q>

    a = Operator(('-', 0))      # Lowering operator of the first particle
    A = Operator(('+', 0))      # Raising operator of the first particle
    b = Operator(('-', 1))      # Lowering operator of the second particle
    B = Operator(('+', 1))      # Raising operator of the second particle

    # States can be added to each other, and multiplied (from the left) by
    # numbers and operators

    chi = 3*psi + a*phi
    chi.show()

    # Operators can be added to each other, multiplied by numbers and other
    # operators, and raised to positive integer powers

    C = a - 2*a*(b + B)
    chi = a*psi + b*chi
    C.show()
    chi.show()
    (C*chi).show()

    D = a*B + (a + b*A*B)**2
    D.show()
    chi = D*psi
    chi.show()

    # Norm of a state and scalar product between two states are defined

    print(psi.norm())
    print((((A*a + B*b))*psi).norm())
    print(scalar_product(psi, phi))
    print(scalar_product(phi, chi))

    # A standard exercise in quantum mechanics: Compute the expectation value
    # of the operator X = (a + a^\dagger)^4 in the state |n>

    X = (a + A)**4
    print(scalar_product(psi, X*psi))

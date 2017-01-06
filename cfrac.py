#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
# Name:       cfrac
# Created:    14.10.16
# Author:     Carlos Esparza Sanchez
#-------------------------------------------------------------------------------
import numbers

import itertools as it
import collections
import math

from typing import List, Iterable
from fractions import Fraction
from gmpy2 import mpz



def euclid_factors(a: numbers.Real, b: numbers.Real):
    """
    An iterator which yield the remainders of the euclidean algorithm applied
    to a and b.
    """
    while b != 0:
        q = math.floor(a / b)
        yield int(q)
        a, b = b, a - q*b


class CFrac(numbers.Real):
    """A continued fraction.

    CFrac(7)                  -> <7>
    CFrac(3.14)               -> <3, 2, 1, 2>
    CFrac(Fraction(123, 456)) -> <0, 3, 1, 2, 2, 2, 2>

        if called with a numeric type it returns the corresponding CF

    CFrac((1 for _ in range(100))) -> <1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ... >

        if called with an iterable it return a CF with the elements of the
        iterable as coefficients

    """

    # determines how many coefficients are compared by ==, <, >, etc.
    DEPTH = 25

    def __init__(self, x):
        """
        :param x: a number or an iterable yielding the CF terms
        """
        self._cached = -1
        self._terms = []

        if isinstance(x, numbers.Real):
            self._gen = euclid_factors(x, 1)

        elif isinstance(x, collections.Iterable):
            self._gen = iter(x)
            try:
                self._more_terms()
            except StopIteration:
                raise ValueError('iterable cannot be empty')

        elif isinstance(x, numbers.Complex):
            raise TypeError("complex numbers are not supported")
        else:
            raise TypeError('argument must be a real number or an iterable')


    def __getitem__(self, index):
        """
        cfrac[n] -> int

            returns the n-th coefficient of the CF

        cfrac[from:to] -> CFrac

            returns the CF with the corresponding coefficients
        """
        if isinstance(index, int):
            try:
                return self._terms[index]
            except IndexError: pass
            try:
                self._more_terms(index - self._cached)
            except StopIteration:
                # We have extracted all coefficients from the iterator. We
                # delete our reference to the iterator so it can be garbage-
                # collected
                del self._gen
                self._more_terms = self._stopper
                raise IndexError('continued fraction not that long')
            return self._terms[index]

        elif isinstance(index, slice):
            if index.start is None:
                return CFrac(it.islice(self, index.stop))
            elif self.longer_than_eq(index.start + 1):
                return CFrac(it.islice(self, index.start, index.stop, index.step))
            else:
                return iter([])

        else:
            raise TypeError('indices must be integers or slices')


    @staticmethod
    def _stopper(self, n=1):
        raise IndexError('continued fraction not that long')

    def _more_terms(self, n=1):
        for _ in range(n):
            # next(self._gen) should already be integral (e. g. 2.0)
            term = int( next(self._gen) )
            self._terms.append(term)
        self._cached += n

    def longer_than_eq(self, n):
        try:
            self[n - 1]
            return True
        except IndexError:
            return False

    def __float__(self):
        A, A_ = self[0], 1
        B, B_ = 1, 0
        for coeff in self[1:self.DEPTH]:
            A, A_ = coeff * A + A_, A
            B, B_ = coeff * B + B_, B
        return A / B

    def __bool__(self):
        return any(x != 0 for x in self[:self.DEPTH])

    def __abs__(self):
        if self < 0: return -self
        else: return self


    def __trunc__(self): # round towards zero
        if self[0] >= 0: return self[0]
        else: return self[0] + 1

    def __round__(self, n=None):
        if n < 0:
            return CFrac(round(self[0], n))
        if n == 1:
            return self[0]
        else:
            return self[:n + 1]

    def __floor__(self):
        return self[0]

    def __ceil__(self):
        if self.longer_than_eq(2):
            return self[0] + 1
        else:
            return self[0]

    def __repr__(self):
        """
        conversion to str
        """
        dots = ', ... ' if self.longer_than_eq(11) else ''
        return '<' + ', '.join(str(x) for x in self[:10]) + dots + '> '


    def __eq__(self, other):
        if isinstance(other, CFrac):
            return all(x == y for x, y in zip(self[:self.DEPTH], other[:self.DEPTH]))
        elif isinstance(object, numbers.Real):
            return bool(self - other)

        else:
            return NotImplemented

    def __lt__(self, other):
        if isinstance(other, CFrac):
            for i in range(self.DEPTH):
                next_self = self[i] if self.longer_than_eq(i + 1) else math.inf
                next_other = other[i] if other.longer_than_eq(i + 1) else math.inf

                if next_self != next_other:
                    return bool((next_self < next_other) ^ (i % 2)) # XOR
                elif next_self == math.inf:
                    return False

            return False # keine unterschied in den ersten DEPTH Gliedern
        else:
            return self - other < CFrac(0)

    def __le__(self, other):
        return self < other or self == other

    def __add__(self, other):
        if isinstance(other, CFrac):
            return BiHomography([0, mpz(1), mpz(1), 0], [0, 0, 0, mpz(1)])(self, other)

        elif isinstance(other, numbers.Real):
            return Homography(mpz(1), other, mpz(0), mpz(1))(self)

        else:
            return NotImplemented

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return Homography(mpz(-1), 0, 0, mpz(1))(self)

    def __pos__(self):
        return self

    def __sub__(self, other):
        if isinstance(other, CFrac):
            return BiHomography([0, mpz(1), mpz(-1), 0], [0, 0, 0, mpz(1)])(self, other)

        elif isinstance(other, numbers.Real):
            return Homography(mpz(1), -other, mpz(0), mpz(1))(self)

        else:
            return NotImplemented

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if isinstance(other, CFrac):
            return BiHomography([mpz(1), 0, 0, 0], [0, 0, 0, mpz(1)])(self, other)

        elif isinstance(other, numbers.Real):
            return Homography(other, mpz(0), mpz(0), mpz(1))(self)

        else:
            return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, CFrac):
            return BiHomography([0, mpz(1), 0, 0], [0, 0, mpz(1), 0])(self, other)

        elif isinstance(other, numbers.Real):
            return Homography(mpz(1), mpz(0), mpz(0), other)(self)

        else:
            return NotImplemented

    def __rtruediv__(self, other):
        if self[0] != 0:
            return CFrac(it.chain((0,), self)) * other
        elif self.longer_than_eq(2):
            return self[1:] * other
        else:
            raise ZeroDivisionError("division or modulo by Zero")

    def __floordiv__(self, other):
        return (self / other)[0]

    def __rfloordiv__(self, other):
        return (other / self)[0]

    def divmod(self, other):
        q = self // other
        r = self - q*other
        return q, r

    def __mod__(self, other):
        q = self // other
        return self - q*other

    def __rmod__(self, other):
        q = other // self
        return other - q*self

    def __pos__(self):
        return self

    def to_frac(self, depth=DEPTH) -> Fraction:
        A, A_ = self[0], 1
        B, B_ = 1, 0

        if depth > 1:
            for coeff in self[1 : depth]:
                A, A_ = coeff * A + A_, A
                B, B_ = coeff * B + B_, B

        return Fraction(A, B)

    def gen_convergents(self) -> Iterable[Fraction]:
        A, A_ = self[0], 1
        B, B_ = 1, 0
        yield Fraction(A, B)

        if self.longer_than_eq(2):
            for coeff in self[1:]:
                A, A_ = coeff * A + A_, A
                B, B_ = coeff * B + B_, B
                yield Fraction(A, B)

    def __pow__(self, power): # exponentiation by squaring
        if isinstance(power, numbers.Integral):
            if power == 1:
                return self
            factor = self ** (power // 2)

            return factor * factor if power % 2 else factor * factor * self
        else:
            return NotImplemented

    def __rpow__(self, base):
        return NotImplemented

    def __divmod__(self, other):
        q = self // other
        return q, self - q*other

    def __rdivmod__(self, other):
        q = other // self
        return q, other - q*self


class Homography:
    """
    Homography(1, 2, 3, 4)

    --->

    1x + 2
    ------
    3x + 4
    """
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def __repr__(self):
        return "{0.a} x + {0.b}\n---------\n{0.c} x + {0.d}".format(self)

    def digit_factory(self, cf: CFrac) -> int:
        """
        The homography

        a x + b
        -------
        c x + d

        if apllied to the continued fraction cf. The result is returned as a CF
        Algorithm from (Gosper, 1972)
        """
        a, b, c, d = self.a, self.b, self.c, self.d
        cfiter = iter(cf)
        while True:
            if (a, b, c, d) == (1, 0, 0, 1): # 1 0 0 1 is the identity homography
                yield from cfiter

            try:
                x = next(cfiter)
                a, b = a * x + b, a
                c, d = c * x + d, c
            except StopIteration:
                yield from euclid_factors(a, c)
                return

            while (c, d) != (0, 0):
                q1 = a // c if c else math.inf

                q2 = (a + b) // (c + d) if c + d else math.inf

                if q1 == q2:
                    yield q2
                    (a, b), (c, d) = (c, d), (a - q2*c, b - q2*d)
                else:
                    break
            else: # nobreak, d h. bb = 0
                raise StopIteration()

    def __call__(self, cf):
        return CFrac(self.digit_factory(cf))


def BiHom_empty(num, denom, yiter):
    a1, a2, a3, a4 = num
    b1, b2, b3, b4 = denom

    yield from Homography(a1, a2, b1, b2)(yiter)


class BiHomography:
    """
    BiHomography([1, 2, 3, 4], [5, 6, 7, 8])

    --->

    1*xy + 2x + 3y + 4
    ------------------
     5xy + 6x + 7y + 8
    """

    def __init__(self, num: List[int], denom: List[int]):
        self.num = num
        self.denom = denom

    def __repr__(self):
        return "{0.num[0]} xy + {0.num[1]} x + {0.num[2]} y + {0.num[3]}\n" \
               "----------------------\n" \
               "{0.denom[0]} xy + {0.denom[1]} x + {0.denom[2]} y + {0.denom[3]}".format(self)

    def digit_factory(self, cfx: CFrac, cfy: CFrac) -> Iterable[int]:
        """
        Die "Bihomographie" wird auf die beiden Kettenbrüche cfx und cfy angewandt

        Algorithmus aus (Gosper, 1972)
        """
        a1, a2, a3, a4 = self.num
        b1, b2, b3, b4 = self.denom

        xiter = iter(cfx)
        yiter = iter(cfy)
        while True:
            try:
                x = next(xiter)
            except StopIteration:
                yield from BiHom_empty([a1, a2, a3, a4], [b1, b2, b3, b4], yiter)
                return

            try:
                y = next(yiter)
            except StopIteration:
                yield from BiHom_empty([a1, a3, a2, a4], [b1, b3, b2, b4],
                                       it.chain(iter([x]), xiter)) # letztes x "zurückstecken"
                return

            a1, a2, a3, a4 = (a1*x*y + a2*x + a3*y + a4, a1*x + a3, a1*y + a2, a1)
            b1, b2, b3, b4 = (b1*x*y + b2*x + b3*y + b4, b1*x + b3, b1*y + b2, b1)

            while (b1, b2, b3, b4) != (0, 0, 0, 0):

                q1 = a1 // b1 if b1 else math.inf

                q2 = (a1 + a2) // (b1 + b2) if (b1 + b2) else math.inf
                if q1 != q2: break

                q3 = (a1 + a3) // (b1 + b3) if (b1 + b3) else math.inf
                if q2 != q3: break

                q4 = (a1 + a2 + a3 + a4) // (b1 + b2 + b3 + b4) if (b1 + b2 + b3 + b4)\
                                                                else math.inf

                if q1 == q4:
                    yield q1
                    a1_, a2_, a3_, a4_ = a1, a2, a3, a4
                    a1, a2, a3, a4 = b1, b2, b3, b4
                    b1, b2, b3, b4 = a1_ - q1*b1, a2_ - q1*b2, a3_ - q1*b3, a4_ - q1*b4

                else:
                    break
            else: # nobreak, d h.  b1 = b2 = b3 = b4 = 0
                raise StopIteration()

    def __call__(self, cfx, cfy):
        return CFrac(self.digit_factory(cfx, cfy))

def e_gen():
    """
    An iterator which yields the coefficients of the continued fraction for e
    """
    yield 2
    i = 2
    while True:
        if i%3:
            yield 1
        else:
            yield i//3 * 2
        i += 1



if __name__ == '__main__':
    # tests / examples

    pi = CFrac(math.pi)
    e = CFrac(e_gen())
    print('pi     = {}'.format(pi))
    assert float(CFrac(math.pi)) - math.pi < 1e-15
    print('e      = {}'.format(e))
    assert float(CFrac(math.e)) - math.e < 1e-15

    print('2pi    = {}'.format(2*pi))
    assert float(2*pi) - 2*math.pi < 1e-15
    print('pi/3   = {}'.format(pi/3))
    assert float(pi/3) - math.pi/3 < 1e-15
    print('1/pi   = {}'.format(1/pi))
    assert float(1/pi) - 1/math.pi < 1e-15


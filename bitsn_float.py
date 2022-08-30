from pymtl3 import *


# from copy import deepcopy
# import traceback
# import sys
def real_to_bits(real, frac_bits):
    # print("---------")
    if real == 0:
        return 0, []
    frac_bits += 1  # for float normal value hidden 1
    real = abs(real)
    bits = []
    exp = 0
    (int_p, frac_p) = divmod(real, 1)
    # print(int_p, frac_p)
    exp = 0
    while int_p > 0:
        (int_p, bit) = divmod(int_p, 2)
        # print(int_p, bit)
        bits.append(int(bit))
        exp += 1
        # print(exp, bits)
    bits.reverse()
    # print("AF int_P", exp, bits)
    if exp < frac_bits:
        if len(bits) == 0:
            while True:
                (bit, frac_p) = divmod(frac_p * 2, 1)
                if (bit == 0):
                    exp -= 1
                else:
                    bits.append(1)
                    break
                # print("zero padding", exp, bits, frac_p)

        # print("AF zero pad", exp, bits, frac_p)
        for i in range(frac_bits - len(bits)):
            (bit, frac_p) = divmod(frac_p * 2, 1)
            # print(frac_p,bit)
            bits.append(int(bit))
            # print(exp, bits, frac_p)
    # shift out hidden 1
    exp -= 1
    bits.pop(0)
    # print(exp, bits)
    return exp, bits


def rte(bits, exp, frac_bitsn):
    rte_deci = ''.join(map(str, bits[:frac_bitsn]))
    rte_carry = 0
    if bits[frac_bitsn]:
        # print("tail_bits",bits[frac_bitsn+1:])
        if 1 in bits[frac_bitsn + 1:]:
            rte_carry = 1
        elif bits[frac_bitsn - 1] == 1:
            rte_carry = 1
    # print("rte_carry:",rte_carry)
    frac_class = mk_bits(frac_bitsn)
    frac = frac_class(int(rte_deci, 2))
    if rte_carry:
        if (frac + 1).uint() == 0:
            exp += 1  # rte carry overflow
            frac = frac_class(0)
        else:
            frac = frac + 1

    return exp, frac


def real_to_bitsn_float(real, exp_bitsn, frac_bitsn, prec_bitsn):
    exp_class = mk_bits(exp_bitsn)
    frac_class = mk_bits(frac_bitsn)
    temp_class = mk_bits(1 + exp_bitsn + frac_bitsn)

    exp, bits = real_to_bits(real, frac_bitsn + prec_bitsn)
    # round to even
    exp, frac = rte(bits, exp, frac_bitsn)

    exp = exp_class(exp + (2 ** (exp_bitsn - 1) - 1))
    temp_ins = temp_class(concat(Bits1(0), exp, frac) if real >= 0 else concat(Bits1(1), exp, frac))
    if (real == 0):
        return bitsn_float(exp_class(0), frac_class(0), temp_class(0), prec_bitsn)
    else:
        return bitsn_float(exp_bitsn, frac_bitsn, temp_ins, prec_bitsn)


class bitsn_float(object):
    def __init__(self, exp_bits_num, frac_bits_num, ini_val=None, cal_prec_ext_bits_num=8):
        self.exp_bitsn = exp_bits_num
        self.frac_bitsn = frac_bits_num
        self.cal_prec_ext_bits_num = cal_prec_ext_bits_num
        self.bits = 1 + exp_bits_num + frac_bits_num
        self.float_class = mk_bitstruct(f"Float{self.bits}", {
            'sign': mk_bits(1),
            'exp' : mk_bits(exp_bits_num),
            'frac': mk_bits(frac_bits_num),
        })
        if (ini_val):
            try:
                self.float_bits = self.float_class.from_bits(ini_val)
            except Exception as e:
                print(e)
        else:
            temp_cl = mk_bits(self.bits)
            self.float_bits = self.float_class.from_bits(temp_cl(0))

    def from_bits(self, bitsval):
        try:
            self.float_bits = self.float_class.from_bits(bitsval)
            # print("DEBUG:", self.float_bits.to_bits())
        except Exception as e:
            print(e)

    def to_bits(self):
        return self.float_bits.to_bits()

    def exp(self):
        return self.float_bits.exp

    def frac(self):
        return self.float_bits.frac

    def sign(self):
        return self.float_bits.sign

    def real(self):
        if self.float_bits.exp > 0:
            exp_real = self.float_bits.exp.uint() - (2 ** (self.exp_bitsn - 1) - 1)
            frac_real = self.float_bits.frac.uint() * (2 ** (-1 * self.frac_bitsn)) + 1
        else:
            exp_real = -1 * (2 ** (self.exp_bitsn - 1) - 2)
            frac_real = self.float_bits.frac.uint() * (2 ** (-1 * self.frac_bitsn))
        return pow(-1, self.float_bits.sign.uint()) * pow(2, exp_real) * frac_real

    def prec_get(self, other):
        exp_bitsn = self.exp_bitsn if self.exp_bitsn > other.exp_bitsn else other.exp_bitsn
        frac_bitsn = self.frac_bitsn if self.frac_bitsn > other.frac_bitsn else other.frac_bitsn
        cal_prec_ext_bits_num = self.cal_prec_ext_bits_num if self.cal_prec_ext_bits_num > other.cal_prec_ext_bits_num else other.cal_prec_ext_bits_num
        return exp_bitsn, frac_bitsn, cal_prec_ext_bits_num

    def __str__(self):
        return "0x{}({})".format(str(self.float_bits.to_bits()), str(self.real()))

    def __add__(self, other):
        temp_add = self.real() + other.real()
        exp_bitsn, frac_bitsn, cal_prec_ext_bits_num = self.prec_get(other)
        return real_to_bitsn_float(temp_add, exp_bitsn, frac_bitsn, cal_prec_ext_bits_num)

    def __mul__(self, other):
        temp_add = self.real() * other.real()
        exp_bitsn, frac_bitsn, cal_prec_ext_bits_num = self.prec_get(other)
        return real_to_bitsn_float(temp_add, exp_bitsn, frac_bitsn, cal_prec_ext_bits_num)

    def __truediv__(self, other):
        temp_add = other.real() * (1.0 / self.real()) if self.frac_bitsn > other.frac_bitsn else self.real() * (
                1.0 / other.real())
        # print(self.real(),other.real(),1.0 / other.real(),temp_add)
        exp_bitsn, frac_bitsn, cal_prec_ext_bits_num = self.prec_get(other)
        return real_to_bitsn_float(temp_add, exp_bitsn, frac_bitsn, cal_prec_ext_bits_num)

    def __neg__(self):
        # print("neg called")
        return self.__mul__(real_to_bitsn_float(-1, self.exp_bitsn, self.frac_bitsn, self.cal_prec_ext_bits_num))


if __name__ == '__main__':
    # print(real_to_bits(0.125 + 0.3125, 10))
    # print(real_to_bits(13.3, 10))
    # print(real_to_bits(0, 10))
    # print(real_to_bits(1, 10))
    # print(real_to_bits(-1, 10))
    print(real_to_bits(0.000000003, 5))

    print("main--begin")
    a = bitsn_float(8, 15)
    bc = mk_bits(24)
    b = bc(0x417000)
    print(type(b))
    a.from_bits(b)
    print(a.to_bits())
    print(type(a.to_bits()))
    print(a.real())
    print(a)
    c = bitsn_float(8, 23, Bits32(0x3e810000))
    print(c.to_bits())
    print(type(c.to_bits()))
    print(c.real())
    print(c)
    print("exp", a.exp())
    print("exp", c.exp())
    print(type(c.exp()))
    print(type(c.frac()))
    d = a + c
    print(type(d))
    print(d)
    print(d.real())
    print(d.to_bits())
    d = a * c
    print(type(d))
    print(d)
    print(d.real())
    print(d.to_bits())
    e = -d
    print(e)
    print(e.real())
    print(e.to_bits())
    f = a / c
    print(f)
    print(f.real())
    print(f.to_bits())

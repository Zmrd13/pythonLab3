from hashlib import md5
import math
import random
import sys
import typing

def fast_pow(base: int, power: int, mod: int = None) -> int:

    if mod == 0:
        raise ValueError("Модуль не может быть равен нулю")
    if power < 0:
        raise ValueError("Показатель не может быть отрицательным")
    result = 1
    if mod:
        base %= mod
        result %= mod  # Для правильного вычисления, когда показатель равен 0, а модуль отрицателен
    while power:
        power, remainder = divmod(power, 2)
        if remainder:
            result *= base
            if mod:
                result %= mod
        base *= base
        if mod:
            base %= mod
    return result


def extgcd(a: int, b: int) -> typing.Tuple[int, int, int]:
    if a <= 0 or b <= 0:
        raise ValueError("Числа могут быть только натуральными")
    # if a < b:
    #     a, b = b, a  # Ломает алгоритм
    u1, u2, u3 = a, 1, 0
    v1, v2, v3 = b, 0, 1
    while v1:
        q = u1 // v1
        t1, t2, t3 = u1 % v1, u2 - q * v2, u3 - q * v3
        u1, u2, u3 = v1, v2, v3
        v1, v2, v3 = t1, t2, t3
    return u1, u2, u3


def is_prime(n, trials=8):

    if n != int(n):
        return False
    n = int(n)
    # Miller-Rabin test for prime
    if n == 0 or n == 1 or n == 4 or n == 6 or n == 8 or n == 9:
        return False

    if n == 2 or n == 3 or n == 5 or n == 7:
        return True
    s = 0
    d = n - 1
    while d % 2 == 0:
        d >>= 1
        s += 1
    assert (2 ** s * d == n - 1)

    def trial_composite(a):
        if pow(a, d, n) == 1:
            return False
        for i in range(s):
            if pow(a, 2 ** i * d, n) == n - 1:
                return False
        return True

    for i in range(trials):  # number of trials
        a = random.randrange(2, n)
        if trial_composite(a):
            return False

    return True


def gen_p(a: int, b: int) -> int:

    while True:
        p = random.randint(a, b)
        if is_prime(p):
            return p


def gen_safe_p(a: int, b: int) -> int:

    while True:
        q = gen_p(a // 2, (b - 1) // 2)
        if is_prime(q * 2 + 1):
            return q * 2 + 1


def gen_g(mod: int) -> int:
    while True:
        g = random.randrange(2, mod)
        if pow(g, (mod - 1) // 2, mod) != 1:
            return g


def gen_public(private_key: int, mod: int):

    return pow(gen_g(mod), private_key, mod)


def gen_common(secret_key: int, public_key: int, mod: int) -> int:

    return pow(public_key, secret_key, mod)


def shanks(y: int, a: int, mod: int) -> typing.Union[int, None]:

    if y >= mod:
        raise ValueError("y не может быть больше или равным mod")
    m = k = math.ceil(math.sqrt(mod))
    seq1 = {pow(a, j, mod) * y % mod: j for j in range(m)}
    seq2 = (pow(a, i * m, mod) for i in range(1, k + 1))
    for i, vel in enumerate(seq2, 1):
        if (j := seq1.get(vel)) is not None:
            return i * m - j
    return None


def gen_mutually_prime(a):
    while True:
        b = random.randrange(2, a)
        if math.gcd(a, b) == 1:
            return b


def write_bytes_to_file(bytes, file):
    for byte in bytes:
        file.write(byte.to_bytes(1, sys.byteorder))


def gen_c_d(p):
    c = gen_mutually_prime(p)
    gcd, d, _ = extgcd(c, p)
    assert gcd == 1
    while d < 0:
        d += p
    return c, d


def inverse(n, p):
    gcd, inv, _ = extgcd(n, p)
    assert gcd == 1
    return inv








print("                             GOST                                   ")

q = gen_p(1 << 255, (1 << 256) - 1)

while True:  # Генерируем p
    b = random.randint(math.ceil((1 << 1023) / q), ((1 << 1024) - 1) // q)
    if is_prime(p := b * q + 1):
        break

while True:  # Находим a
    g = random.randrange(2, p - 1)
    if (a := pow(g, b, p)) > 1:
        break

x = random.randrange(1, q)  # Закрытый ключ
y = pow(a, x, p)  # Открытый ключ

h = int.from_bytes(md5(open('test.txt', 'rb').read()).digest(), byteorder=sys.byteorder)
assert 0 < h < q
while True:
    k = random.randrange(1, q)
    if (r := pow(a, k, p) % q) == 0:
        continue
    if (s := (k * h % q + x * r % q) % q) != 0:
        break
open('gen/gost.txt', 'w').write(f'{r = }\n{s = }')

# Проверка подписи
assert 0 < r < q
assert 0 < s < q
gcd, hh, _ = extgcd(h, q)
u1 = s * hh % q
u2 = -r * hh % q
v = pow(a, u1, p) * pow(y, u2, p) % p % q
assert v == r
print(f'{q = }\n{p = }\n{a = }\n{x = }\n{y = }\n{h = }\n{k = }\n{s = }\n{u1 = }\n{u2 = }\n{r = }\n{v = }')


print("                             GAMAL                                   ")

h = md5(open('test.txt', 'rb').read())
h_int = int.from_bytes(h.digest(), byteorder=sys.byteorder)  # Хеш файла

p = gen_safe_p(1 << h.digest_size * 8 + 1, 1 << h.digest_size * 8 + 2)

g = gen_g(p)
x = random.randrange(2, p - 1)  # Закрытый ключ
y = pow(g, x, p)  # Открытый ключ
k = gen_mutually_prime(p - 1)
r = pow(g, k, p)
u = (h_int - x * r) % (p - 1)
gcd, kk, _ = extgcd(k, p - 1)
assert gcd == 1
s = kk * u % (p - 1)
open('gen/elgamal.txt', 'w').write(str(s))

# Проверка подписи
yr = pow(y, r, p) * pow(r, s, p) % p
gh = pow(g, h_int, p)
assert yr == gh
print(f'{h_int = }\n{p = }\n{g = }\n{x = }\n{y = }\n{k = }\n{r = }\n{u = }\n{s = }\n{yr = }\n{gh = }')
print("                             RSA                                   ")

h = md5(open('test.txt', 'rb').read())
h_int = int.from_bytes(h.digest(), byteorder=sys.byteorder)  # Хеш файла

p = gen_p(1 << h.digest_size * 8 // 2 + 1, 1 << h.digest_size * 8 // 2 + 2)
q = gen_p(1 << h.digest_size * 8 // 2 + 1, 1 << h.digest_size * 8 // 2 + 2)
assert p != q
n = p * q
phi = (p - 1) * (q - 1)
d, c = gen_c_d(phi)  # Открытый и закрытый ключ

s = pow(h_int, c, n)  # Подпись
open('gen/rsa.txt', 'w').write(str(s))
# Проверка подписи
e = pow(s, d, n)
assert e == h_int
print(f'{p = }\n{q = }\n{n = }\n{phi = }\n{d = }\n{c = }\n{h_int = }\n{s = }\n{e = }')


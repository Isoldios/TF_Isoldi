# %%
#-- Definiciones de parámetros --#
#-- y^2 = x^3 + a*x + b mod(p) --#
p = 433
a = 0
b = 7

# %%
#-- Operaciones modulares --#
def inv_mod(x, p):
    return pow(x % p, p - 2, p)

def legendre_symbol(a, p):
    return pow(a % p, (p - 1) // 2, p)

def modular_sqrt(n, p):
    n %= p
    if n == 0:
        return 0
    if p == 2:
        return n
    ls = legendre_symbol(n, p)
    if ls == p - 1:
        return None
    if p % 4 == 3:
        r = pow(n, (p + 1) // 4, p)
        return r
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre_symbol(z, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    while True:
        if t == 0:
            return 0
        if t == 1:
            return r
        t2i = t
        i = 0
        for i_check in range(1, m):
            t2i = pow(t2i, 2, p)
            if t2i == 1:
                i = i_check
                break
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p

# %%
#-- Operaciones de curvas elípticas --#
def on_curve(P):
    if P is None:
        return True
    x, y = P
    return (y*y - (x*x*x + a*x + b)) % p == 0

def point_neg(P):
    if P is None:
        return None
    x, y = P
    return (x, (-y) % p)

def point_add(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P != Q:
        den = (x2 - x1) % p
        m = ((y2 - y1) * inv_mod(den, p)) % p
    else:
        if y1 % p == 0:
            return None
        den = (2 * y1) % p
        m = ((3 * x1 * x1 + a) * inv_mod(den, p)) % p
    x3 = (m*m - x1 - x2) % p
    y3 = (m*(x1 - x3) - y1) % p
    return (x3, y3)

def scalar_mult(k, P):
    R = None
    N = P
    kk = k
    while kk > 0:
        if kk & 1:
            R = point_add(R, N)
        N = point_add(N, N)
        kk >>= 1
    return R

# %%
#-- Funciones principales --#
def enumerate_points():
    pts = []
    for x in range(p):
        rhs = (x*x*x + a*x + b) % p
        y = modular_sqrt(rhs, p)
        if y is None:
            continue
        y2 = (-y) % p
        pts.append((x, y))
        if y2 != y:
            pts.append((x, y2))
    return pts

def factor_trial(n):
    # devuelve dict {primo: exponente}
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors

def order_from_group_order(P, group_order, group_factors):
    if P is None:
        return 1
    ord_candidate = group_order
    for q, exp in group_factors.items():
        for _ in range(exp):
            if ord_candidate % q != 0:
                break
            test = ord_candidate // q
            if scalar_mult(test, P) is None:
                ord_candidate = test
            else:
                break
    return ord_candidate

def find_max_order_generator(stop_if_full_generator=True):
    points = enumerate_points()
    group_order = len(points) + 1
    group_factors = factor_trial(group_order)
    max_ord = 0
    best_pts = []
    for idx, P in enumerate(points):
        if P[1] == 0:
            ordP = order_from_group_order(P, group_order, group_factors)
        else:
            ordP = order_from_group_order(P, group_order, group_factors)
        if ordP is None:
            continue
        if ordP > max_ord:
            max_ord = ordP
            best_pts = [P]
        elif ordP == max_ord:
            best_pts.append(P)
        if stop_if_full_generator and max_ord == group_order:
            break
    return group_order, max_ord, best_pts

# %%
#-- Ejecutable --#
if __name__ == "__main__":
    group_order, max_order, gens = find_max_order_generator()
    print(f"#E(F_{p}) = {group_order}")
    print(f"Mayor orden de punto encontrado: {max_order}")
    print(f"Cantidad de puntos con ese orden: {len(gens)}")
    for i, P in enumerate(gens[:10], 1):
        print(f"{i:2d}. G = {P}")
    if len(gens) > 10:
        print(f"... y {len(gens)-10} más")

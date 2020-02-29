grav = 981 *m/s**2 / 100

mass, a, mu = var("mass, a, mu")

sub_list = [
    (mass, 800 *kg),
    (a, S(5)/10 *m/s**2),
    (mu, S(1)/10),
    ]

G = mass*grav
c = sqrt(2)/2

Nv, S = var("Nv, S")
eq1 = Eq(mass*a, S*c - mu*Nv)
eq2 = Eq(G, Nv + S*c)
sol = solve([eq1, eq2],[Nv, S])

pprint("\nN / Newton:")
Nv = sol[Nv]
tmp = Nv
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= Newton
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\nS / Newton:")
S = sol[S]
tmp = S
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp /= Newton
tmp = iso_round(tmp,0.1)
pprint(tmp)

pprint("\nphi / deg:")
phi = asin(mu * Nv / S)
tmp = phi
tmp = tmp.subs(sub_list)
tmp = tmp.simplify()
tmp = tmp * 180 / pi
tmp = iso_round(tmp,0.001)
pprint(tmp)

pprint("\ntheta / deg:")
tmp = 45 - tmp
tmp = iso_round(tmp,0.001)
pprint(tmp)

# N / Newton:
# 6770.9
#
# S / Newton:
# 1523.2
#
# phi / deg:
# 26.392
#
# theta / deg:
# 18.608

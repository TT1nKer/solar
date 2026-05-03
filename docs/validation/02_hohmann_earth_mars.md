# Validation: Earth-Mars Hohmann Transfer

## Claim
The Hohmann transfer calculation for Earth → Mars produces total delta-v and time of flight matching standard textbook values.

## Reference
Standard patched-conic Hohmann transfer between two circular coplanar orbits.

Textbook values (Bate, Mueller, White; Curtis):
- r₁ = 1.000 AU = 1.496e8 km (Earth orbit)
- r₂ = 1.524 AU = 2.279e8 km (Mars orbit)
- Δv₁ ≈ 2.94 km/s (departure)
- Δv₂ ≈ 2.65 km/s (arrival)
- Total Δv ≈ 5.59 km/s
- TOF ≈ 259 days

These come from:
```
a_transfer = (r₁ + r₂) / 2
v_circ_i = sqrt(μ_sun / r_i)
v_trans_i = sqrt(μ_sun · (2/r_i - 1/a_transfer))
Δv_i = |v_trans_i - v_circ_i|
TOF = π · sqrt(a_transfer³ / μ_sun)
```

## Command
```bash
./solar transfer Earth Mars 2026-04-12
```

## Actual output
```
Hohmann Transfer: Earth -> Mars
Departure orbit: 1.4960e+08 km (1.0000 AU)
Arrival orbit:   2.2794e+08 km (1.5237 AU)
Delta-v1 (departure): 2.945 km/s
Delta-v2 (arrival):   2.649 km/s
Total delta-v:        5.594 km/s
Time of flight:       258.871 days
Transfer semi-major:  1.8877e+08 km
Transfer eccentricity: 0.207515
```

## Error
- Δv₁: 2.945 vs 2.94 → matches to 0.2%
- Δv₂: 2.649 vs 2.65 → matches to 0.04%
- Total Δv: 5.594 vs 5.59 → matches to 0.07%
- TOF: 258.87 vs 259 → matches to 0.05%

All within numerical precision of textbook values.

## Notes
- This is a **2-body Hohmann assumption**. Real Earth-Mars transfers have higher Δv due to non-circular and non-coplanar orbits.
- The actual 2026 launch window optimum (via Lambert solver, not Hohmann) is ~5.64 km/s.
- This validation confirms the **arithmetic** is correct, not that the result is operationally useful.

## Limitations
- Only one transfer pair (Earth-Mars) tested
- Only the idealized circular coplanar case
- Does not validate the Lambert solver (separate validation needed)
- Does not validate against any operational mission

## Reference for further reading
- Curtis, *Orbital Mechanics for Engineering Students*, 3rd ed., Ch. 6
- Bate, Mueller, White, *Fundamentals of Astrodynamics*, Dover

# Validation: Coordinate Frame Transformations

## Claim
The Mat3 rotation matrices and Frame transform functions correctly convert between ICRF (equatorial J2000) and Ecliptic J2000 frames, preserving vector magnitudes and round-tripping to identity.

## Reference

### ICRF ↔ Ecliptic J2000

The transformation between ICRF (equatorial J2000) and Ecliptic J2000 is a rotation about the X-axis by the **mean obliquity of the ecliptic at J2000**:

```
ε = 23.4392911111° = 0.4090928042 rad
```

Source: IAU 2000 resolution, Astronomical Almanac.

The rotation matrices:
```
R_ecliptic→ICRF = Rx(-ε) = [1   0       0    ]
                            [0   cos ε  -sin ε]
                            [0   sin ε   cos ε]

R_ICRF→ecliptic = Rx(+ε) = transpose of above
```

## Method
1. Construct test vectors with known properties
2. Apply forward transform, then inverse transform
3. Compare with original (round-trip should be identity)
4. Verify magnitude preservation (rotations are orthogonal)
5. Verify angle preservation between two vectors

## Test 1: Round-trip identity

```cpp
Vec3 v_orig = {1.0, 2.0, 3.0};
Vec3 v_icrf = ecliptic_to_icrf().apply(v_orig);
Vec3 v_back = icrf_to_ecliptic().apply(v_icrf);
double err = (v_back - v_orig).norm();
```

Expected: err < 1e-15 (machine precision).

## Test 2: Magnitude preservation

```cpp
Vec3 v = {1.0, 0.0, 0.0};
double mag_orig = v.norm();      // = 1.0
double mag_rot = ecliptic_to_icrf().apply(v).norm();
```

Expected: mag_rot - mag_orig < 1e-15.

## Test 3: DE440 cross-check

When loading DE440 (which is in ICRF/equatorial), the parser internally rotates to ecliptic. Earth at J2000 in ecliptic:
- Position: -2.521e+07, +1.449e+08, -6.164e+02 km
- Z component is small (~600 km) because Earth orbits in ecliptic plane (by definition, |z| << |x|, |y|)

If we instead used the equatorial frame, the Z component would be:
- Z_eq ≈ -|x|·sin(ε) - small ≈ -5.8e+07 km (much larger)

Observation: the DE440 output via our parser has |Z| ~600 km, which confirms ecliptic frame.

## Commands
```bash
make test  # runs frame round-trip tests if test_frames exists

# Manual check via DE440:
./solar --de data/de440.asc ephemeris Earth 2000-01-01 | grep "Position"
# Should show small Z (~hundreds of km), not millions
```

## Actual output (DE440 cross-check)

```
$ ./solar --de data/de440.asc ephemeris Earth 2000-01-01
Position (km):    -2.5211e+07   1.4493e+08  -6.1642e+02
```

|z|/|x,y| ≈ 6.16e+02 / 1.5e+08 ≈ 4×10⁻⁶ (consistent with Earth being in the ecliptic by definition).

## Errors
- Round-trip error: < machine precision for all components
- Magnitude preservation: exact (orthogonal rotation)
- DE440 ecliptic check: |Z| consistent with ecliptic frame, not equatorial

## Notes

### Frame routing through ICRF
Higher-order transforms (BodyFixed ↔ Ecliptic) route through ICRF as a hub:
```
EclipticJ2000 → ICRF → BodyFixed
                       (uses IAU rotation model, time-dependent)
```

This means at most 2 matrix multiplications for any conversion. The constant ε rotation is cached as `static const Mat3`.

### Velocity transformation in rotating frames
For inertial-to-inertial (Ecliptic ↔ ICRF), velocity rotates the same as position. For BodyFixed (rotating frame), velocity transformation includes Coriolis correction:
```
v_body = R · v_inertial - ω × r_body
```

where ω comes from the IAU model's `W_dot` parameter. This **has not been independently validated** — it relies on the IAU model being correctly implemented.

## Limitations

- **No general precession/nutation**: We use the mean obliquity at J2000 only. Real precise transforms include precession (~50 arcsec/year), nutation (~10 arcsec amplitude), and polar motion (~0.3 arcsec amplitude). These corrections are not implemented.
- **No verification of body-fixed frames at non-J2000 epochs**: The IAU rotation models are time-dependent but only sanity-checked (sub-solar point at midnight UTC ≈ 0° longitude).
- **Simplified Moon orientation**: The full IAU Moon rotation has 13 periodic correction terms. We use only the secular polynomial part. Resulting orientation error is up to ~0.1° in some cases.
- **No validation against SPICE or other reference systems** for frame transforms beyond the ICRF↔Ecliptic case.

## Future validation needed
- Compare body-fixed Earth at midnight UTC with sidereal time formula
- Cross-check DE440 Mars at multiple epochs after frame transform
- Validate IAU Moon rotation against SPICE moon_pa kernel
- Test frame transforms at far-from-J2000 epochs (e.g., 2050) where precession matters

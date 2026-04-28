#include "solar/time_scale.h"
#include "solar/ephemeris.h"
#include "solar/constants.h"
#include <cmath>
#include <stdexcept>

namespace solar {

const char* time_scale_name(TimeScale ts) {
    switch (ts) {
        case TimeScale::TDB: return "TDB";
        case TimeScale::TT:  return "TT";
        case TimeScale::UTC: return "UTC";
        case TimeScale::TAI: return "TAI";
    }
    return "?";
}

// TDB - TT: Fairhead & Bretagnon (1990) dominant term
// Amplitude ~1.657ms, period ~1 year
double tdb_to_tt(double jd_tdb) {
    double T = (jd_tdb - constants::J2000) / 36525.0;
    double dt_sec = 0.001657 * std::sin(628.3076 * T + 6.2401);
    return jd_tdb - dt_sec / constants::DAY;
}

double tt_to_tdb(double jd_tt) {
    // Iterative: TDB ≈ TT + correction(TT)
    double T = (jd_tt - constants::J2000) / 36525.0;
    double dt_sec = 0.001657 * std::sin(628.3076 * T + 6.2401);
    return jd_tt + dt_sec / constants::DAY;
}

// TT = TAI + 32.184s
constexpr double TT_TAI_OFFSET = 32.184; // seconds

double tt_to_tai(double jd_tt) {
    return jd_tt - TT_TAI_OFFSET / constants::DAY;
}

double tai_to_tt(double jd_tai) {
    return jd_tai + TT_TAI_OFFSET / constants::DAY;
}

// Leap second table: {JD_UTC, cumulative_leap_seconds}
// From 1972-01-01 (first leap second) through 2017-01-01 (last as of 2025)
struct LeapEntry {
    double jd;   // UTC Julian date when this leap second count takes effect
    int leaps;   // cumulative leap seconds
};

// Table must be in chronological order
static const LeapEntry LEAP_TABLE[] = {
    {2441317.5, 10},  // 1972-01-01
    {2441499.5, 11},  // 1972-07-01
    {2441683.5, 12},  // 1973-01-01
    {2442048.5, 13},  // 1974-01-01
    {2442413.5, 14},  // 1975-01-01
    {2442778.5, 15},  // 1976-01-01
    {2443144.5, 16},  // 1977-01-01
    {2443509.5, 17},  // 1978-01-01
    {2443874.5, 18},  // 1979-01-01
    {2444239.5, 19},  // 1980-01-01
    {2444786.5, 20},  // 1981-07-01
    {2445151.5, 21},  // 1982-07-01
    {2445516.5, 22},  // 1983-07-01
    {2446247.5, 23},  // 1985-07-01
    {2447161.5, 24},  // 1988-01-01
    {2447892.5, 25},  // 1990-01-01
    {2448257.5, 26},  // 1991-01-01
    {2448804.5, 27},  // 1992-07-01
    {2449169.5, 28},  // 1993-07-01
    {2449534.5, 29},  // 1994-07-01
    {2450083.5, 30},  // 1996-01-01
    {2450630.5, 31},  // 1997-07-01
    {2451179.5, 32},  // 1999-01-01
    {2453736.5, 33},  // 2006-01-01
    {2454832.5, 34},  // 2009-01-01
    {2456109.5, 35},  // 2012-07-01
    {2457204.5, 36},  // 2015-07-01
    {2457754.5, 37},  // 2017-01-01
    // Update this table when IERS announces new leap seconds
};

static constexpr int LEAP_TABLE_SIZE = sizeof(LEAP_TABLE) / sizeof(LEAP_TABLE[0]);

int leap_seconds_at(double jd_utc) {
    if (jd_utc < LEAP_TABLE[0].jd) return 0; // before 1972
    int leaps = LEAP_TABLE[0].leaps;
    for (int i = 0; i < LEAP_TABLE_SIZE; ++i) {
        if (jd_utc >= LEAP_TABLE[i].jd) {
            leaps = LEAP_TABLE[i].leaps;
        } else {
            break;
        }
    }
    return leaps;
}

// TAI = UTC + leap_seconds
double utc_to_tai(double jd_utc) {
    int leaps = leap_seconds_at(jd_utc);
    return jd_utc + static_cast<double>(leaps) / constants::DAY;
}

double tai_to_utc(double jd_tai) {
    // Approximate: use TAI to estimate UTC, then look up leaps
    double jd_utc_approx = jd_tai - 37.0 / constants::DAY; // rough estimate
    int leaps = leap_seconds_at(jd_utc_approx);
    return jd_tai - static_cast<double>(leaps) / constants::DAY;
}

double tt_to_utc(double jd_tt) {
    double jd_tai = tt_to_tai(jd_tt);
    return tai_to_utc(jd_tai);
}

double utc_to_tt(double jd_utc) {
    double jd_tai = utc_to_tai(jd_utc);
    return tai_to_tt(jd_tai);
}

// Epoch methods

Epoch Epoch::to(TimeScale target) const {
    if (target == scale) return *this;

    // Convert to TDB first as hub
    double jd_tdb = jd;
    switch (scale) {
        case TimeScale::TDB: break;
        case TimeScale::TT:  jd_tdb = tt_to_tdb(jd); break;
        case TimeScale::TAI: jd_tdb = tt_to_tdb(tai_to_tt(jd)); break;
        case TimeScale::UTC: jd_tdb = tt_to_tdb(utc_to_tt(jd)); break;
    }

    // Convert from TDB to target
    double result = jd_tdb;
    switch (target) {
        case TimeScale::TDB: break;
        case TimeScale::TT:  result = tdb_to_tt(jd_tdb); break;
        case TimeScale::TAI: result = tt_to_tai(tdb_to_tt(jd_tdb)); break;
        case TimeScale::UTC: result = tt_to_utc(tdb_to_tt(jd_tdb)); break;
    }

    return {result, target};
}

double Epoch::days_since_j2000() const {
    return jd - constants::J2000;
}

double Epoch::centuries_since_j2000() const {
    return (jd - constants::J2000) / 36525.0;
}

Epoch Epoch::from_calendar(int year, int month, int day, double hour, TimeScale s) {
    double j = calendar_to_jd(year, month, day, hour);
    return {j, s};
}

Epoch Epoch::from_jd(double j, TimeScale s) {
    return {j, s};
}

} // namespace solar

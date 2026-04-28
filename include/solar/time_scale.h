#pragma once

#include <string>

namespace solar {

enum class TimeScale {
    TDB,    // Barycentric Dynamical Time (DE440 native)
    TT,     // Terrestrial Time (TDB - periodic ~1.7ms)
    UTC,    // Coordinated Universal Time (TT - leap seconds - 32.184s)
    TAI,    // International Atomic Time (TT - 32.184s)
};

struct Epoch {
    double jd;
    TimeScale scale;

    // Convert to another time scale
    Epoch to(TimeScale target) const;

    // Days and centuries since J2000.0
    double days_since_j2000() const;
    double centuries_since_j2000() const;

    // Construct from calendar date
    static Epoch from_calendar(int year, int month, int day,
                               double hour = 0.0,
                               TimeScale scale = TimeScale::UTC);
    static Epoch from_jd(double jd, TimeScale scale = TimeScale::TDB);
};

const char* time_scale_name(TimeScale ts);

// Core conversions (Julian dates)
double tdb_to_tt(double jd_tdb);
double tt_to_tdb(double jd_tt);
double tt_to_tai(double jd_tt);
double tai_to_tt(double jd_tai);
double tai_to_utc(double jd_tai);
double utc_to_tai(double jd_utc);
double tt_to_utc(double jd_tt);
double utc_to_tt(double jd_utc);

// Leap seconds accumulated at a given UTC Julian date
int leap_seconds_at(double jd_utc);

} // namespace solar

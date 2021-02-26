import Foundation

/// Convert to degrees from radians
func toDegrees<T: FloatingPoint>(_ rad: T) -> T {
    rad * 180 / T.pi
}

/// Convert to radians from degrees
func toRadians<T: FloatingPoint>(_ deg: T) -> T {
    deg * T.pi / 180
}

/// Astronomical Unit in km. As defined by JPL
let AU: Double = 149_597_870.691

/// Earth equatorial radius in km. IERS 2003 Conventions
let EARTH_RADIUS: Double = 6378.1366

/// Length of a sidereal day in days according to IERS Conventions
let SIDEREAL_DAY_LENGTH: Double = 1.00273781191135448

/// Julian century conversion constant = 100 * days per year
let JULIAN_DAYS_PER_CENTURY: Double = 36525

/// Seconds in one day
let SECONDS_PER_DAY: Double = 86400

/// Our default epoch.
/// The Julian Day which represents noon on 2000-01-01
let J2000: Double = 2_451_545

/// Lunar cycle length in days
let LUNAR_CYCLE_DAYS: Double = 29.530588853

/// The set of twilights to calculate (types of rise/set events)
public enum Twilight {
    /// Event ID for calculation of rising and setting times for astronomical
    /// twilight. In this case, the calculated time will be the time when the
    /// center of the object is at -18 degrees of geometrical elevation below the
    /// astronomical horizon. At this time astronomical observations are possible
    /// because the sky is dark enough
    case Astronomical

    /// Event ID for calculation of rising and setting times for nautical
    /// twilight. In this case, the calculated time will be the time when the
    /// center of the object is at -12 degrees of geometric elevation below the
    /// astronomical horizon
    case Nautical

    /// Event ID for calculation of rising and setting times for civil twilight.
    /// In this case, the calculated time will be the time when the center of the
    /// object is at -6 degrees of geometric elevation below the astronomical
    /// horizon
    case Civil

    /// The standard value of 34' for the refraction at the local horizon
    case Horizon34arcmin
}

/// Error types which can be thrown
public enum Errors: Error {
    /// Julian day does not exists
    case invalidJulianDay(_ jd: Double)

    /// Invalid location
    case invalidLocation(longitude: Double, latitude: Double)
}

/// Create instance of Sun/Moon Calculator
/// - Parameters:
///   - date: The date/time of observations (local timezone)
///   - longitude: Longitude of observation (degrees)
///   - latitude: Latitude of observation (degrees)
///   - twilight: twilight configuration
/// - Throws: invalidJulianDay if the date does not exists
/// - Returns: Sun and moon information
public func calcSunAndMoon(date: Date, latitude: Double, longitude: Double, twilight: Twilight = .Horizon34arcmin) throws -> (sun: Sun, moon: Moon) {
    let c = try SunMoonCalculator(date: date, longitude: longitude, latitude: latitude, twilight: twilight)
    c.calcSunAndMoon()

    return (
        sun: Sun(ephemeris: Ephemeris(azimuth: Measurement(value: c.sunAzimuth, unit: .radians), elevation: Measurement(value: c.sunElevation, unit: .radians), rise: Date(julianDay: c.sunRise)!, set: Date(julianDay: c.sunSet)!, transit: Date(julianDay: c.sunTransit)!, transitElevation: Measurement(value: c.sunTransitElevation, unit: .radians), distance: Measurement(value: c.sunDistance, unit: .astronomicalUnits))),
        moon: Moon(ephemeris: Ephemeris(azimuth: Measurement(value: c.moonAzimuth, unit: .radians), elevation: Measurement(value: c.moonElevation, unit: .radians), rise: Date(julianDay: c.moonRise)!, set: Date(julianDay: c.moonSet)!, transit: Date(julianDay: c.moonTransit)!, transitElevation: Measurement(value: c.moonTransitElevation, unit: .radians), distance: Measurement(value: c.moonDistance, unit: .astronomicalUnits)), age: c.moonAge, illumination: c.moonIllumination, positionAngleOfAxis: Measurement(value: c.moonP, unit: .radians), brightLimbAngle: Measurement(value: c.moonBL, unit: .radians), paralacticAngle: Measurement(value: c.moonPar, unit: .radians))
    )
}

public struct Sun: Equatable {
    let ephemeris: Ephemeris
}

public struct Moon: Equatable {
    let ephemeris: Ephemeris

    /// Age (days: 0-29.5)
    let age: Double

    /// Illumination (percentage)
    let illumination: Double

    /// Position angle of axis (radians)
    let positionAngleOfAxis: Measurement<UnitAngle>

    /// Bright limb angle (radians)
    let brightLimbAngle: Measurement<UnitAngle>

    /// Paralactic angle (radians)
    let paralacticAngle: Measurement<UnitAngle>
}

public struct Ephemeris: Equatable {
    /// Azimuth (radians)
    let azimuth: Measurement<UnitAngle>

    /// Elevation (radians)
    let elevation: Measurement<UnitAngle>

    /// Rise (Julian days as per UTC)
    let rise: Date

    /// Set (Date per UTC)
    let set: Date

    /// Transit (Date per UTC)
    let transit: Date

    /// Transit elevation (radians)
    let transitElevation: Measurement<UnitAngle>

    /// Sun distance (AUs)
    let distance: Measurement<UnitLength>
}

/// A very simple Sun/Moon calculator without using JPARSEC library
/// - note: Swift port of the excellent [ephemerides (in Java)](http://conga.oan.es/~alonso/doku.php?id=blog:sun_moon_position) by Tomás Alonso Albi
class SunMoonCalculator {
    /// Create instance of Sun/Moon Calculator
    /// - Parameters:
    ///   - date: The date/time of observations (local timezone)
    ///   - longitude: Longitude of observation (degrees)
    ///   - latitude: Latitude of observation (degrees)
    ///   - twilight: Twilight type
    /// - Throws: invalidJulianDay if the date does not exists
    init(date: Date, longitude: Double, latitude: Double, twilight: Twilight = .Horizon34arcmin) throws {
        if longitude.isNaN || latitude.isNaN || abs(longitude) > 180 || abs(latitude) > 90 {
            throw Errors.invalidLocation(longitude: longitude, latitude: latitude)
        }
        var calendar = Calendar(identifier: .gregorian)
        calendar.timeZone = TimeZone(abbreviation: "UTC")!
        let dc: DateComponents = calendar.dateComponents([.year, .month, .day, .hour, .minute, .second], from: date),
            year: Int = dc.year!,
            month: Int = dc.month!,
            day: Int = dc.day!,
            h: Int = dc.hour!,
            m: Int = dc.minute!,
            s: Int = dc.second!

        // The conversion formulas are from Meeus, chapter 7.
        var julian: Bool = false
        if year < 1582 || (year == 1582 && month <= 10) || (year == 1582 && month == 10 && day < 15) {
            julian = true
        }
        let D: Int = day
        var M: Int = month,
            Y: Int = year
        if M < 3 {
            Y -= 1
            M += 12
        }
        let A: Int = Y / 100,
            B: Int = julian ? 0 : 2 - A + A / 4,
            dayFraction: Double = (Double(h) + (Double(m) + (Double(s) / 60)) / 60) / 24,
            jd: Double = dayFraction + Double(Int(365.25 * Double(Y + 4716)) + Int(30.6001 * Double(M + 1))) + Double(D + B) - 1524.5

        try validateJulianDay(jd)

        TTminusUT = 0
        if year > -600, year < 2200 {
            let x = Double(year) + (Double(month) - 1 + Double(day) / 30) / 12
            let x2: Double = x * x, x3: Double = x2 * x, x4: Double = x3 * x
            if year < 1600 {
                TTminusUT = 10535.328003326353 - 9.995238627481024 * x + 0.003067307630020489 * x2 - 7.76340698361363e-6 * x3 + 3.1331045394223196e-9 * x4 +
                    8.225530854405553e-12 * x2 * x3 - 7.486164715632051e-15 * x4 * x2 + 1.9362461549678834e-18 * x4 * x3 - 8.489224937827653e-23 * x4 * x4
            } else {
                TTminusUT = -1_027_175.3477559977 + 2523.256625418965 * x - 1.885686849058459 * x2 + 5.869246227888417e-5 * x3 + 3.3379295816475025e-7 * x4 +
                    1.7758961671447929e-10 * x2 * x3 - 2.7889902806153024e-13 * x2 * x4 + 1.0224295822336825e-16 * x3 * x4 - 1.2528102370680435e-20 * x4 * x4
            }
        }
        obsLon = toRadians(longitude)
        obsLat = toRadians(latitude)
        self.twilight = twilight
        setUTDate(jd)
    }

    /// Calculates everything for the Sun and the Moon
    func calcSunAndMoon() {
        let jd: Double = jd_UT

        // First the Sun
        var out: [Double] = doCalc(getSun())
        sunAzimuth = out[0]
        sunElevation = out[1]
        sunRise = out[2]
        sunSet = out[3]
        sunTransit = out[4]
        sunTransitElevation = out[5]
        let sunRA: Double = out[6], sunDec: Double = out[7]
        sunDistance = out[8]
        let sa: Double = sanomaly, sl: Double = slongitude, lst: Double = out[9]

        var niter: Int = 3 // Number of iterations to get accurate rise/set/transit times
        sunRise = obtainAccurateRiseSetTransit(riseSetJD: sunRise, index: 2, niter: niter, sun: true)
        sunSet = obtainAccurateRiseSetTransit(riseSetJD: sunSet, index: 3, niter: niter, sun: true)
        sunTransit = obtainAccurateRiseSetTransit(riseSetJD: sunTransit, index: 4, niter: niter, sun: true)
        if sunTransit == -1 {
            sunTransitElevation = 0
        } else {
            // Update Sun's maximum elevation
            setUTDate(sunTransit)
            out = doCalc(getSun())
            sunTransitElevation = out[5]
        }

        // Now Moon
        setUTDate(jd)
        sanomaly = sa
        slongitude = sl
        out = doCalc(getMoon())
        moonAzimuth = out[0]
        moonElevation = out[1]
        moonRise = out[2]
        moonSet = out[3]
        moonTransit = out[4]
        moonTransitElevation = out[5]
        let moonRA = out[6],
            moonDec = out[7]
        moonDistance = out[8]
        moonIllumination = (1 - cos(acos(sin(sunDec) * sin(moonDec) + cos(sunDec) * cos(moonDec) * cos(moonRA - sunRA)))) * 0.5
        let ma: Double = moonAge

        niter = 5 // Number of iterations to get accurate rise/set/transit times
        moonRise = obtainAccurateRiseSetTransit(riseSetJD: moonRise, index: 2, niter: niter, sun: false)
        moonSet = obtainAccurateRiseSetTransit(riseSetJD: moonSet, index: 3, niter: niter, sun: false)
        moonTransit = obtainAccurateRiseSetTransit(riseSetJD: moonTransit, index: 4, niter: niter, sun: false)
        if moonTransit == -1 {
            moonTransitElevation = 0
        } else {
            // Update Moon's maximum elevation
            setUTDate(moonTransit)
            _ = getSun()
            out = doCalc(getMoon())
            moonTransitElevation = out[5]
        }
        setUTDate(jd)
        sanomaly = sa
        slongitude = sl
        moonAge = ma

        out = getMoonDiskOrientationAngles(lst: lst, sunRA: sunRA, sunDec: sunDec,
                                           moonLon: toRadians(moonAzimuth), moonLat: toRadians(moonElevation), moonRA: moonRA, moonDec: moonDec)
        moonP = out[2]
        moonBL = out[3]
        moonPar = out[4]
    }

    private func setUTDate(_ jd: Double) {
        jd_UT = jd
        t = (jd + TTminusUT / SECONDS_PER_DAY - J2000) / JULIAN_DAYS_PER_CENTURY
    }

    private func getSun() -> [Double] {
        // SUN PARAMETERS (Formulae from "Calendrical Calculations")
        let lon: Double = (280.46645 + 36000.76983 * t + 0.0003032 * t * t),
            anom: Double = (357.5291 + 35999.0503 * t - 0.0001559 * t * t - 4.8e-07 * t * t * t)
        sanomaly = toRadians(anom)
        var c: Double = (1.9146 - 0.004817 * t - 0.000014 * t * t) * sin(sanomaly)
        c = c + (0.019993 - 0.000101 * t) * sin(2 * sanomaly)
        c = c + 0.00029 * sin(3 * sanomaly) // Correction to the mean ecliptic longitude
        // Now, let calculate nutation and aberration
        let M1: Double = toRadians(124.90 - 1934.134 * t + 0.002063 * t * t),
            M2: Double = toRadians(201.11 + 72001.5377 * t + 0.00057 * t * t),
            d: Double = -0.00569 - 0.0047785 * sin(M1) - 0.0003667 * sin(M2)

        slongitude = lon + c + d // apparent longitude (error<0.003 deg)
        let slatitude: Double = 0, // Sun's ecliptic latitude is always negligible
            ecc: Double = 0.016708617 - 4.2037e-05 * t - 1.236e-07 * t * t, // Eccentricity
            v: Double = sanomaly + toRadians(c), // True anomaly
            sdistance: Double = 1.000001018 * (1 - ecc * ecc) / (1 + ecc * cos(v)) // In UA
        return [slongitude, slatitude, sdistance, atan(696_000 / (AU * sdistance))]
    }

    private func getMoon() -> [Double] {
        // MOON PARAMETERS (Formulae from "Calendrical Calculations")
        let phase: Double = normalizeRadians(toRadians(297.8502042 + 445_267.1115168 * t - 0.00163 * t * t + t * t * t / 538_841 - t * t * t * t / 65_194_000))

        // Anomalistic phase
        var anomaly: Double = (134.9634114 + 477_198.8676313 * t + 0.008997 * t * t + t * t * t / 69699 - t * t * t * t / 14_712_000)
        anomaly = toRadians(anomaly)

        // Degrees from ascending node
        var node: Double = (93.2720993 + 483_202.0175273 * t - 0.0034029 * t * t - t * t * t / 3_526_000 + t * t * t * t / 863_310_000)
        node = toRadians(node)

        let E: Double = 1 - (0.002495 + 7.52e-06 * (t + 1)) * (t + 1)

        // Now longitude, with the three main correcting terms of evection,
        // variation, and equation of year, plus other terms (error<0.01 deg)
        // P. Duffet's MOON program taken as reference
        var l: Double = (218.31664563 + 481_267.8811958 * t - 0.00146639 * t * t + t * t * t / 540_135.03 - t * t * t * t / 65_193_770.4)
        l += 6.28875 * sin(anomaly) + 1.274018 * sin(2 * phase - anomaly) + 0.658309 * sin(2 * phase)
        l += 0.213616 * sin(2 * anomaly) - E * 0.185596 * sin(sanomaly) - 0.114336 * sin(2 * node)
        l += 0.058793 * sin(2 * phase - 2 * anomaly) + 0.057212 * E * sin(2 * phase - anomaly - sanomaly) + 0.05332 * sin(2 * phase + anomaly)
        l += 0.045874 * E * sin(2 * phase - sanomaly) + 0.041024 * E * sin(anomaly - sanomaly) - 0.034718 * sin(phase) - E * 0.030465 * sin(sanomaly + anomaly)
        l += 0.015326 * sin(2 * (phase - node)) - 0.012528 * sin(2 * node + anomaly) - 0.01098 * sin(2 * node - anomaly) + 0.010674 * sin(4 * phase - anomaly)
        l += 0.010034 * sin(3 * anomaly) + 0.008548 * sin(4 * phase - 2 * anomaly)
        l += -E * 0.00791 * sin(sanomaly - anomaly + 2 * phase) - E * 0.006783 * sin(2 * phase + sanomaly) + 0.005162 * sin(anomaly - phase) + E * 0.005 * sin(sanomaly + phase)
        l += 0.003862 * sin(4 * phase) + E * 0.004049 * sin(anomaly - sanomaly + 2 * phase) + 0.003996 * sin(2 * (anomaly + phase)) + 0.003665 * sin(2 * phase - 3 * anomaly)
        l += E * 2.695e-3 * sin(2 * anomaly - sanomaly) + 2.602e-3 * sin(anomaly - 2 * (node + phase))
        l += E * 2.396e-3 * sin(2 * (phase - anomaly) - sanomaly) - 2.349e-3 * sin(anomaly + phase)
        l += E * E * 2.249e-3 * sin(2 * (phase - sanomaly)) - E * 2.125e-3 * sin(2 * anomaly + sanomaly)
        l += -E * E * 2.079e-3 * sin(2 * sanomaly) + E * E * 2.059e-3 * sin(2 * (phase - sanomaly) - anomaly)
        l += -1.773e-3 * sin(anomaly + 2 * (phase - node)) - 1.595e-3 * sin(2 * (node + phase))
        l += E * 1.22e-3 * sin(4 * phase - sanomaly - anomaly) - 1.11e-3 * sin(2 * (anomaly + node))
        var longitude: Double = l

        // Let's add nutation here also
        let M1: Double = toRadians(124.90 - 1934.134 * t + 0.002063 * t * t),
            M2: Double = toRadians(201.11 + 72001.5377 * t + 0.00057 * t * t),
            d: Double = -0.0047785 * sin(M1) - 0.0003667 * sin(M2)
        longitude += d

        // Get accurate Moon age
        moonAge = normalizeRadians(toRadians(longitude - slongitude)) * LUNAR_CYCLE_DAYS / (2 * Double.pi)

        // Now Moon parallax
        var parallax: Double = 0.950724 + 0.051818 * cos(anomaly) + 0.009531 * cos(2 * phase - anomaly)
        parallax += 0.007843 * cos(2 * phase) + 0.002824 * cos(2 * anomaly)
        parallax += 0.000857 * cos(2 * phase + anomaly) + E * 0.000533 * cos(2 * phase - sanomaly)
        parallax += E * 0.000401 * cos(2 * phase - anomaly - sanomaly) + E * 0.00032 * cos(anomaly - sanomaly) - 0.000271 * cos(phase)
        parallax += -E * 0.000264 * cos(sanomaly + anomaly) - 0.000198 * cos(2 * node - anomaly)
        parallax += 1.73e-4 * cos(3 * anomaly) + 1.67e-4 * cos(4 * phase - anomaly)

        // So Moon distance in Earth radii is, more or less,
        let distance: Double = 1 / sin(toRadians(parallax))

        // Ecliptic latitude with nodal phase (error<0.01 deg)
        l = 5.128189 * sin(node) + 0.280606 * sin(node + anomaly) + 0.277693 * sin(anomaly - node)
        l += 0.173238 * sin(2 * phase - node) + 0.055413 * sin(2 * phase + node - anomaly)
        l += 0.046272 * sin(2 * phase - node - anomaly) + 0.032573 * sin(2 * phase + node)
        l += 0.017198 * sin(2 * anomaly + node) + 0.009267 * sin(2 * phase + anomaly - node)
        l += 0.008823 * sin(2 * anomaly - node) + E * 0.008247 * sin(2 * phase - sanomaly - node) + 0.004323 * sin(2 * (phase - anomaly) - node)
        l += 0.0042 * sin(2 * phase + node + anomaly) + E * 0.003372 * sin(node - sanomaly - 2 * phase)
        l += E * 2.472e-3 * sin(2 * phase + node - sanomaly - anomaly)
        l += E * 2.222e-3 * sin(2 * phase + node - sanomaly)
        l += E * 2.072e-3 * sin(2 * phase - node - sanomaly - anomaly)
        let latitude: Double = l

        return [longitude, latitude, distance * EARTH_RADIUS / AU, atan(1737.4 / (distance * EARTH_RADIUS))]
    }

    private func doCalc(_ pos: [Double]) -> [Double] {
        var pos: [Double] = pos
        // Ecliptic to equatorial coordinates
        let t2: Double = t / 100
        var tmp: Double = t2 * (27.87 + t2 * (5.79 + t2 * 2.45))
        tmp = t2 * (-249.67 + t2 * (-39.05 + t2 * (7.12 + tmp)))
        tmp = t2 * (-1.55 + t2 * (1999.25 + t2 * (-51.38 + tmp)))
        tmp = (t2 * (-4680.93 + tmp)) / 3600
        var angle: Double = toRadians(23.4392911111111 + tmp) // obliquity
        // Add nutation in obliquity
        let M1: Double = toRadians(124.90 - 1934.134 * t + 0.002063 * t * t),
            M2: Double = toRadians(201.11 + 72001.5377 * t + 0.00057 * t * t),
            d: Double = 0.002558 * cos(M1) - 0.00015339 * cos(M2)
        angle += toRadians(d)

        pos[0] = toRadians(pos[0])
        pos[1] = toRadians(pos[1])
        let cl: Double = cos(pos[1]),
            x: Double = pos[2] * cos(pos[0]) * cl
        var y: Double = pos[2] * sin(pos[0]) * cl,
            z: Double = pos[2] * sin(pos[1])
        tmp = y * cos(angle) - z * sin(angle)
        z = y * sin(angle) + z * cos(angle)
        y = tmp

        // Obtain local apparent sidereal time
        let jd0: Double = floor(jd_UT - 0.5) + 0.5,
            T0: Double = (jd0 - J2000) / JULIAN_DAYS_PER_CENTURY,
            secs: Double = (jd_UT - jd0) * SECONDS_PER_DAY
        var gmst: Double = (((((-6.2e-6 * T0) + 9.3104e-2) * T0) + 8_640_184.812866) * T0) + 24110.54841
        let msday: Double = 1 + (((((-1.86e-5 * T0) + 0.186208) * T0) + 8_640_184.812866) / (SECONDS_PER_DAY * JULIAN_DAYS_PER_CENTURY))
        gmst = (gmst + msday * secs) * toRadians(15 / 3600)
        let lst: Double = gmst + obsLon

        // Obtain topocentric rectangular coordinates
        // Set radiusAU = 0 for geocentric calculations
        // (rise/set/transit will have no sense in this case)
        let radiusAU: Double = EARTH_RADIUS / AU
        let correction: [Double] = [
            radiusAU * cos(obsLat) * cos(lst),
            radiusAU * cos(obsLat) * sin(lst),
            radiusAU * sin(obsLat),
        ]
        let xtopo: Double = x - correction[0],
            ytopo: Double = y - correction[1],
            ztopo: Double = z - correction[2]

        // Obtain topocentric equatorial coordinates
        var ra: Double = 0,
            dec: Double = Double.pi / 2
        if ztopo < 0 {
            dec = -dec
        }
        if ytopo != 0 || xtopo != 0 {
            ra = atan2(ytopo, xtopo)
            dec = atan2(ztopo / sqrt(xtopo * xtopo + ytopo * ytopo), 1)
        }
        let dist: Double = sqrt(xtopo * xtopo + ytopo * ytopo + ztopo * ztopo)

        // Hour angle
        let angh: Double = lst - ra

        // Obtain azimuth and geometric alt
        let sinlat: Double = sin(obsLat),
            coslat: Double = cos(obsLat),
            sindec: Double = sin(dec), cosdec: Double = cos(dec),
            h: Double = sinlat * sindec + coslat * cosdec * cos(angh)
        var alt: Double = asin(h)
        let azy: Double = sin(angh),
            azx: Double = cos(angh) * sinlat - sindec * coslat / cosdec,
            azi: Double = Double.pi + atan2(azy, azx) // 0 = north
        // Get apparent elevation
        if alt > toRadians(-3) {
            let r: Double = toRadians(0.016667) * abs(tan(Double.pi / 2 - toRadians(toDegrees(alt) + 7.31 / (toDegrees(alt) + 4.4))))
            let refr: Double = r * (0.28 * 1010 / (10 + 273)) // Assuming pressure of 1010 mb and T = 10 C
            alt = min(alt + refr, Double.pi / 2) // This is not accurate, but acceptable
        }

        switch twilight {
        case .Horizon34arcmin:
            // Rise, set, transit times, taking into account Sun/Moon angular radius (pos[3]).
            // The 34' factor is the standard refraction at horizon.
            // Removing angular radius will do calculations for the center of the disk instead
            // of the upper limb.
            tmp = -toRadians(34 / 60) - pos[3]
        case .Civil:
            tmp = toRadians(-6)
        case .Nautical:
            tmp = toRadians(-12)
        case .Astronomical:
            tmp = toRadians(-18)
        }

        // Compute cosine of hour angle
        tmp = (sin(tmp) - sin(obsLat) * sin(dec)) / (cos(obsLat) * cos(dec))
        let celestialHoursToEarthTime: Double = 180 / (15 * Double.pi) / 24 / SIDEREAL_DAY_LENGTH

        // Make calculations for the meridian
        let transit_time1: Double = celestialHoursToEarthTime * normalizeRadians(ra - lst),
            transit_time2: Double = celestialHoursToEarthTime * (normalizeRadians(ra - lst) - 2 * Double.pi)
        var transit_alt: Double = asin(sin(dec) * sin(obsLat) + cos(dec) * cos(obsLat))
        if transit_alt > toRadians(-3) {
            let r: Double = toRadians(0.016667) * abs(tan(Double.pi / 2 - (toDegrees(transit_alt) + 7.31 / toRadians(toDegrees(transit_alt) + 4.4)))),
                refr: Double = r * (0.28 * 1010 / (10 + 273)) // Assuming pressure of 1010 mb and T = 10 C
            transit_alt = min(transit_alt + refr, Double.pi / 2) // This is not accurate, but acceptable
        }

        // Obtain the current event in time
        var transit_time: Double = transit_time1
        let jdToday: Double = floor(jd_UT - 0.5) + 0.5,
            transitToday2: Double = floor(jd_UT + transit_time2 - 0.5) + 0.5
        // Obtain the transit time. Preference should be given to the closest event
        // in time to the current calculation time
        if jdToday == transitToday2, abs(transit_time2) < abs(transit_time1) {
            transit_time = transit_time2
        }
        let transit: Double = jd_UT + transit_time

        // Make calculations for rise and set
        var rise: Double = -1, set: Double = -1
        if abs(tmp) <= 1 {
            let ang_hor: Double = abs(acos(tmp)),
                rise_time1: Double = celestialHoursToEarthTime * normalizeRadians(ra - ang_hor - lst),
                set_time1: Double = celestialHoursToEarthTime * normalizeRadians(ra + ang_hor - lst),
                rise_time2: Double = celestialHoursToEarthTime * (normalizeRadians(ra - ang_hor - lst) - 2 * Double.pi),
                set_time2: Double = celestialHoursToEarthTime * (normalizeRadians(ra + ang_hor - lst) - 2 * Double.pi)

            // Obtain the current events in time. Preference should be given to the closest event
            // in time to the current calculation time (so that iteration in other method will converge)
            var rise_time: Double = rise_time1
            let riseToday2: Double = floor(jd_UT + rise_time2 - 0.5) + 0.5
            if jdToday == riseToday2, abs(rise_time2) < abs(rise_time1) {
                rise_time = rise_time2
            }

            var set_time: Double = set_time1
            let setToday2: Double = floor(jd_UT + set_time2 - 0.5) + 0.5
            if jdToday == setToday2, abs(set_time2) < abs(set_time1) {
                set_time = set_time2
            }
            rise = jd_UT + rise_time
            set = jd_UT + set_time
        }

        return [azi, alt, rise, set, transit, transit_alt, ra, dec, dist, lst]
    }

    private func obtainAccurateRiseSetTransit(riseSetJD: Double, index: Int, niter: Int, sun: Bool) -> Double {
        var riseSetJD: Double = riseSetJD,
            step: Double = -1
        for _ in 0 ..< niter {
            if riseSetJD == -1 {
                return riseSetJD // -1 means no rise/set from that location
            }
            setUTDate(riseSetJD)
            var out: [Double]
            if sun {
                out = doCalc(getSun())
            } else {
                _ = getSun()
                out = doCalc(getMoon())
            }
            step = abs(riseSetJD - out[index])
            riseSetJD = out[index]
        }
        if step > 1 / SECONDS_PER_DAY {
            return -1 // did not converge => without rise/set/transit in this date
        }
        return riseSetJD
    }

    /// Method to calculate values of Moon Disk
    /// - Returns: [optical librations (lp), lunar coordinates of the centre of the disk (bp), position angle of axis (p), bright limb angle (bl), paralactic angle (par)]
    private func getMoonDiskOrientationAngles(lst: Double, sunRA: Double, sunDec: Double,
                                              moonLon: Double, moonLat: Double, moonRA: Double, moonDec: Double) -> [Double] {
        // Moon's argument of latitude
        let F: Double = toRadians(93.2720993 + 483_202.0175273 * t - 0.0034029 * t * t
            - t * t * t / 3_526_000 + t * t * t * t / 863_310_000)
        // Moon's inclination
        let I: Double = toRadians(1.54242)
        // Moon's mean ascending node longitude
        let omega: Double = toRadians(125.0445550 - 1934.1361849 * t + 0.0020762 * t * t
            + t * t * t / 467_410 - t * t * t * t / 18_999_000)
        // Obliquity of ecliptic (approx, better formulae up)
        let eps: Double = toRadians(23.43929)

        // Obtain optical librations lp and bp
        let W: Double = moonLon - omega,
            sinA: Double = sin(W) * cos(moonLat) * cos(I) - sin(moonLat) * sin(I),
            cosA: Double = cos(W) * cos(moonLat),
            A: Double = atan2(sinA, cosA),
            lp: Double = normalizeRadians(A - F),
            sinbp: Double = -sin(W) * cos(moonLat) * sin(I) - sin(moonLat) * cos(I),
            bp: Double = asin(sinbp)

        // Obtain position angle of axis p
        var x: Double = sin(I) * sin(omega),
            y: Double = sin(I) * cos(omega) * cos(eps) - cos(I) * sin(eps)
        let w: Double = atan2(x, y),
            sinp: Double = sqrt(x * x + y * y) * cos(moonRA - w) / cos(bp),
            p: Double = asin(sinp)

        // Compute bright limb angle bl
        let bl: Double = (Double.pi + atan2(cos(sunDec) * sin(moonRA - sunRA), cos(sunDec)
                * sin(moonDec) * cos(moonRA - sunRA) - sin(sunDec) * cos(moonDec)))

        // Paralactic angle par
        y = sin(lst - moonRA)
        x = tan(obsLat) * cos(moonDec) - sin(moonDec) * cos(lst - moonRA)
        let par: Double = x != 0 ? atan2(y, x) : (y / abs(y)) * Double.pi / 2

        return [lp, bp, p, bl, par]
    }

    /// Input values
    private var jd_UT: Double = 0,
                t: Double = 0,
                obsLon: Double = 0,
                obsLat: Double = 0,
                TTminusUT: Double = 0,
                twilight: Twilight = .Horizon34arcmin,
                slongitude: Double = 0,
                sanomaly: Double = 0

    /// Sun azimuth (radians)
    private(set) var sunAzimuth: Double = .nan

    /// Sun elevation (radians)
    private(set) var sunElevation: Double = .nan

    /// Sun rise (Julian days as per UTC)
    private(set) var sunRise: Double = .nan

    /// Sun set (Julian days as per UTC)
    private(set) var sunSet: Double = .nan

    /// Sun transit (Julian days as per UTC)
    private(set) var sunTransit: Double = .nan

    /// Sun transit elevation (radians)
    private(set) var sunTransitElevation: Double = .nan

    /// Sun distance (AUs)
    private(set) var sunDistance: Double = .nan

    /// Moon azimuth (radians)
    private(set) var moonAzimuth: Double = .nan

    /// Moon elevation (radians)
    private(set) var moonElevation: Double = .nan

    /// Moon rise (Julian days as per UTC)
    private(set) var moonRise: Double = .nan

    /// Moon set (Julian days as per UTC)
    private(set) var moonSet: Double = .nan

    /// Moon transit (Julian days as per UTC)
    private(set) var moonTransit: Double = .nan

    /// Moon transit elevation (radians)
    private(set) var moonTransitElevation: Double = .nan

    /// Moon distance (AUs)
    private(set) var moonDistance: Double = .nan

    /// Moon age (days: 0-29.5)
    private(set) var moonAge: Double = .nan

    /// Moon illumination (percentage)
    private(set) var moonIllumination: Double = .nan

    /// Moon position angle of axis (radians)
    private(set) var moonP: Double = .nan

    /// Moon bright limb angle (radians)
    private(set) var moonBL: Double = .nan

    /// Moon paralactic angle (radians)
    private(set) var moonPar: Double = .nan
}

///  Transforms a Julian day (rise/set/transit fields) to a Date
/// - Parameter jd: Julian day
/// - Returns: Date for specified Julian day
public extension Date {
    init?(julianDay jd: Double) {
        do {
            try validateJulianDay(jd)
        } catch {
            return nil
        }

        // The conversion formulas are from Meeus, Chapter 7.
        let Z: Double = floor(jd + 0.5),
            F: Double = jd + 0.5 - Z
        var A: Double = Z
        if Z >= 2_299_161 {
            let a = Int((Z - 1_867_216.25) / 36524.25)
            A += 1 + Double(a) - Double(a) / 4
        }
        let B: Double = A + 1524,
            C: Int = Int((B - 122.1) / 365.25),
            D: Int = Int(Double(C) * 365.25),
            E: Int = Int((B - Double(D)) / 30.6001),
            exactDay: Double = F + B - Double(D) - Double(Int(30.6001 * Double(E))),
            day: Int = Int(exactDay),
            month: Int = (E < 14) ? E - 1 : E - 13
        var year: Int = C - 4715
        if month > 2 {
            year -= 1
        }
        let h: Double = ((exactDay - Double(day)) * SECONDS_PER_DAY) / 3600,
            hour: Int = Int(h),
            m: Double = (h - Double(hour)) * 60,
            minute: Int = Int(m),
            second: Int = Int((m - Double(minute)) * 60)

        var calendar = Calendar(identifier: .gregorian)
        calendar.timeZone = TimeZone(abbreviation: "UTC")!
        let date = calendar.date(from: DateComponents(year: year, month: month, day: day, hour: hour, minute: minute, second: second))!
        self.init(timeIntervalSince1970: date.timeIntervalSince1970)
    }
}

/// Reduce an angle in radians to the range (0 - 2 Pi)
/// - Parameter r: Angle in radians
/// - Returns: Reduced angle in radians
private func normalizeRadians(_ r: Double) -> Double {
    var r: Double = r
    if r < 0, r >= -2 * Double.pi {
        return r + 2 * Double.pi
    }
    if r >= 2 * Double.pi, r < 4 * Double.pi {
        return r - 2 * Double.pi
    }
    if r >= 0, r < 2 * Double.pi {
        return r
    }

    r -= 2 * Double.pi * floor(r * 1 / (2 * Double.pi))
    if r < 0 {
        r += 2 * Double.pi
    }

    return r
}

func validateJulianDay(_ jd: Double) throws {
    if jd < 2_299_160, jd >= 2_299_150 {
        throw Errors.invalidJulianDay(jd)
    }
}

import Foundation

/// Original Java code: http://conga.oan.es/~alonso/doku.php?id=blog:sun_moon_position

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

let J1970: Double = 2_440_588

/// Lunar cycle length in days
let LUNAR_CYCLE_DAYS = 29.530588853

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

typealias ObjectLocation = (latitude: Double, longitude: Double, distance: Double, angularRadius: Double)


public struct Location {
    let latitude: Measurement<UnitAngle>
    let longitude: Measurement<UnitAngle>
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

    /// INPUT VARIABLES

    let jd = JulianDate(date: date)
    let obsLon = toRadians(longitude)
    let obsLat = toRadians(latitude)
    
    /// OUTPUT VARIABLES
    
    let sunData = calculateEphemerisData(locationProvider: getSunObjectLocation, niter: 3, jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight)
    let moonData = calculateEphemerisData(locationProvider: getMoonObjectLocation, niter: 5, jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight)

    let moonAge = getMoonAge(jd: jd)
    let moonIllumination = (1 - cos(acos(sin(sunData.declination) * sin(moonData.declination) + cos(sunData.declination) * cos(moonData.declination) * cos(moonData.rightAscension - sunData.rightAscension)))) * 0.5
    let moonVisualAngles = getMoonVisualAngles(moonData: moonData, sunData: sunData, jd: jd, obsLat: obsLat)

    return (
        sun: Sun(ephemeris: Ephemeris(data: sunData)),
        moon: Moon(ephemeris: Ephemeris(data: moonData), age: moonAge, illumination: moonIllumination, visualAngles: moonVisualAngles)
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

    /// Angles for drawing
    let visualAngles: VisualAngles
}

public struct Ephemeris: Equatable {
    /// Azimuth (radians)
    let azimuth: Measurement<UnitAngle>

    /// Elevation (radians)
    let elevation: Measurement<UnitAngle>

    /// Rise (Julian days as per UTC)
    let rise: Date?

    /// Set (Date per UTC)
    let set: Date?

    /// Transit (Date per UTC)
    let transit: Date?

    /// Transit elevation (radians)
    let transitElevation: Measurement<UnitAngle>

    /// Sun distance (AUs)
    let distance: Measurement<UnitLength>
}

public extension Ephemeris {
    static func == (lhs: Ephemeris, rhs: Ephemeris) -> Bool {
        return lhs.azimuth == rhs.azimuth
            && lhs.elevation == rhs.elevation
            && floor(lhs.rise?.timeIntervalSince1970 ?? 0) == floor(rhs.rise?.timeIntervalSince1970 ?? 0)
            && floor(lhs.set?.timeIntervalSince1970 ?? 0) == floor(rhs.set?.timeIntervalSince1970 ?? 0)
            && floor(lhs.transit?.timeIntervalSince1970 ?? 0) == floor(rhs.transit?.timeIntervalSince1970 ?? 0)
            && lhs.transitElevation == rhs.transitElevation
            && lhs.distance == rhs.distance
    }
}

/// Reduce an angle in radians to the range (0 - 2 Pi)
/// - Parameter r: Angle in radians
/// - Returns: Reduced angle in radians
private func normalizeRadians(_ r: Double) -> Double {
    var r: Double = r
    if r < 0, r >= -2 * .pi {
        return r + 2 * .pi
    }
    if r >= 2 * Double.pi, r < 4 * .pi {
        return r - 2 * .pi
    }
    if r >= 0, r < 2 * .pi {
        return r
    }

    r -= 2 * .pi * floor(r * 1 / (2 * .pi))
    if r < 0 {
        r += 2 * .pi
    }

    return r
}

struct EphemerisData {
    
    /// azimuth (radians)
    let azimuth: Double

    /// elevation (radians)
    let elevation: Double

    /// rise (Julian days as per UTC)
    var rise: JulianDate?

    /// set (Julian days as per UTC)
    var set: JulianDate?

    /// transit (Julian days as per UTC)
    var transit: JulianDate?

    /// transit elevation (radians)
    var transitElevation: Double

    /// topographic right ascension  (radians)
    let rightAscension: Double
    
    /// topographic declination  (radians)
    let declination: Double
    
    /// distance (AUs)
    let distance: Double
    
    /// local apparent sidereal time (Julian days as per UTC)
    let localApparentSiderealTime: Double
}


typealias ObjectLocationProvider = (JulianDate) -> ObjectLocation

func doCalc(locationProvider: ObjectLocationProvider, jd: JulianDate, obsLat: Double, obsLon: Double, twilight: Twilight) -> EphemerisData {
    let data = locationProvider(jd)

    let t = jd.timeFactor
    
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

    let lat = toRadians(data.latitude)
    let lon = toRadians(data.longitude)

    let cl: Double = cos(lat),
        x: Double = data.distance * cos(lon) * cl
    var y: Double = data.distance * sin(lon) * cl,
        z: Double = data.distance * sin(lat)
    tmp = y * cos(angle) - z * sin(angle)
    z = y * sin(angle) + z * cos(angle)
    y = tmp

    // Obtain local apparent sidereal time
    let jd0 = floor(Double(jd) - 0.5) + 0.5,
        T0 = (jd0 - J2000) / JULIAN_DAYS_PER_CENTURY,
        secs = (Double(jd) - jd0) * SECONDS_PER_DAY
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
        // Rise, set, transit times, taking into account Sun/Moon angular radius
        // The 34' factor is the standard refraction at horizon.
        // Removing angular radius will do calculations for the center of the disk instead
        // of the upper limb.
        tmp = -toRadians(34 / 60) - data.angularRadius
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
    let transitTime1 = celestialHoursToEarthTime * normalizeRadians(ra - lst),
        transitTime2 = celestialHoursToEarthTime * (normalizeRadians(ra - lst) - 2 * Double.pi)
    var transitAltitude = asin(sin(dec) * sin(obsLat) + cos(dec) * cos(obsLat))
    if transitAltitude > toRadians(-3) {
        let r: Double = toRadians(0.016667) * abs(tan(Double.pi / 2 - (toDegrees(transitAltitude) + 7.31 / toRadians(toDegrees(transitAltitude) + 4.4)))),
            refr: Double = r * (0.28 * 1010 / (10 + 273)) // Assuming pressure of 1010 mb and T = 10 C
        transitAltitude = min(transitAltitude + refr, Double.pi / 2) // This is not accurate, but acceptable
    }

    // Obtain the current event in time
    var transitTime = transitTime1
    let jdToday = floor(Double(jd) - 0.5) + 0.5,
        transitToday2 = floor(Double(jd) + transitTime2 - 0.5) + 0.5
    // Obtain the transit time. Preference should be given to the closest event
    // in time to the current calculation time
    if jdToday == transitToday2, abs(transitTime2) < abs(transitTime1) {
        transitTime = transitTime2
    }
    let transit = JulianDate(Double(jd) + transitTime)

    // Make calculations for rise and set
    var rise: JulianDate?, set: JulianDate?
    if abs(tmp) <= 1 {
        let ang_hor = abs(acos(tmp)),
            riseTime1 = celestialHoursToEarthTime * normalizeRadians(ra - ang_hor - lst),
            setTime1 = celestialHoursToEarthTime * normalizeRadians(ra + ang_hor - lst),
            riseTime2 = celestialHoursToEarthTime * (normalizeRadians(ra - ang_hor - lst) - 2 * Double.pi),
            setTime2 = celestialHoursToEarthTime * (normalizeRadians(ra + ang_hor - lst) - 2 * Double.pi)

        // Obtain the current events in time. Preference should be given to the closest event
        // in time to the current calculation time (so that iteration in other method will converge)
        var riseTime = riseTime1
        let riseToday2 = floor(Double(jd) + riseTime2 - 0.5) + 0.5
        if jdToday == riseToday2, abs(riseTime2) < abs(riseTime1) {
            riseTime = riseTime2
        }

        var setTime = setTime1
        let setToday2 = floor(Double(jd) + setTime2 - 0.5) + 0.5
        if jdToday == setToday2, abs(setTime2) < abs(setTime1) {
            setTime = setTime2
        }
        rise = JulianDate(Double(jd) + riseTime)
        set = JulianDate(Double(jd) + setTime)
    }

    return EphemerisData(azimuth: azi, elevation: alt, rise: rise, set: set, transit: transit, transitElevation: transitAltitude, rightAscension: ra, declination: dec, distance: dist, localApparentSiderealTime: lst)
}

func obtainAccurateRiseSetTransit(data: EphemerisData, keyPath: KeyPath<EphemerisData, JulianDate?>, niter: Int, locationProvider: ObjectLocationProvider, obsLat: Double, obsLon: Double, twilight: Twilight) -> JulianDate? {
    guard var jd = data[keyPath: keyPath] else { return nil } // nil means  no rise/set from that location
    var step: Double = -1
    
    for _ in 0 ..< niter {
        let data = doCalc(locationProvider: locationProvider, jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight)
        guard let newValue = data[keyPath: keyPath] else { return nil }
        step = abs(Double(jd) - Double(newValue))
        jd = newValue
    }
    
    // did not converge => without rise/set/transit in this date
    guard step <= 1 / SECONDS_PER_DAY else { return nil }

    return jd
}

class JulianDate {
    let julianDays: Double
    
    lazy var date: Date = {
        let timeInterval = (julianDays + 0.5 - J1970) * SECONDS_PER_DAY
        return Date(timeIntervalSince1970: timeInterval)
    }()
    
    init(_ julianDays: Double) {
        self.julianDays = julianDays
    }

    init(date: Date) {
        julianDays = date.timeIntervalSince1970 / SECONDS_PER_DAY - 0.5 + J1970
        self.date = date
    }

    private lazy var TTminusUT: Double = {
        var calendar = Calendar(identifier: .gregorian)
        calendar.timeZone = TimeZone(abbreviation: "UTC")!
        let dc = calendar.dateComponents([.year, .month, .day], from: date)
        
        guard let year = dc.year, let month = dc.month, let day = dc.day, year > -600, year < 2200 else { return 0 }
        
        let x = Double(year) + (Double(month) - 1 + Double(day) / 30) / 12
        let x2 = x * x, x3 = x2 * x, x4 = x3 * x
        if year < 1600 {
            return 10535.328003326353 - 9.995238627481024 * x + 0.003067307630020489 * x2 - 7.76340698361363e-6 * x3 + 3.1331045394223196e-9 * x4 +
                8.225530854405553e-12 * x2 * x3 - 7.486164715632051e-15 * x4 * x2 + 1.9362461549678834e-18 * x4 * x3 - 8.489224937827653e-23 * x4 * x4
        }
        
        return -1_027_175.3477559977 + 2523.256625418965 * x - 1.885686849058459 * x2 + 5.869246227888417e-5 * x3 + 3.3379295816475025e-7 * x4 +
                1.7758961671447929e-10 * x2 * x3 - 2.7889902806153024e-13 * x2 * x4 + 1.0224295822336825e-16 * x3 * x4 - 1.2528102370680435e-20 * x4 * x4
    }()

    lazy var timeFactor: Double = {
        return (julianDays + TTminusUT / SECONDS_PER_DAY - J2000) / JULIAN_DAYS_PER_CENTURY
    }()
}

extension Double {
    init(_ jd: JulianDate) {
        self = jd.julianDays
    }
}

// SUN PARAMETERS (Formulae from "Calendrical Calculations")


/// Correction to the mean ecliptic longitude
func getSunLongitudeCorrection(jd: JulianDate) -> Double {
    let t = jd.timeFactor
    let sanomaly = getSunAnomaly(jd: jd)
    var c: Double = (1.9146 - 0.004817 * t - 0.000014 * t * t) * sin(sanomaly)
    c = c + (0.019993 - 0.000101 * t) * sin(2 * sanomaly)
    c = c + 0.00029 * sin(3 * sanomaly)
    
    return c
}

func getSunLongitude(jd: JulianDate) -> Double {
    let t = jd.timeFactor
    let longitudeCorrection = getSunLongitudeCorrection(jd: jd)
    
    // Now, let calculate nutation and aberration
    let M1 = toRadians(124.90 - 1934.134 * t + 0.002063 * t * t),
        M2 = toRadians(201.11 + 72001.5377 * t + 0.00057 * t * t),
        aberration = -0.00569 - 0.0047785 * sin(M1) - 0.0003667 * sin(M2)

    return 280.46645 + 36000.76983 * t + 0.0003032 * t * t + longitudeCorrection + aberration
}

func getSunAnomaly(jd: JulianDate) -> Double {
    let t = jd.timeFactor
    return toRadians(357.5291 + 35999.0503 * t - 0.0001559 * t * t - 4.8e-07 * t * t * t)
}

func getSunObjectLocation(jd: JulianDate) -> ObjectLocation {
    let t = jd.timeFactor
    
    let sanomaly = getSunAnomaly(jd: jd)
    let slongitude = getSunLongitude(jd: jd)
    
    let slatitude: Double = 0, // Sun's ecliptic latitude is always negligible
        ecc: Double = 0.016708617 - 4.2037e-05 * t - 1.236e-07 * t * t, // Eccentricity
        v: Double = sanomaly + toRadians(getSunLongitudeCorrection(jd: jd)), // True anomaly
        sdistance: Double = 1.000001018 * (1 - ecc * ecc) / (1 + ecc * cos(v)) // In UA
    return (slatitude, slongitude, sdistance, atan(696_000 / (AU * sdistance)))
}




func getMoonObjectLocation(jd: JulianDate) -> ObjectLocation {
    let t = jd.timeFactor

    // MOON PARAMETERS (Formulae from "Calendrical Calculations")
    let phase: Double = normalizeRadians(toRadians(297.8502042 + 445_267.1115168 * t - 0.00163 * t * t + t * t * t / 538_841 - t * t * t * t / 65_194_000))

    // Anomalistic phase
    var anomaly: Double = (134.9634114 + 477_198.8676313 * t + 0.008997 * t * t + t * t * t / 69699 - t * t * t * t / 14_712_000)
    anomaly = toRadians(anomaly)

    // Degrees from ascending node
    var node: Double = (93.2720993 + 483_202.0175273 * t - 0.0034029 * t * t - t * t * t / 3_526_000 + t * t * t * t / 863_310_000)
    node = toRadians(node)

    let E: Double = 1 - (0.002495 + 7.52e-06 * (t + 1)) * (t + 1)

    let sanomaly = getSunAnomaly(jd: jd)
    
    // Now longitude, with the three main correcting terms of evection,
    // variation, and equation of year, plus other terms (error<0.01 deg)
    // P. Duffet's MOON program taken as reference
    var longitude = (218.31664563 + 481_267.8811958 * t - 0.00146639 * t * t + t * t * t / 540_135.03 - t * t * t * t / 65_193_770.4)
    longitude += 6.28875 * sin(anomaly) + 1.274018 * sin(2 * phase - anomaly) + 0.658309 * sin(2 * phase)
    longitude += 0.213616 * sin(2 * anomaly) - E * 0.185596 * sin(sanomaly) - 0.114336 * sin(2 * node)
    longitude += 0.058793 * sin(2 * phase - 2 * anomaly) + 0.057212 * E * sin(2 * phase - anomaly - sanomaly) + 0.05332 * sin(2 * phase + anomaly)
    longitude += 0.045874 * E * sin(2 * phase - sanomaly) + 0.041024 * E * sin(anomaly - sanomaly) - 0.034718 * sin(phase) - E * 0.030465 * sin(sanomaly + anomaly)
    longitude += 0.015326 * sin(2 * (phase - node)) - 0.012528 * sin(2 * node + anomaly) - 0.01098 * sin(2 * node - anomaly) + 0.010674 * sin(4 * phase - anomaly)
    longitude += 0.010034 * sin(3 * anomaly) + 0.008548 * sin(4 * phase - 2 * anomaly)
    longitude += -E * 0.00791 * sin(sanomaly - anomaly + 2 * phase) - E * 0.006783 * sin(2 * phase + sanomaly) + 0.005162 * sin(anomaly - phase) + E * 0.005 * sin(sanomaly + phase)
    longitude += 0.003862 * sin(4 * phase) + E * 0.004049 * sin(anomaly - sanomaly + 2 * phase) + 0.003996 * sin(2 * (anomaly + phase)) + 0.003665 * sin(2 * phase - 3 * anomaly)
    longitude += E * 2.695e-3 * sin(2 * anomaly - sanomaly) + 2.602e-3 * sin(anomaly - 2 * (node + phase))
    longitude += E * 2.396e-3 * sin(2 * (phase - anomaly) - sanomaly) - 2.349e-3 * sin(anomaly + phase)
    longitude += E * E * 2.249e-3 * sin(2 * (phase - sanomaly)) - E * 2.125e-3 * sin(2 * anomaly + sanomaly)
    longitude += -E * E * 2.079e-3 * sin(2 * sanomaly) + E * E * 2.059e-3 * sin(2 * (phase - sanomaly) - anomaly)
    longitude += -1.773e-3 * sin(anomaly + 2 * (phase - node)) - 1.595e-3 * sin(2 * (node + phase))
    longitude += E * 1.22e-3 * sin(4 * phase - sanomaly - anomaly) - 1.11e-3 * sin(2 * (anomaly + node))

    // Let's add nutation here also
    let M1: Double = toRadians(124.90 - 1934.134 * t + 0.002063 * t * t),
        M2: Double = toRadians(201.11 + 72001.5377 * t + 0.00057 * t * t),
        d: Double = -0.0047785 * sin(M1) - 0.0003667 * sin(M2)
    longitude += d

    // Now Moon parallax
    var parallax = 0.950724 + 0.051818 * cos(anomaly) + 0.009531 * cos(2 * phase - anomaly)
    parallax += 0.007843 * cos(2 * phase) + 0.002824 * cos(2 * anomaly)
    parallax += 0.000857 * cos(2 * phase + anomaly) + E * 0.000533 * cos(2 * phase - sanomaly)
    parallax += E * 0.000401 * cos(2 * phase - anomaly - sanomaly) + E * 0.00032 * cos(anomaly - sanomaly) - 0.000271 * cos(phase)
    parallax += -E * 0.000264 * cos(sanomaly + anomaly) - 0.000198 * cos(2 * node - anomaly)
    parallax += 1.73e-4 * cos(3 * anomaly) + 1.67e-4 * cos(4 * phase - anomaly)

    // So Moon distance in Earth radii is, more or less,
    let distance: Double = 1 / sin(toRadians(parallax))

    // Ecliptic latitude with nodal phase (error<0.01 deg)
    var latitude = 5.128189 * sin(node) + 0.280606 * sin(node + anomaly) + 0.277693 * sin(anomaly - node)
    latitude += 0.173238 * sin(2 * phase - node) + 0.055413 * sin(2 * phase + node - anomaly)
    latitude += 0.046272 * sin(2 * phase - node - anomaly) + 0.032573 * sin(2 * phase + node)
    latitude += 0.017198 * sin(2 * anomaly + node) + 0.009267 * sin(2 * phase + anomaly - node)
    latitude += 0.008823 * sin(2 * anomaly - node) + E * 0.008247 * sin(2 * phase - sanomaly - node) + 0.004323 * sin(2 * (phase - anomaly) - node)
    latitude += 0.0042 * sin(2 * phase + node + anomaly) + E * 0.003372 * sin(node - sanomaly - 2 * phase)
    latitude += E * 2.472e-3 * sin(2 * phase + node - sanomaly - anomaly)
    latitude += E * 2.222e-3 * sin(2 * phase + node - sanomaly)
    latitude += E * 2.072e-3 * sin(2 * phase - node - sanomaly - anomaly)

    return (latitude, longitude, distance * EARTH_RADIUS / AU, atan(1737.4 / (distance * EARTH_RADIUS)))
}

// Get accurate Moon age
func getMoonAge(jd: JulianDate) -> Double {
    let data = getMoonObjectLocation(jd: jd)
    return normalizeRadians(toRadians(data.longitude - getSunLongitude(jd: jd))) * LUNAR_CYCLE_DAYS / (2 * .pi)
}

struct VisualAngles : Equatable {
    static func == (lhs: VisualAngles, rhs: VisualAngles) -> Bool {
        return lhs.axisPosition == rhs.axisPosition && lhs.brightLimb == rhs.brightLimb && lhs.paralactic == rhs.paralactic && lhs.opticalLibration.bandPass == rhs.opticalLibration.bandPass && lhs.opticalLibration.longPass == rhs.opticalLibration.longPass
    }
    
    /// Position angle of axis (radians)
    let axisPosition: Measurement<UnitAngle>
    
    /// Bright limb angle (radians)
    let brightLimb: Measurement<UnitAngle>
    
    /// Paralactic angle (radians)
    let paralactic: Measurement<UnitAngle>
    
    let opticalLibration: (longPass: Measurement<UnitAngle>, bandPass: Measurement<UnitAngle>)
}

/// Method to calculate values of Moon Disk
/// - Returns: [optical librations (lp), lunar coordinates of the centre of the disk (bp), position angle of axis (p), bright limb angle (bl), paralactic angle (par)]
func getMoonVisualAngles(moonData: EphemerisData, sunData: EphemerisData, jd: JulianDate, obsLat: Double) -> VisualAngles {
    let lst = sunData.localApparentSiderealTime
    let sunRA = sunData.rightAscension
    let sunDec = sunData.declination
    let moonLon = toRadians(moonData.azimuth)
    let moonLat = toRadians(moonData.elevation)
    let moonRA = moonData.rightAscension
    let moonDec = moonData.declination
    
    let t = jd.timeFactor
    
    // Moon's argument of latitude
    let F = toRadians(93.2720993 + 483_202.0175273 * t - 0.0034029 * t * t - t * t * t / 3_526_000 + t * t * t * t / 863_310_000)
    
    // Moon's inclination
    let I = toRadians(1.54242)
    
    // Moon's mean ascending node longitude
    let omega = toRadians(125.0445550 - 1934.1361849 * t + 0.0020762 * t * t + t * t * t / 467_410 - t * t * t * t / 18_999_000)
    
    // Obliquity of ecliptic (approx, better formulae up)
    let eps = toRadians(23.43929)

    // Obtain optical librations lp and bp
    let W = moonLon - omega,
        sinA = sin(W) * cos(moonLat) * cos(I) - sin(moonLat) * sin(I),
        cosA = cos(W) * cos(moonLat),
        A = atan2(sinA, cosA),
        lp = normalizeRadians(A - F),
        sinbp = -sin(W) * cos(moonLat) * sin(I) - sin(moonLat) * cos(I),
        bp = asin(sinbp)

    // Obtain position angle of axis p
    var x = sin(I) * sin(omega),
        y = sin(I) * cos(omega) * cos(eps) - cos(I) * sin(eps)
    let w = atan2(x, y),
        sinp = sqrt(x * x + y * y) * cos(moonRA - w) / cos(bp),
        p = asin(sinp)

    // Compute bright limb angle bl
    let bl = (Double.pi + atan2(cos(sunDec) * sin(moonRA - sunRA), cos(sunDec) * sin(moonDec) * cos(moonRA - sunRA) - sin(sunDec) * cos(moonDec)))

    // Paralactic angle par
    y = sin(lst - moonRA)
    x = tan(obsLat) * cos(moonDec) - sin(moonDec) * cos(lst - moonRA)
    let par: Double = x != 0 ? atan2(y, x) : (y / abs(y)) * Double.pi / 2

    return VisualAngles(
        axisPosition: .init(value: p, unit: .radians), brightLimb: .init(value: bl, unit: .radians), paralactic: .init(value: par, unit: .radians),
        opticalLibration: (.init(value: lp, unit: .radians), .init(value: lp, unit: .radians))
    )
}

extension Double {
    var toRadiansMeasurement: Measurement<UnitAngle> {
        Measurement(value: self, unit: .radians)
    }
}

func calculateEphemerisData(locationProvider: ObjectLocationProvider, niter: Int, jd: JulianDate, obsLat: Double, obsLon: Double, twilight: Twilight) -> EphemerisData {

    var data = doCalc(locationProvider: locationProvider, jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight)

    for keyPath in [\EphemerisData.rise, \EphemerisData.set, \EphemerisData.transit] {
        data[keyPath: keyPath] = obtainAccurateRiseSetTransit(data: data, keyPath: keyPath, niter: niter, locationProvider: locationProvider, obsLat: obsLat, obsLon: obsLon, twilight: twilight)
    }
    
    // Update Sun's maximum elevation
    if let sunTransit = data.transit {
        let calculationData = doCalc(locationProvider: locationProvider, jd: sunTransit, obsLat: obsLat, obsLon: obsLon, twilight: twilight)
        data.transitElevation = calculationData.transitElevation
    } else {
        data.transitElevation = 0
    }
    
    return data
}

extension Ephemeris {
    init(data: EphemerisData) {
        self.init(azimuth: data.azimuth.toRadiansMeasurement, elevation: data.elevation.toRadiansMeasurement, rise: data.rise?.date, set: data.set?.date, transit: data.transit?.date, transitElevation: data.transitElevation.toRadiansMeasurement, distance: Measurement<UnitLength>(value: data.distance, unit: .astronomicalUnits))
    }
}

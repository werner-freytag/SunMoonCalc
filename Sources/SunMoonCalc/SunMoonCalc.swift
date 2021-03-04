import Foundation

/** Radians to degrees. */
let RAD_TO_DEG = 180.0 / .pi

/** Degrees to radians. */
let DEG_TO_RAD = 1.0 / RAD_TO_DEG

/* Arcseconds to radians */
let ARCSEC_TO_RAD = (DEG_TO_RAD / 3600.0)

/// Astronomical Unit in km. As defined by JPL
let AU = 149_597_870.691

/// Earth equatorial radius in km. IERS 2003 Conventions
let EARTH_RADIUS = 6378.1366

/// Two times Pi.
let TWO_PI = 2.0 * .pi

/// Pi divided by two.
let PI_OVER_TWO = Double.pi / 2

/// Length of a sidereal day in days according to IERS Conventions
let SIDEREAL_DAY_LENGTH = 1.00273781191135448

/// Julian century conversion constant = 100 * days per year
let JULIAN_DAYS_PER_CENTURY: Double = 36525

/// Seconds in one day
let SECONDS_PER_DAY: Double = 86400

/// Light time in days for 1 AU. DE405 definition.
let LIGHT_TIME_DAYS_PER_AU = 0.00577551833109

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
    case astronomical

    /// Event ID for calculation of rising and setting times for nautical
    /// twilight. In this case, the calculated time will be the time when the
    /// center of the object is at -12 degrees of geometric elevation below the
    /// astronomical horizon
    case nautical

    /// Event ID for calculation of rising and setting times for civil twilight.
    /// In this case, the calculated time will be the time when the center of the
    /// object is at -6 degrees of geometric elevation below the astronomical
    /// horizon
    case civil

    /// The standard value of 34' for the refraction at the local horizon
    case horizon34arcmin
}

/// The set of bodies to compute ephemerides.
enum Body {
    case mercury
    case venus
    case mars
    case jupiter
    case saturn
    case uranus
    case neptune
    case moon
    case sun
    case emb

    private typealias Data = (index: Int, eqRadius: Double)

    private static let dataMapping: [Self: Data] = [
        .mercury: (0, 2439.7),
        .venus: (1, 6051.8),
        .mars: (3, 3396.19),
        .jupiter: (4, 71492),
        .saturn: (5, 60268),
        .uranus: (6, 25559),
        .neptune: (7, 24764),
        .moon: (-2, 1737.4),
        .sun: (-1, 696_000),
        .emb: (2, 0),
    ]

    var eqRadius: Double { Self.dataMapping[self]!.eqRadius }
}

public struct Location {
    let latitude: Measurement<UnitAngle>
    let longitude: Measurement<UnitAngle>

    init(latitude: Measurement<UnitAngle>, longitude: Measurement<UnitAngle>) {
        self.latitude = latitude
        self.longitude = longitude
    }

    init(_ latitude: Double, _ longitude: Double) {
        self.init(latitude: .init(value: latitude, unit: .degrees), longitude: .init(value: longitude, unit: .degrees))
    }
}

public class Sun {
    /// Create instance of Sun
    /// - Parameters:
    ///   - date: The date of observations
    ///   - location: Location of observation
    ///   - twilight: twilight configuration
    init(location: Location, date: Date = .init(), twilight: Twilight = .horizon34arcmin) {
        calculation = .init(date: date, location: location, twilight: twilight)
    }

    let calculation: SunCalculation

    lazy var ephemeris: Ephemeris = { .init(data: calculation.ephemerisData) }()
}

public class Moon {
    /// Create instance of Moon
    /// - Parameters:
    ///   - date: The date of observations
    ///   - location: Location of observation
    ///   - twilight: twilight configuration
    init(location: Location, date: Date = .init(), twilight: Twilight = .horizon34arcmin) {
        calculation = .init(date: date, location: location, twilight: twilight)
    }

    let calculation: MoonCalculation

    lazy var ephemeris: Ephemeris = { .init(data: calculation.ephemerisData) }()

    /// Age (days: 0-29.5)
    lazy var age: Double = { calculation.age }()

    /// Illumination (percentage)
    lazy var illumination: Double = { calculation.illumination }()
}

public struct Ephemeris {
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

/// Reduce an angle in radians to the range (0 - 2 Pi)
/// - Parameter r: Angle in radians
/// - Returns: Reduced angle in radians
func normalizeRadians(_ r: Double) -> Double {
    switch r {
    case ..<(-TWO_PI):
        return r + TWO_PI * (floor(-r / TWO_PI) + 1)
    case -TWO_PI ..< 0:
        return r + TWO_PI
    case 0 ..< TWO_PI:
        return r
    case TWO_PI ..< (4 * .pi):
        return r - TWO_PI
    default:
        return r - TWO_PI * floor(r / TWO_PI)
    }
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

class JulianDate {
    let jd: Double

    lazy var date: Date = {
        let timeInterval = (jd + 0.5 - J1970) * SECONDS_PER_DAY
        return Date(timeIntervalSince1970: timeInterval)
    }()

    init(_ jd: Double) {
        self.jd = jd
    }

    init(date: Date) {
        jd = date.timeIntervalSince1970 / SECONDS_PER_DAY - 0.5 + J1970
        self.date = date
    }

    private lazy var TTminusUT: Double = {
        var calendar = Calendar(identifier: .gregorian)
        calendar.timeZone = TimeZone(abbreviation: "UTC")!
        let dc = calendar.dateComponents([.year, .month, .day], from: date)

        let year = Double(dc.year!), month = Double(dc.month!), day = Double(dc.day!)

        var TTMinusUT: Double = 0
        let ndot = -25.858
        var c0 = 0.91072 * (ndot + 26.0)
        if year < -500 || year >= 2200 {
            let u = (jd - 2_385_800.5) / 36525.0 // centuries since J1820
            TTminusUT = -20 + 32.0 * u * u
        } else {
            let x = year + (month - 1 + (day - 1) / 30.0) / 12.0
            let x2 = x * x, x3 = x2 * x, x4 = x3 * x
            if year < 1600 {
                TTminusUT = 10535.328003 - 9.9952386275 * x + 0.00306730763 * x2 - 7.7634069836e-6 * x3 + 3.1331045394e-9 * x4 +
                    8.2255308544e-12 * x2 * x3 - 7.4861647156e-15 * x4 * x2 + 1.936246155e-18 * x4 * x3 - 8.4892249378e-23 * x4 * x4
            } else {
                TTminusUT = -1_027_175.34776 + 2523.2566254 * x - 1.8856868491 * x2 + 5.8692462279e-5 * x3 + 3.3379295816e-7 * x4 +
                    1.7758961671e-10 * x2 * x3 - 2.7889902806e-13 * x2 * x4 + 1.0224295822e-16 * x3 * x4 - 1.2528102371e-20 * x4 * x4
            }
            c0 = 0.91072 * (ndot + 25.858)
        }

        let c = -c0 * pow((jd - 2_435_109.0) / 36525.0, 2)
        if year < 1955 || year > 2005 { TTminusUT += c }

        return TTMinusUT
    }()

    lazy var timeFactor: Double = {
        (jd + TTminusUT / SECONDS_PER_DAY - J2000) / JULIAN_DAYS_PER_CENTURY
    }()
}

extension Double {
    init(_ jd: JulianDate) {
        self = jd.jd
    }
}

private typealias ObjectLocation = (latitude: Double, longitude: Double, distance: Double, angularRadius: Double)

class ObjectCalculation {
    /// INPUT VARIABLES
    let jd: JulianDate
    let obsLat: Double
    let obsLon: Double
    let twilight: Twilight

    required init(jd: JulianDate, obsLat: Double, obsLon: Double, twilight: Twilight) {
        self.jd = jd
        self.obsLat = obsLat
        self.obsLon = obsLon
        self.twilight = twilight
    }

    required init(date: Date, location: Location, twilight: Twilight = .horizon34arcmin) {
        jd = JulianDate(date: date)
        obsLon = location.longitude.converted(to: .radians).value
        obsLat = location.latitude.converted(to: .radians).value
        self.twilight = twilight
    }

    fileprivate lazy var t = jd.timeFactor

    private class var accuracyIterationsOfRiseSetTransit: Int { 15 }

    fileprivate var objectLocation: ObjectLocation {
        preconditionFailure("Must be implemented in child classes")
    }

    lazy var ephemerisData: EphemerisData = {
        var ephemerisData = calculateEphemerisData()

        for keyPath in [\EphemerisData.rise, \EphemerisData.set, \EphemerisData.transit] {
            ephemerisData[keyPath: keyPath] = obtainAccurateRiseSetTransit(data: ephemerisData, keyPath: keyPath)
        }

        // Update maximum elevation
        if let transit = ephemerisData.transit {
            let calculationData = Self(jd: transit, obsLat: obsLat, obsLon: obsLon, twilight: twilight).calculateEphemerisData()
            ephemerisData.transitElevation = calculationData.transitElevation
        } else {
            ephemerisData.transitElevation = 0
        }

        return ephemerisData
    }()

    private func calculateEphemerisData() -> EphemerisData {
        let data = objectLocation
        let t = jd.timeFactor

        // Ecliptic to equatorial coordinates
        let t2: Double = t / 100
        var tmp: Double = t2 * (27.87 + t2 * (5.79 + t2 * 2.45))
        tmp = t2 * (-249.67 + t2 * (-39.05 + t2 * (7.12 + tmp)))
        tmp = t2 * (-1.55 + t2 * (1999.25 + t2 * (-51.38 + tmp)))
        tmp = (t2 * (-4680.93 + tmp)) / 3600
        var angle = (23.4392911111111 + tmp) * DEG_TO_RAD // obliquity
        // Add nutation in obliquity
        let M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * DEG_TO_RAD,
            M2 = (201.11 + 72001.5377 * t + 0.00057 * t * t) * DEG_TO_RAD,
            d = 0.002558 * cos(M1) - 0.00015339 * cos(M2)
        angle += d * DEG_TO_RAD

        let lat = data.latitude
        let lon = data.longitude

        let cl = cos(lat),
            x = data.distance * cos(lon) * cl
        var y = data.distance * sin(lon) * cl,
            z = data.distance * sin(lat)
        tmp = y * cos(angle) - z * sin(angle)
        z = y * sin(angle) + z * cos(angle)
        y = tmp

        // Obtain local apparent sidereal time
        let jd0 = floor(Double(jd) - 0.5) + 0.5,
            T0 = (jd0 - J2000) / JULIAN_DAYS_PER_CENTURY,
            secs = (Double(jd) - jd0) * SECONDS_PER_DAY
        var gmst: Double = (((((-6.2e-6 * T0) + 9.3104e-2) * T0) + 8_640_184.812866) * T0) + 24110.54841
        let msday: Double = 1 + (((((-1.86e-5 * T0) + 0.186208) * T0) + 8_640_184.812866) / (SECONDS_PER_DAY * JULIAN_DAYS_PER_CENTURY))
        gmst = (gmst + msday * secs) * (15 / 3600) * DEG_TO_RAD
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
            dec = PI_OVER_TWO
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
        let sinlat = sin(obsLat),
            coslat = cos(obsLat),
            sindec = sin(dec), cosdec: Double = cos(dec),
            h = sinlat * sindec + coslat * cosdec * cos(angh)
        var alt = asin(h)
        let azy = sin(angh),
            azx = cos(angh) * sinlat - sindec * coslat / cosdec,
            azi = .pi + atan2(azy, azx) // 0 = north
        // Get apparent elevation
        if alt > -3 * DEG_TO_RAD {
            let r = 0.016667 * DEG_TO_RAD * abs(tan(PI_OVER_TWO - (alt * RAD_TO_DEG + 7.31 / (alt * RAD_TO_DEG + 4.4)) * DEG_TO_RAD))
            let refr = r * (0.28 * 1010 / (10 + 273)) // Assuming pressure of 1010 mb and T = 10 C
            alt = min(alt + refr, PI_OVER_TWO) // This is not accurate, but acceptable
        }

        switch twilight {
        case .horizon34arcmin:
            // Rise, set, transit times, taking into account Sun/Moon angular radius
            // The 34' factor is the standard refraction at horizon.
            // Removing angular radius will do calculations for the center of the disk instead
            // of the upper limb.
            tmp = -(34 / 60) * DEG_TO_RAD - data.angularRadius
        case .civil:
            tmp = -6 * DEG_TO_RAD
        case .nautical:
            tmp = -12 * DEG_TO_RAD
        case .astronomical:
            tmp = -18 * DEG_TO_RAD
        }

        // Compute cosine of hour angle
        tmp = (sin(tmp) - sin(obsLat) * sin(dec)) / (cos(obsLat) * cos(dec))
        let celestialHoursToEarthTime: Double = 180 / (15 * .pi) / 24 / SIDEREAL_DAY_LENGTH

        // Make calculations for the meridian
        let transitTime1 = celestialHoursToEarthTime * normalizeRadians(ra - lst),
            transitTime2 = celestialHoursToEarthTime * (normalizeRadians(ra - lst) - TWO_PI)
        var transitAltitude = asin(sin(dec) * sin(obsLat) + cos(dec) * cos(obsLat))
        if transitAltitude > -3 * DEG_TO_RAD {
            let r: Double = 0.016667 * DEG_TO_RAD * abs(tan(PI_OVER_TWO - (transitAltitude * RAD_TO_DEG + 7.31 / (transitAltitude * RAD_TO_DEG + 4.4) * DEG_TO_RAD))),
                refr: Double = r * (0.28 * 1010 / (10 + 273)) // Assuming pressure of 1010 mb and T = 10 C
            transitAltitude = min(transitAltitude + refr, PI_OVER_TWO) // This is not accurate, but acceptable
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
                riseTime2 = celestialHoursToEarthTime * (normalizeRadians(ra - ang_hor - lst) - TWO_PI),
                setTime2 = celestialHoursToEarthTime * (normalizeRadians(ra + ang_hor - lst) - TWO_PI)

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

    private func obtainAccurateRiseSetTransit(data: EphemerisData, keyPath: KeyPath<EphemerisData, JulianDate?>) -> JulianDate? {
        guard var jd = data[keyPath: keyPath] else { return nil } // nil means  no rise/set from that location
        var step: Double = -1

        for _ in 0 ..< Self.accuracyIterationsOfRiseSetTransit {
            let ephemerisData = Self(jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight).calculateEphemerisData()
            guard let newValue = ephemerisData[keyPath: keyPath] else { return nil }
            step = abs(Double(jd) - Double(newValue))
            jd = newValue
        }

        // did not converge => without rise/set/transit in this date
        guard step <= 1 / SECONDS_PER_DAY else { return nil }

        return jd
    }
}

// SUN PARAMETERS (Formulae from "Calendrical Calculations")

class SunCalculation: ObjectCalculation {
    /// Correction to the mean ecliptic longitude
    private lazy var longitudeCorrection: Double = {
        var c = (1.9146 - 0.004817 * t - 0.000014 * t * t) * sin(anomaly)
        c = c + (0.019993 - 0.000101 * t) * sin(2 * anomaly)
        c = c + 0.00029 * sin(3 * anomaly)

        return c
    }()

    lazy var longitude: Double = {
        // Now, let calculate nutation and aberration
        let M1 = (124.90 - 1934.134 * t + 0.002063 * t * t) * DEG_TO_RAD,
            M2 = (201.11 + 72001.5377 * t + 0.00057 * t * t) * DEG_TO_RAD,
            aberration = -0.00569 - 0.0047785 * sin(M1) - 0.0003667 * sin(M2)

        return 280.46645 + 36000.76983 * t + 0.0003032 * t * t + longitudeCorrection + aberration
    }()

    lazy var anomaly: Double = {
        (357.5291 + 35999.0503 * t - 0.0001559 * t * t - 4.8e-07 * t * t * t) * DEG_TO_RAD
    }()

    // Sun data from the expansion "Planetary Programs
    // and Tables" by Pierre Bretagnon and Jean-Louis
    // Simon, Willman-Bell, 1986
    private static let data: [[Double]] = [
        [403_406.0, 0.0, 4.721964, 1.621043],
        [195_207.0, -97597.0, 5.937458, 62830.348067],
        [119_433.0, -59715.0, 1.115589, 62830.821524],
        [112_392.0, -56188.0, 5.781616, 62829.634302],
        [3891.0, -1556.0, 5.5474, 125_660.5691],
        [2819.0, -1126.0, 1.512, 125_660.9845],
        [1721.0, -861.0, 4.1897, 62832.4766],
        [0.0, 941.0, 1.163, 0.813],
        [660.0, -264.0, 5.415, 125_659.31],
        [350.0, -163.0, 4.315, 57533.85],
        [334.0, 0.0, 4.553, -33.931],
        [314.0, 309.0, 5.198, 777_137.715],
        [268.0, -158.0, 5.989, 78604.191],
        [242.0, 0.0, 2.911, 5.412],
        [234.0, -54.0, 1.423, 39302.098],
        [158.0, 0.0, 0.061, -34.861],
        [132.0, -93.0, 2.317, 115_067.698],
        [129.0, -20.0, 3.193, 15774.337],
        [114.0, 0.0, 2.828, 5296.67],
        [99.0, -47.0, 0.52, 58849.27],
        [93.0, 0.0, 4.65, 5296.11],
        [86.0, 0.0, 4.35, -3980.7],
        [78.0, -33.0, 2.75, 52237.69],
        [72.0, -32.0, 4.5, 55076.47],
        [68.0, 0.0, 3.23, 261.08],
        [64.0, -10.0, 1.22, 15773.85],
        [46.0, -16.0, 0.14, 188_491.03],
        [38.0, 0.0, 3.44, -7756.55],
        [37.0, 0.0, 4.37, 264.89],
        [32.0, -24.0, 1.14, 117_906.27],
        [29.0, -13.0, 2.84, 55075.75],
        [28.0, 0.0, 5.96, -7961.39],
        [27.0, -9.0, 5.09, 188_489.81],
        [27.0, 0.0, 1.72, 2132.19],
        [25.0, -17.0, 2.56, 109_771.03],
        [24.0, -11.0, 1.92, 54868.56],
        [21.0, 0.0, 0.09, 25443.93],
        [21.0, 31.0, 5.98, -55731.43],
        [20.0, -10.0, 4.03, 60697.74],
        [18.0, 0.0, 4.27, 2132.79],
        [17.0, -12.0, 0.79, 109_771.63],
        [14.0, 0.0, 4.24, -7752.82],
        [13.0, -5.0, 2.01, 188_491.91],
        [13.0, 0.0, 2.65, 207.81],
        [13.0, 0.0, 4.98, 29424.63],
        [12.0, 0.0, 0.93, -7.99],
        [10.0, 0.0, 2.21, 46941.14],
        [10.0, 0.0, 3.59, -68.29],
        [10.0, 0.0, 1.5, 21463.25],
        [10.0, -9.0, 2.55, 157_208.4],
    ]

    override fileprivate var objectLocation: ObjectLocation {
        var L = 0.0, R = 0.0, t2 = t * 0.01
        var Lp = 0.0, deltat = 0.5, t2p = (t + deltat / JULIAN_DAYS_PER_CENTURY) * 0.01
        for elements in Self.data {
            let v = elements[2] + elements[3] * t2
            let u = normalizeRadians(v)
            L = L + elements[0] * sin(u)
            R = R + elements[1] * cos(u)

            let vp = elements[2] + elements[3] * t2p
            let up = normalizeRadians(vp)
            Lp = Lp + elements[0] * sin(up)
        }

        var lon = normalizeRadians(4.9353929 + normalizeRadians(62833.196168 * t2) + L / 10_000_000.0) * RAD_TO_DEG
        let distance = 1.0001026 + R / 10_000_000.0

        // Now subtract aberration
        let dlon = ((Lp - L) / 10_000_000.0 + 62833.196168 * (t2p - t2)) / deltat
        let aberration = dlon * distance * LIGHT_TIME_DAYS_PER_AU
        lon -= aberration * RAD_TO_DEG

        let longitude = lon * DEG_TO_RAD // apparent longitude (error<0.001 deg)
        let latitude: Double = 0 // Sun's ecliptic latitude is always negligible

        return (latitude, longitude, distance, atan(Body.sun.eqRadius / (AU * distance)))
    }
}

class MoonCalculation: ObjectCalculation {
    private lazy var sun = SunCalculation(jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight)

    // Anomalistic phase
    private lazy var anomaly: Double = {
        (134.9634114 + 477_198.8676313 * t + 0.008997 * t * t + t * t * t / 69699 - t * t * t * t / 14_712_000) * DEG_TO_RAD
    }()

    override fileprivate var objectLocation: ObjectLocation {
        // These expansions up to t^7 for the mean elements are taken from S. L. Moshier, see program cmoon
        /* Mean elongation of moon = D */
        var x = (1.6029616009939659e+09 * t + 1.0722612202445078e+06)
        x += (((((-3.207663637426e-013 * t + 2.555243317839e-011) * t + 2.560078201452e-009) * t - 3.702060118571e-005) * t + 6.9492746836058421e-03) * t /* D, t^3 */
            - 6.7352202374457519e+00) * t * t /* D, t^2 */
        let phase = normalizeRadians(ARCSEC_TO_RAD * x)

        /* Mean distance of moon from its ascending node = F */
        x = (1.7395272628437717e+09 * t + 3.3577951412884740e+05)
        x += (((((4.474984866301e-013 * t + 4.189032191814e-011) * t - 2.790392351314e-009) * t - 2.165750777942e-006) * t - 7.5311878482337989e-04) * t /* F, t^3 */
            - 1.3117809789650071e+01) * t * t /* F, t^2 */
        let node = normalizeRadians(ARCSEC_TO_RAD * x)

        /* Mean anomaly of sun = l' (J. Laskar) */
        x = (1.2959658102304320e+08 * t + 1.2871027407441526e+06)
        x += ((((((((1.62e-20 * t - 1.0390e-17) * t - 3.83508e-15) * t + 4.237343e-13) * t + 8.8555011e-11) * t - 4.77258489e-8) * t - 1.1297037031e-5) * t + 8.7473717367324703e-05) * t - 5.5281306421783094e-01) * t * t
        let sanomaly = normalizeRadians(ARCSEC_TO_RAD * x)

        /* Mean anomaly of moon = l */
        x = (1.7179159228846793e+09 * t + 4.8586817465825332e+05)
        x += (((((-1.755312760154e-012 * t + 3.452144225877e-011) * t - 2.506365935364e-008) * t - 2.536291235258e-004) * t + 5.2099641302735818e-02) * t /* l, t^3 */
            + 3.1501359071894147e+01) * t * t /* l, t^2 */
        let anomaly = normalizeRadians(ARCSEC_TO_RAD * x)

        /* Mean longitude of moon, re mean ecliptic and equinox of date = L */
        x = (1.7325643720442266e+09 * t + 7.8593980921052420e+05)
        x += (((((7.200592540556e-014 * t + 2.235210987108e-010) * t - 1.024222633731e-008) * t - 6.073960534117e-005) * t + 6.9017248528380490e-03) * t /* L, t^3 */
            - 5.6550460027471399e+00) * t * t /* L, t^2 */
        var longitude = normalizeRadians(ARCSEC_TO_RAD * x) * RAD_TO_DEG

        // Now longitude, with the three main correcting terms of evection,
        // variation, and equation of year, plus other terms (error<0.01 deg)
        // P. Duffet's MOON program taken as reference for the periodic terms
        let E = 1.0 - (0.002495 + 7.52e-06 * (t + 1.0)) * (t + 1.0)

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

        // Now Moon parallax
        var parallax = 0.950724 + 0.051818 * cos(anomaly) + 0.009531 * cos(2 * phase - anomaly)
        parallax += 0.007843 * cos(2 * phase) + 0.002824 * cos(2 * anomaly)
        parallax += 0.000857 * cos(2 * phase + anomaly) + E * 0.000533 * cos(2 * phase - sanomaly)
        parallax += E * 0.000401 * cos(2 * phase - anomaly - sanomaly) + E * 0.00032 * cos(anomaly - sanomaly) - 0.000271 * cos(phase)
        parallax += -E * 0.000264 * cos(sanomaly + anomaly) - 0.000198 * cos(2 * node - anomaly)
        parallax += 1.73e-4 * cos(3 * anomaly) + 1.67e-4 * cos(4 * phase - anomaly)

        // So Moon distance in Earth radii is, more or less,
        let distance = 1.0 / sin(parallax * DEG_TO_RAD)

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

        return (latitude * DEG_TO_RAD, longitude * DEG_TO_RAD, distance * EARTH_RADIUS / AU, atan(Body.moon.eqRadius / (distance * EARTH_RADIUS)))
    }

    // Get accurate Moon age
    var age: Double {
        normalizeRadians(objectLocation.longitude - sun.longitude) * LUNAR_CYCLE_DAYS / TWO_PI
    }

    lazy var illumination: Double = {
        let sun = sun.ephemerisData
        let moon = ephemerisData

        return (1 - cos(acos(sin(sun.declination) * sin(moon.declination) + cos(sun.declination) * cos(moon.declination) * cos(moon.rightAscension - sun.rightAscension)))) * 0.5
    }()
}

private extension Double {
    var toRadiansMeasurement: Measurement<UnitAngle> {
        Measurement(value: self, unit: .radians)
    }
}

private extension Ephemeris {
    init(data: EphemerisData) {
        self.init(azimuth: data.azimuth.toRadiansMeasurement, elevation: data.elevation.toRadiansMeasurement, rise: data.rise?.date, set: data.set?.date, transit: data.transit?.date, transitElevation: data.transitElevation.toRadiansMeasurement, distance: Measurement<UnitLength>(value: data.distance, unit: .astronomicalUnits))
    }
}

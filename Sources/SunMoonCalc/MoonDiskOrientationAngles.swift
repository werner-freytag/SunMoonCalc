import Foundation

public struct DiskOrientationAngles: Equatable {
    public static func == (lhs: DiskOrientationAngles, rhs: DiskOrientationAngles) -> Bool {
        return lhs.axisPosition == rhs.axisPosition && lhs.brightLimb == rhs.brightLimb && lhs.paralactic == rhs.paralactic && lhs.opticalLibration.bandPass == rhs.opticalLibration.bandPass && lhs.opticalLibration.longPass == rhs.opticalLibration.longPass
    }

    /// Position angle of axis (radians)
    public let axisPosition: Measurement<UnitAngle>

    /// Bright limb angle (radians)
    public let brightLimb: Measurement<UnitAngle>

    /// Paralactic angle (radians)
    public let paralactic: Measurement<UnitAngle>

    public let opticalLibration: (longPass: Measurement<UnitAngle>, bandPass: Measurement<UnitAngle>)
}

/// Calculates the orientation angles of the lunar disk figure. Illumination fraction
/// is returned in the main program. Simplification of the method presented by
/// Eckhardt, D. H., "Theory of the Libration of the Moon", Moon and planets 25, 3
/// (1981), without the physical librations of the Moon. Accuracy around 0.5 deg
/// for each value.
extension MoonCalculation {
    var diskOrientationAngles: DiskOrientationAngles {
        let sunCalculation = SunCalculation(jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight, twilightMode: twilightMode)
        let sunData = sunCalculation.ephemerisData
        let moonData = ephemerisData

        let lst = sunCalculation.localApparentSiderealTime
        let sunRA = sunData.rightAscension
        let sunDec = sunData.declination
        let moonLon = objectLocation.longitude
        let moonLat = objectLocation.latitude
        let moonRA = moonData.rightAscension
        let moonDec = moonData.declination

        let t = jd.t

        // Obliquity of ecliptic
        let eps = meanObliquity + nutation.obl

        // Moon's argument of latitude
        let F = (93.2720993 + 483_202.0175273 * t - 0.0034029 * t * t - t * t * t / 3_526_000.0 + t * t * t * t / 863_310_000.0) * DEG_TO_RAD

        // Moon's inclination
        let I = 1.54242 * DEG_TO_RAD

        // Moon's mean ascending node longitude
        let omega = (125.0445550 - 1934.1361849 * t + 0.0020762 * t * t + t * t * t / 467_410.0 - t * t * t * t / 18_999_000.0) * DEG_TO_RAD

        let cosI = cos(I), sinI = sin(I)
        let cosMoonLat = cos(moonLat), sinMoonLat = sin(moonLat)
        let cosMoonDec = cos(moonDec), sinMoonDec = sin(moonDec)

        // Obtain optical librations lp and bp
        let W = moonLon - omega
        let sinA = sin(W) * cosMoonLat * cosI - sinMoonLat * sinI
        let cosA = cos(W) * cosMoonLat
        let A = atan2(sinA, cosA)
        let lp = normalizeRadians(A - F)
        let sinbp = -sin(W) * cosMoonLat * sinI - sinMoonLat * cosI
        let bp = asin(sinbp)

        // Obtain position angle of axis p
        var x = sinI * sin(omega)
        var y = sinI * cos(omega) * cos(eps) - cosI * sin(eps)
        let w = atan2(x, y)
        let sinp = hypot(x, y) * cos(moonRA - w) / cos(bp)
        let p = asin(sinp)

        // Compute bright limb angle bl
        let bl = (.pi + atan2(cos(sunDec) * sin(moonRA - sunRA), cos(sunDec) * sinMoonDec * cos(moonRA - sunRA) - sin(sunDec) * cosMoonDec))

        // Paralactic angle par
        y = sin(lst - moonRA)
        x = tan(obsLat) * cosMoonDec - sinMoonDec * cos(lst - moonRA)
        var par = 0.0
        if x != 0.0 {
            par = atan2(y, x)
        } else {
            par = (y / abs(y)) * PI_OVER_TWO
        }

        return .init(
            axisPosition: .init(value: p, unit: .radians), brightLimb: .init(value: bl, unit: .radians), paralactic: .init(value: par, unit: .radians),
            opticalLibration: (.init(value: lp, unit: .radians), .init(value: bp, unit: .radians))
        )
    }
}

public extension Moon {
    var diskOrientationAngles: DiskOrientationAngles { calculation.diskOrientationAngles }
}

public struct ViewingAngles {
    /// Phase angle which is related to the illumination
    public var phase: Measurement<UnitAngle>

    /// Angle of moon rotation
    public var moon: Measurement<UnitAngle>

    /// Angle of shadow rotation
    public var shadow: Measurement<UnitAngle>
}

public extension Moon {
    var diskOrientationViewingAngles: ViewingAngles {
        .init(
            phase: Measurement(value: acos(-cos(age / LUNAR_CYCLE_DAYS * 2 * .pi)), unit: .radians).converted(to: .degrees),
            moon: (-1 * (diskOrientationAngles.brightLimb - diskOrientationAngles.paralactic)).converted(to: .degrees),
            shadow: (-1 * (diskOrientationAngles.axisPosition - diskOrientationAngles.paralactic)).converted(to: .degrees)
        )
    }
}

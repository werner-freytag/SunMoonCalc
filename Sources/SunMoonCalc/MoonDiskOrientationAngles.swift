import Foundation

public struct DiskOrientationAngles: Equatable {
    public static func == (lhs: DiskOrientationAngles, rhs: DiskOrientationAngles) -> Bool {
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
        let moonLon = moonData.azimuth * DEG_TO_RAD
        let moonLat = moonData.elevation * DEG_TO_RAD
        let moonRA = moonData.rightAscension
        let moonDec = moonData.declination

        let t = jd.t

        // Moon's argument of latitude
        let F = (93.2720993 + 483_202.0175273 * t - 0.0034029 * t * t - t * t * t / 3_526_000 + t * t * t * t / 863_310_000) * DEG_TO_RAD

        // Moon's inclination
        let I = 1.54242 * DEG_TO_RAD

        // Moon's mean ascending node longitude
        let omega = (125.0445550 - 1934.1361849 * t + 0.0020762 * t * t + t * t * t / 467_410 - t * t * t * t / 18_999_000) * DEG_TO_RAD

        // Obliquity of ecliptic (approx, better formulae up)
        let eps = 23.43929 * DEG_TO_RAD

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
    var phase: Measurement<UnitAngle>

    /// Angle of moon rotation
    var moon: Measurement<UnitAngle>

    /// Angle of shadow rotation
    var shadow: Measurement<UnitAngle>
}

extension Moon {
    var diskOrientationViewingAngles: ViewingAngles {
        .init(
            phase: .init(value: acos(-cos(age / LUNAR_CYCLE_DAYS * 2 * .pi)), unit: .radians).converted(to: .degrees),
            moon: (-1 * (diskOrientationAngles.brightLimb - diskOrientationAngles.paralactic)).converted(to: .degrees),
            shadow: (-1 * (diskOrientationAngles.axisPosition - diskOrientationAngles.paralactic)).converted(to: .degrees)
        )
    }
}

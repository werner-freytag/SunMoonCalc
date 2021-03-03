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

/// Method to calculate values of Moon Disk
/// - Returns: [optical librations (lp), lunar coordinates of the centre of the disk (bp), position angle of axis (p), bright limb angle (bl), paralactic angle (par)]
extension MoonCalculation {
    var diskOrientationAngles: DiskOrientationAngles {
        let sunData = SunCalculation(jd: jd, obsLat: obsLat, obsLon: obsLon, twilight: twilight).ephemerisData
        let moonData = ephemerisData

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

        return .init(
            axisPosition: .init(value: p, unit: .radians), brightLimb: .init(value: bl, unit: .radians), paralactic: .init(value: par, unit: .radians),
            opticalLibration: (.init(value: lp, unit: .radians), .init(value: lp, unit: .radians))
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

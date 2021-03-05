@testable import SunMoonCalc
import XCTest

final class SunMoonCalcTests: XCTestCase {
    private lazy var dateFormatter: DateFormatter = {
        let dateFormatter = DateFormatter()
        dateFormatter.dateFormat = "yyyy-MM-dd HH:mm:ss Z"
        return dateFormatter
    }()

    func testSun() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let sun = Sun(location: location, date: date)

        XCTAssertEqual(sun.ephemeris, Ephemeris(
            azimuth: Measurement(value: 0.22479019677783363, unit: .radians),
            elevation: Measurement(value: -0.8718445884699833, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 06:00:25 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 16:53:26 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 11:26:32 +0000")!,
            transitElevation: Measurement(value: 0.5820695076646274, unit: .radians),
            distance: Measurement(value: 0.9900659837432739, unit: .astronomicalUnits)
        ))
    }

    func testMoon() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let moon = Moon(location: location, date: date)

        XCTAssertEqual(moon.ephemeris, Ephemeris(
            azimuth: Measurement(value: 3.8685629085102007, unit: .radians),
            elevation: Measurement(value: 0.953334406991476, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 15:44:09 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 06:02:35 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 23:15:00 +0000")!,
            transitElevation: Measurement(value: 0.9784585056858743, unit: .radians),
            distance: Measurement(value: 0.002479871337300747, unit: .astronomicalUnits)
        ))

        XCTAssertEqual(moon.age, 13.325133319225687)
        XCTAssertEqual(moon.illumination, 0.9753461326608158)
    }

    func testMoonDiskOrientation() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let result = Moon(location: location, date: date).diskOrientationAngles

        XCTAssertEqual(result, .init(
            axisPosition: Measurement(value: 0.3391984841184668, unit: .radians),
            brightLimb: Measurement(value: 5.252281717032904, unit: .radians),
            paralactic: Measurement(value: 0.48687916720474905, unit: .radians),
            opticalLibration: (longPass: Measurement(value: 6.182599205266342, unit: .radians),
                               bandPass: Measurement(value: -0.10274916497635754, unit: .radians))
        ))
    }
}

extension Ephemeris: Equatable {
    public static func == (lhs: Ephemeris, rhs: Ephemeris) -> Bool {
        return lhs.azimuth == rhs.azimuth
            && lhs.elevation == rhs.elevation
            && floor(lhs.rise?.timeIntervalSince1970 ?? 0) == floor(rhs.rise?.timeIntervalSince1970 ?? 0)
            && floor(lhs.set?.timeIntervalSince1970 ?? 0) == floor(rhs.set?.timeIntervalSince1970 ?? 0)
            && floor(lhs.transit?.timeIntervalSince1970 ?? 0) == floor(rhs.transit?.timeIntervalSince1970 ?? 0)
            && lhs.transitElevation == rhs.transitElevation
            && lhs.distance == rhs.distance
    }
}

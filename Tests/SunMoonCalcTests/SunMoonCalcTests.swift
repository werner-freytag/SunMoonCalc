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
            azimuth: Measurement(value: 0.2250523624268812, unit: .radians),
            elevation: Measurement(value: -0.8718563214579441, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 06:00:24 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 16:53:23 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 11:26:30 +0000")!,
            transitElevation: Measurement(value: 0.5822069829475769, unit: .radians),
            distance: Measurement(value: 0.9901094623869873, unit: .astronomicalUnits)
        ))
    }

    func testMoon() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let moon = Moon(location: location, date: date)

        XCTAssertEqual(moon.ephemeris, Ephemeris(
            azimuth: Measurement(value: 3.869137520615342, unit: .radians),
            elevation: Measurement(value: 0.9532487104797652, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 15:44:03 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 06:02:31 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 23:14:56 +0000")!,
            transitElevation: Measurement(value: 0.980261403394202, unit: .radians),
            distance: Measurement(value: 0.0024798996480156228, unit: .astronomicalUnits)
        ))

        XCTAssertEqual(moon.age, 13.324321196424968)
        XCTAssertEqual(moon.illumination, 0.9753196030933413)
    }

    func testMoonDiskOrientation() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let result = Moon(location: location, date: date).diskOrientationAngles

        XCTAssertEqual(result, .init(
            axisPosition: Measurement(value: 0.3373163619617799, unit: .radians),
            brightLimb: Measurement(value: 5.252056936039397, unit: .radians),
            paralactic: Measurement(value: 0.48723530437531215, unit: .radians),
            opticalLibration: (longPass: Measurement(value: 3.806970357998125, unit: .radians),
                               bandPass: Measurement(value: 3.806970357998125, unit: .radians))
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

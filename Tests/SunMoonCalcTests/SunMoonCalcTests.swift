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
            azimuth: Measurement(value: 0.2248122999945994, unit: .radians),
            elevation: Measurement(value: -0.871847920500529, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 06:00:25 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 16:53:26 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 11:26:32 +0000")!,
            transitElevation: Measurement(value: 0.5820640491655357, unit: .radians),
            distance: Measurement(value: 0.9900657900481319, unit: .astronomicalUnits)
        ))
    }

    func testMoon() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let moon = Moon(location: location, date: date)

        XCTAssertEqual(moon.ephemeris, Ephemeris(
            azimuth: Measurement(value: 3.8689066714962914, unit: .radians),
            elevation: Measurement(value: 0.9532944455050057, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 15:44:05 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 06:02:33 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 23:14:58 +0000")!,
            transitElevation: Measurement(value: 0.9785343943804162, unit: .radians),
            distance: Measurement(value: 0.00247989875343263, unit: .astronomicalUnits)
        ))

        XCTAssertEqual(moon.age, 13.18057635505499)
        XCTAssertEqual(moon.illumination, 0.9753178310696914)
    }

    func testMoonDiskOrientation() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let result = Moon(location: location, date: date).diskOrientationAngles

        XCTAssertEqual(result, .init(
            axisPosition: Measurement(value: 0.3373337857141483, unit: .radians),
            brightLimb: Measurement(value: 5.252072231687541, unit: .radians),
            paralactic: Measurement(value: 0.48709403178431765, unit: .radians),
            opticalLibration: (longPass: Measurement(value: 3.806966322825474, unit: .radians),
                               bandPass: Measurement(value: 0.008969020286691026, unit: .radians))
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

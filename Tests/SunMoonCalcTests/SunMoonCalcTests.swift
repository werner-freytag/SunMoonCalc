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
            azimuth: Measurement(value: 0.22491823240857967, unit: .radians),
            elevation: Measurement(value: -0.8718361145706249, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 06:00:24 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 16:53:25 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 11:26:31 +0000")!,
            transitElevation: Measurement(value: 0.5817304606714634, unit: .radians),
            distance: Measurement(value: 0.9900657897242581, unit: .astronomicalUnits)
        ))
    }

    func testMoon() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let moon = Moon(location: location, date: date)

        XCTAssertEqual(moon.ephemeris, Ephemeris(
            azimuth: Measurement(value: 3.869007768509518, unit: .radians),
            elevation: Measurement(value: 0.9532608338997101, unit: .radians),
            rise: dateFormatter.date(from: "2021-02-26 15:44:04 +0000")!,
            set: dateFormatter.date(from: "2021-02-26 06:02:32 +0000")!,
            transit: dateFormatter.date(from: "2021-02-26 23:14:57 +0000")!,
            transitElevation: Measurement(value: 0.9788817812131069, unit: .radians),
            distance: Measurement(value: 0.0024798995772362955, unit: .astronomicalUnits)
        ))

        XCTAssertEqual(moon.age, 13.18057635505499)
        XCTAssertEqual(moon.illumination, 0.9753177382119795)
    }

    func testMoonDiskOrientation() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let location = Location(48.137154, 11.576124)

        let result = Moon(location: location, date: date).diskOrientationAngles

        XCTAssertEqual(result, .init(
            axisPosition: Measurement(value: 0.337333364797035, unit: .radians),
            brightLimb: Measurement(value: 5.252067047999782, unit: .radians),
            paralactic: Measurement(value: 0.48715374686552315, unit: .radians),
            opticalLibration: (longPass: Measurement(value: 3.8069680919395523, unit: .radians),
                               bandPass: Measurement(value: 3.8069680919395523, unit: .radians))
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

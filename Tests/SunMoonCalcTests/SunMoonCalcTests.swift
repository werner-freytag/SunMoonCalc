@testable import SunMoonCalc
import XCTest

final class SunMoonCalcTests: XCTestCase {
    private lazy var dateFormatter: DateFormatter = {
        let dateFormatter = DateFormatter()
        dateFormatter.dateFormat = "yyyy-MM-dd HH:mm:ss Z"
        return dateFormatter
    }()
    
    func testCalculation() {
        let date = dateFormatter.date(from: "2021-02-26 00:00:00 +0000")!
        let result = try! calcSunAndMoon(date: date, latitude: 48.137154, longitude: 11.576124)

        XCTAssertEqual(result.sun, Sun(
            ephemeris: Ephemeris(
                azimuth: Measurement(value: 0.22503024947665207, unit: .radians),
                elevation: Measurement(value: -0.8718529905677839, unit: .radians),
                rise: dateFormatter.date(from: "2021-02-26 06:00:24 +0000")!,
                set: dateFormatter.date(from: "2021-02-26 16:53:24 +0000")!,
                transit: dateFormatter.date(from: "2021-02-26 11:26:30 +0000")!,
                transitElevation: Measurement(value: 0.5822120999129049, unit: .radians),
                distance: Measurement(value: 0.9901096566337206, unit: .astronomicalUnits)
            )
        ))

        XCTAssertEqual(result.moon, Moon(
            ephemeris: Ephemeris(
                azimuth: Measurement(value: 3.8687937729816104, unit: .radians),
                elevation: Measurement(value: 0.9532887065772591, unit: .radians),
                rise: dateFormatter.date(from: "2021-02-26 15:44:07 +0000")!,
                set: dateFormatter.date(from: "2021-02-26 06:02:33 +0000")!,
                transit: dateFormatter.date(from: "2021-02-26 23:14:58 +0000")!,
                transitElevation: Measurement(value: 0.9802353344085494, unit: .radians),
                distance: Measurement(value: 0.002479872229507751, unit: .astronomicalUnits)
            ),
            age: 13.325196027873513,
            illumination: 0.975347904390796,
            positionAngleOfAxis: Measurement(value: 0.3373622533190958, unit: .radians),
            brightLimbAngle: Measurement(value: 5.25226644381023, unit: .radians),
            paralacticAngle: Measurement(value: 0.48702047010300514, unit: .radians)
        ))
    }
}

@testable import SunMoonCalc
import XCTest

final class AngleHelpersTests: XCTestCase {
    func testNormalizeAngular() {
        XCTAssertEqual(3 - 1, normalizeAngle(value: -3.0 - 1, max: 3))
        XCTAssertEqual(3 - 1, normalizeAngle(value: -1.0, max: 3))
        XCTAssertEqual(3 - 0.1, normalizeAngle(value: -0.1, max: 3))
        XCTAssertEqual(0, normalizeAngle(value: 0, max: 3))
        XCTAssertEqual(3 - 0.1, normalizeAngle(value: 3 - 0.1, max: 3))
        XCTAssertEqual(0, normalizeAngle(value: 3, max: 3))
        XCTAssertEqual(3 - 1, normalizeAngle(value: 2 * 3.0 - 1, max: 3))
    }

    func testNormalizeDegrees() {
        XCTAssertEqual(350, normalizeAngle(value: -10.0, max: 360))
    }
}

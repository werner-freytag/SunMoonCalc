@testable import SunMoonCalc
import XCTest

final class SkyDirectionTests: XCTestCase {
    func testDirections() {
        XCTAssertEqual(.N, SkyDirection(angle: Measurement<UnitAngle>(value: 0, unit: .degrees)))
        XCTAssertEqual(.E, SkyDirection(angle: Measurement<UnitAngle>(value: 90, unit: .degrees)))
        XCTAssertEqual(.S, SkyDirection(angle: Measurement<UnitAngle>(value: 180, unit: .degrees)))
        XCTAssertEqual(.W, SkyDirection(angle: Measurement<UnitAngle>(value: 270, unit: .degrees)))

        XCTAssertEqual(.N, SkyDirection(angle: Measurement<UnitAngle>(value: 0, unit: .radians)))
        XCTAssertEqual(.E, SkyDirection(angle: Measurement<UnitAngle>(value: .pi / 2, unit: .radians)))
        XCTAssertEqual(.S, SkyDirection(angle: Measurement<UnitAngle>(value: .pi, unit: .radians)))
        XCTAssertEqual(.W, SkyDirection(angle: Measurement<UnitAngle>(value: .pi * 1.5, unit: .radians)))

        XCTAssertEqual(.NE, SkyDirection(angle: Measurement<UnitAngle>(value: 45, unit: .degrees)))
        XCTAssertEqual(.NNE, SkyDirection(angle: Measurement<UnitAngle>(value: 22.5, unit: .degrees)))

        // border tests
        XCTAssertEqual(.NNE, SkyDirection(angle: Measurement<UnitAngle>(value: 11.25, unit: .degrees)))
        XCTAssertEqual(.N, SkyDirection(angle: Measurement<UnitAngle>(value: 11.249, unit: .degrees)))
        XCTAssertEqual(.N, SkyDirection(angle: Measurement<UnitAngle>(value: -11.25, unit: .degrees)))
        XCTAssertEqual(.NNW, SkyDirection(angle: Measurement<UnitAngle>(value: -11.251, unit: .degrees)))
    }
}

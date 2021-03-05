import Foundation

public enum SkyDirection: CaseIterable {
    case N, NNE, NE, ENE, E, ESE, SE, SSE, S, SSW, SW, WSW, W, WNW, NW, NNW
    public init(angle: Measurement<UnitAngle>) {
        let index = Int(floor(normalizeDegrees(angle.converted(to: .degrees).value + 360 / 32) / 360 * 16)) % 16
        self = SkyDirection.allCases[index]
    }
}

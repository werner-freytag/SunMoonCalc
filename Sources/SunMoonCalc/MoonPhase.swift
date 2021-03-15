/// Enum describing the phase of the moon
public extension Moon {
    enum Phase {
        case newMoon
        case waxingCrescent
        case firstQuarter
        case waxingGibbous
        case fullMoon
        case waningGibbous
        case lastQuarter
        case waningCrescent

        /// Init per moon age (days)
        /// - Parameter age: Lunar age in days
        public init?(age: Double) {
            guard (0 ... LUNAR_CYCLE_DAYS).contains(age) else { return nil }

            switch age {
            case 1 ..< 6.4:
                self = .waxingCrescent
            case 6.4 ..< 8.4:
                self = .firstQuarter
            case 8.4 ..< 13.8:
                self = .waxingGibbous
            case 13.8 ..< 15.8:
                self = .fullMoon
            case 15.8 ..< 21.1:
                self = .waningGibbous
            case 21.1 ..< 23.1:
                self = .lastQuarter
            case 23.1 ..< LUNAR_CYCLE_DAYS - 1:
                self = .waningCrescent
            default:
                self = .newMoon
            }
        }
    }
}

public extension Moon {
    var phase: Phase { .init(age: phaseAge)! }
}

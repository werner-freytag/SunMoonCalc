/** Method to return the phase of the Moon as per Moon Age (days)
 - parameter lunarAge:    Lunar age in days
 - returns:    Phase of the Moon */
public func getMoonPhaseName(lunarAge: Double) -> String {
    if lunarAge >= 0, lunarAge <= LUNAR_CYCLE_DAYS, lunarAge < 1 || lunarAge > LUNAR_CYCLE_DAYS - 1 {
        return "New Moon"
    } else if lunarAge >= 1, lunarAge < 6.4 {
        return "Waxing Crescent"
    } else if lunarAge >= 6.4, lunarAge < 8.4 {
        return "First Quarter"
    } else if lunarAge >= 8.4, lunarAge < 13.8 {
        return "Waxing Gibbous"
    } else if lunarAge >= 13.8, lunarAge < 15.8 {
        return "Full Moon"
    } else if lunarAge >= 15.8, lunarAge < 21.1 {
        return "Waning Gibbous"
    } else if lunarAge >= 21.1, lunarAge < 23.1 {
        return "Last/Third Quarter"
    } else if lunarAge >= 23.1, lunarAge <= LUNAR_CYCLE_DAYS - 1 {
        return "Waning Crescent"
    } else {
        return "-"
    }
}

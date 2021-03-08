// swift-tools-version:5.3
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "SunMoonCalc",
    platforms: [
        .macOS(.v10_12),
    ],
    products: [
        // Products define the executables and libraries a package produces, and make them visible to other packages.
        .library(
            name: "SunMoonCalc",
            targets: ["SunMoonCalc"]
        ),
    ],
    dependencies: [
        .package(url: "https://github.com/werner-freytag/SwiftToolbox.git", .branch("develop")),
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages this package depends on.
        .target(
            name: "SunMoonCalc",
            dependencies: ["SwiftToolbox"]
        ),
        .testTarget(
            name: "SunMoonCalcTests",
            dependencies: ["SunMoonCalc"]
        ),
    ]
)

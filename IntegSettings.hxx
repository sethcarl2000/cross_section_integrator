#ifndef IntegSettings_HH
#define IntegSettings_HH

namespace Setting {
    enum Bit : int {
        kNone = 0,
        kConserveEnergy     = 1 << 0, //enforce external energy conservation for each update
        kConserveMomenta    = 1 << 1, //enforce external momentum conservation for each update 
        kAutoAdjustScan     = 1 << 2, //automatically adjust the scan rate to keep the acceptance prob around ~50% 
        kVerbose            = 1 << 3  //if true, then report output
    };
}

#endif
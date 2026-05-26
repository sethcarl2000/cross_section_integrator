#ifndef IntegSettings_HH
#define IntegSettings_HH

namespace Setting {
    enum Bit : int {
        kNone = 0,
        kConserveEnergy     = 1 << 0, //e nforce external energy conservation for each update
        kConserveMomenta    = 1 << 1, //enforce external momentum conservation for each update 
        kAutoAdjustScan     = 1 << 2, //automatically adjust the scan rate to keep the acceptance prob around ~50% 
        kVerbose            = 1 << 3, //if true, then report output
    }; 

    //specifcy which algorithm to use from the GSL monte-carlo integration library
    enum MCStrategy : int {
        kPLAIN      = 1 << 4,  
        kMISER      = 1 << 5, 
        kVEGAS      = 1 << 6, 
        kADAPTIVE   = 1 << 7, 
        kMETROPOLIS = 1 << 8
    };
    //ZZZZZZZZZZZZZZZ
    //lkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk,566666666666666666666 bnnnnnnnnnnnnnnnnn
    //   -- muon's suggestion (20 May 26)
}

#endif
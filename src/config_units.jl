# define Sverdrups
module UnitfulOcean; using Unitful; 
@unit sverdrup "Sv" Sverdrup (10^6)u"m^3/s" false;
end

Unitful.register(UnitfulOcean)

const kg = u"kg"
const m = u"m"
const yr = u"yr"
const Sv = u"sverdrup"
const km = u"km"

ENV["UNITFUL_FANCY_EXPONENTS"] = true

const lam=268e-9;#NOTE: choosing an 'average' between d1 and d2.  sorry if it winds up mattering; normally these are close enough taht it's fine.
const gam = 2 * pi * 26e6;#NOTE: choosing an 'average' between d1 and d2.  sorry if it winds up mattering; normally these are close enough taht it's fine.
const normalizedBohrMag = 0.054; #\mu_{B}/\hbar\Gamma in units 1/Gauss
const mass = (197) * 1.67e-27;#mass of Au
const excitedEnergyD1 = [-15.3,9.2];#Nuclear Physics A580 (1994) 173-212; ISOLDE Collaboration: D1 relevant hyperfine energies (note: will only have 2 entries for molecules where this matters (Ag, K, Li, etc.)).  Always ordered such that 'higher' F corresponds to last value
const excitedEnergyD2 = [-5.0,3.0];#Not accurate: Real info in PRA 50, 209 (1994) Svanberg paper (also see Ref [4] of that paper)
# also see PRA 47, 4725 (1993) Matthias paper for how to use A & B hyperfine params to calculate energy levels
const nucSpin = 3/2;

const kA = 2 * pi / lam; #wavenumber
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2

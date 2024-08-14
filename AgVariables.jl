const lam=333e-9;#NOTE: choosing an 'average' between d1 and d2.  sorry if it winds up mattering; normally these are close enough taht it's fine.
const gam = 2 * pi * 22.5e6;#NOTE: choosing an 'average' between d1 and d2.  sorry if it winds up mattering; normally these are close enough taht it's fine.
const normalizedBohrMag = 0.061; #\mu_{B}/\hbar\Gamma in units 1/Gauss
const mass = (107) * 1.67e-27;#mass of Ag
const excitedEnergyD1 = [7.2,-2.4];#D1 relevant hyperfine energies (note: will only have 2 entries for molecules where this matters (Ag, K, Li, etc.)).  Always ordered such that 'higher' F corresponds to last value
const excitedEnergyD2 = [2.4,-.8];#D2 relevant hyperfine energies (note: will only have 2 entries for molecules where this matters (Ag, K, Li, etc.)) Always ordered such that 'higher' F corresponds to last value
const nucSpin = 1/2;

const kA = 2 * pi / lam; #wavenumber
const velFactor = (gam/kA);
const hbar=1.05e-34;
const accelFactor = (1e-3*hbar*kA*gam/mass);#normalized force units in program are 1e-3\hbar*k*\gam.  So, the factor converts this to m/s^2
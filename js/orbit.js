//Based on pykep.propagate_lagrangian
function propagate(r0, v0, tof, mu) {
	let r = r0.norm;
	let v = v0.norm;

	let energy = (v * v / 2 - mu / r);
	let a = -mu / 2 / energy;

	let dM = Math.sqrt(Math.abs(mu / (a * a * a))) * tof;

	let C1 = dot3(r0, v0) / Math.sqrt(Math.abs(mu * a));
	let C2 = 1 - r / a;

	let F, G, Ft, Gt;

	if (energy < 0) {
		// Elliptical case
		// Kepler equation: M       = E - e * sin(E)
		// diff:            M1 - M0 = E1 - e * sin(E1) - (E0 - e * sin(E0))
		// diff:            M1 - M0 = E1 - e * sin(E1) - E0 + e * sin(E0)
		// diff:            dM      = dE - e * sin(E1) + e * sin(E0)
		// diff:            dM      = dE + e * (sin(E0) - sin(E1))
		// diff:            dM      = dE + e * (sin(E0) - sin(E0 + dE))
		// diff:            dM      = dE + e * (sin(E0) - sin(E0)cos(dE) - cos(E0)sin(dE))
		// diff:            dM      = dE + e * (sin(E0)(1 - cos(dE)) - cos(E0)sin(dE))
		// diff:            dM      = dE + e * sin(E0)(1 - cos(dE)) - e * cos(E0)sin(dE)
		// Since r = a(1 - e * cos(E0)), e * cos(E0) = 1 - r / a
		// diff:            dM      = dE + e * sin(E0)(1 - cos(dE)) - (1 - r / a)sin(dE)
		// e * sin(E0) = dot(r0, v0) / sqrt(mu * a), derivation might be posted later
		//
		// +============= Difference form of Kepler equation for ellipses =============+
		// | dM = dE + dot(r0, v0) * (1 - cos(dE)) / sqrt(mu * a) - (1 - r / a)sin(dE) |
		// +===========================================================================+
		//
		// Setting up Newton's method:
		// f(dE) = dE + dot(r0, v0) * (1 - cos(dE)) / sqrt(mu * a) - (1 - r / a)sin(dE) - dM
		// f'(dE) = 1 + dot(r0, v0) * sin(dE) / sqrt(mu * a) - (1 - r / a)cos(dE)

		let dE = dM;

		// Newton's method
		for (let i = 0; i < 5; i++) {
			let f = dE + C1 * (1 - Math.cos(dE)) - C2 * Math.sin(dE) - dM;
			let df = 1 + C1 * Math.sin(dE) - C2 * Math.cos(dE);

			dE -= f / df;
		}

		let r2 = a + (r - a) * Math.cos(dE) + dot3(r0, v0) * Math.sqrt(a / mu) * Math.sin(dE);

		F = 1 - a / r * (1 - Math.cos(dE));
		G = a * dot3(r0, v0) / mu * (1 - Math.cos(dE)) + r * Math.sqrt(a / mu) * Math.sin(dE);
		Ft = -Math.sqrt(mu * a) / (r2 * r) * Math.sin(dE);
		Gt = 1 - a / r2 * (1 - Math.cos(dE));
	} else {
		// Hyperbolic case
		// Kepler equation: M       = e * sinh(H) - H
		// diff:            M1 - M0 = e * sinh(H1) - H1 - e * sinh(H0) + H0
		// diff:            dM      = e * (sinh(H1) - sinh(H0)) - dH
		// diff:            dM      = e * (sinh(H0 + dH) - sinh(H0)) - dH
		// diff:            dM      = e * (sinh(H0)cosh(dH) + cosh(H0)sinh(dH) - sinh(H0)) - dH
		// diff:            dM      = e * (sinh(H0)(cosh(dH) - 1) + cosh(H0)sinh(dH)) - dH
		// diff:            dM      = e * sinh(H0)(cosh(dH) - 1) + e * cosh(H0)sinh(dH) - dH
		// Since r = -a(e * cosh(H0) - 1), e * cosh(H0) = 1 - r / a
		// diff:            dM      = e * sinh(H0)(cosh(dH) - 1) + (1 - r / a)sinh(dH) - dH
		// Similarly to the elliptical orbit, 
		// e * sinh(H0) = dot(r0, v0) / sqrt(-mu * a).
		//
		// +============= Difference form of Kepler equation for hyperbolas ==============+
		// | dM = dot(r0, v0) * (cosh(dH) - 1) / sqrt(-mu * a) + (1 - r / a)sinh(dH) - dH |
		// +==============================================================================+
		//
		// Setting up Newton's method:
		// f(dH)  = dot(r0, v0) * (cosh(dH) - 1) / sqrt(-mu * a) + (1 - r / a)sinh(dH) - dH - dM
		// f'(dH) = dot(r0, v0) * sinh(dH) / sqrt(-mu * a)       + (1 - r / a)cosh(dH) - 1

		let dH = (tof > 0) ? 1 : -1;

		//Newton's method
		for (let i = 0; i < 5; i++) {
			let f = C1 * (Math.cosh(dH) - 1) + C2 * Math.sinh(dH) - dH - dM;
			let df = C1 * Math.sinh(dH) + C2 * Math.cosh(dH) - 1;

			dH -= f / df;
		}

		let r2 = a + (r - a) * Math.cosh(dH) + dot3(r0, v0) * Math.sqrt(-a / mu) * Math.sinh(dH);

		F = 1 - a / r * (1 - Math.cosh(dH));
		G = a * dot3(r0, v0) / mu * (1 - Math.cosh(dH)) + r * Math.sqrt(-a / mu) * Math.sinh(dH);
		Ft = -Math.sqrt(-mu * a) / (r2 * r) * Math.sinh(dH);
		Gt = 1 - a / r2 * (1 - Math.cosh(dH));

		console.log(F, G, Ft, Gt);
	}

	return {
		r: add3(mult3(r0, F),  mult3(v0, G)),
		v: add3(mult3(r0, Ft), mult3(v0, Gt))
	};
}
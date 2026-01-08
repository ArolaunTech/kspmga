import { Vector3 } from './vector.js';

export class Orbit {
	constructor(mu, sma, eccentricity, inclination, lan, argp, epoch, meananomalyatepoch) {
		this.mu = mu;
		this.sma = sma;
		this.eccentricity = eccentricity;
		this.inclination = inclination;
		this.longascendingnode = lan;
		this.argperiapsis = argp;
		this.epoch = epoch;
		this.meananomalyatepoch = meananomalyatepoch;
	}

	get meanmotion() {
		return Math.sqrt(this.mu / Math.abs(this.sma * this.sma * this.sma));
	}

	get period() {
		if (this.sma < 0) return Infinity;
		return 2 * Math.PI / this.meanmotion;
	}

	get semilatusrectum() {
		return this.sma * (1 - this.eccentricity * this.eccentricity);
	}

	get specificangularmomentum() {
		return Math.sqrt(this.semilatusrectum * this.mu);
	}

	get semiminoraxis() {
		return this.sma * Math.sqrt(Math.abs(this.eccentricity * this.eccentricity - 1));
	}

	get periapsis() {
		return this.sma * (1 - this.eccentricity);
	}

	get apoapsis() {
		return this.sma * (1 + this.eccentricity);
	}

	fromRV(t, r0, v0, mu) {
		// RV to COE
		this.mu = mu;
		this.epoch = t;

		let dist0 = r0.norm;
		let vel0 = v0.norm;
		let rv0 = Vector3.dot(r0, v0) / dist0;

		let energy = vel0 * vel0 / 2 - mu / dist0;

		this.sma = -mu / 2 / energy;

		let hvector = Vector3.cross(r0, v0);
		let h = hvector.norm;
		
		let evector = Vector3.sub(Vector3.mult(Vector3.cross(v0, hvector), 1 / mu), Vector3.normalize(r0));

		this.eccentricity = evector.norm;
		this.inclination = Math.acos(hvector.z / h);

		let an = Vector3.cross(new Vector3(0, 0, 1), hvector);

		this.longascendingnode = Math.atan2(an.y, an.x);
		this.argperiapsis = Math.acos(Vector3.dot(evector, an) / this.eccentricity / an.norm);

		if (evector.z < 0) this.argperiapsis = 2 * Math.PI - this.argperiapsis;

		let cosE = (1 - dist0 / this.sma) / this.eccentricity;
		let E;

		if (this.sma > 0) {
			E = Math.acos(cosE);
		} else {
			E = Math.acosh(cosE);
		}

		if (rv0 < 0) {
			E *= -1;
		}

		if (this.sma > 0) {
			this.meananomalyatepoch = E - this.eccentricity * Math.sin(E);
		} else {
			this.meananomalyatepoch = this.eccentricity * Math.sinh(E) - E;
		}
	}

	meananomaly(time) {
		let anomaly = this.meananomalyatepoch + (time - this.epoch) * this.meanmotion;

		if (this.sma > 0) {
			anomaly = anomaly % (2 * Math.PI);

			if (anomaly < 0) anomaly += 2 * Math.PI;
		}

		return anomaly;
	}

	eccentricanomaly(time) {
		let meananomaly = this.meananomaly(time);

		if (this.sma > 0) {
			// Elliptical case
			// Kepler equation: M = E - e * sin(E)
			// f(E)  = E - e * sin(E) - M
			// f'(E) = 1 - e * cos(E)

			let E = meananomaly;

			for (let i = 0; i < 5; i++) {
				let f = E - this.eccentricity * Math.sin(E) - meananomaly;
				let df = 1 - this.eccentricity * Math.cos(E);

				E -= f / df;
			}

			return E;
		} else {
			// Hyperbolic case
			// Kepler equation: M = e * sinh(H) - H
			// f(H)  = e * sinh(H) - H - M
			// f'(H) = e * cosh(H) - 1

			let H = Math.cbrt(meananomaly); //Extremely good guess

			for (let i = 0; i < 5; i++) {
				let f = this.eccentricity * Math.sinh(H) - H - meananomaly;
				let df = this.eccentricity * Math.cosh(H) - 1;

				H -= f / df;
			}

			return H;
		}
	}

	trueanomaly(time) {
		let E = this.eccentricanomaly(time);
		let sqrtquotient = Math.sqrt(Math.abs(
			(this.eccentricity + 1) / (this.eccentricity - 1)
		));

		let tanE;

		if (this.sma > 0) {
			//Elliptical case
			tanE = Math.tan(E / 2);
		} else {
			//Hyperbolic case
			tanE = Math.tanh(E / 2);
		}

		return 2 * Math.atan(tanE * sqrtquotient);
	}

	distance(time) {
		let E = this.eccentricanomaly(time);

		if (this.sma > 0) {
			return this.sma * (1 - this.eccentricity * Math.cos(E));
		}
		return this.sma * (1 - this.eccentricity * Math.cosh(E));
	}

	speed(time) {
		let dist = this.distance(time);
		return Math.sqrt(this.mu * (2 / dist - 1 / this.sma));
	}

	hvel(time) {
		let dist = this.distance(time);
		return this.specificangularmomentum / dist;
	}

	vvel(time) {
		let absvvel = Math.sqrt(
			Math.pow(this.speed(time), 2) - Math.pow(this.hvel(time), 2)
		);

		let meananomaly = this.meananomaly(time);
		if (meananomaly < 0) return -absvvel;
		if (this.sma < 0) return absvvel;
		if (meananomaly < Math.PI) return absvvel;
		return -absvvel;
	}

	getState(time) {
		let E = this.eccentricanomaly(time);
		let f = this.trueanomaly(time);

		let latusspeed = Math.sqrt(this.mu / this.semilatusrectum);
		let P = latusspeed * this.eccentricity;

		let pos;
		let vel = new Vector3(
			-latusspeed * Math.sin(f),
			latusspeed * Math.cos(f) + P,
			0
		);

		if (this.sma > 0) {
			// Elliptical case - perifocal coordinate system
			pos = new Vector3(
				this.sma * (Math.cos(E) - this.eccentricity), 
				this.semiminoraxis * Math.sin(E),
				0
			);
		} else {
			// Hyperbolic case - perifocal coordinate system
			pos = new Vector3(
				this.sma * (Math.cosh(E) - this.eccentricity), 
				this.semiminoraxis * Math.sinh(E), 
				0
			);
		}

		// Rotate by argument of periapsis
		let cosargp = Math.cos(this.argperiapsis);
		let sinargp = Math.sin(this.argperiapsis);

		pos = new Vector3(
			pos.x * cosargp - pos.y * sinargp,
			pos.y * cosargp + pos.x * sinargp,
			0
		);

		vel = new Vector3(
			vel.x * cosargp - vel.y * sinargp,
			vel.y * cosargp + vel.x * sinargp,
			0
		);
		// Rotate by inclination
		let cosinc = Math.cos(this.inclination);
		let sininc = Math.sin(this.inclination);

		pos = new Vector3(
			pos.x,
			pos.y * cosinc,
			pos.y * sininc
		);

		vel = new Vector3(
			vel.x,
			vel.y * cosinc,
			vel.y * sininc
		);
		// Rotate by longitude of ascending node
		let coslan = Math.cos(this.longascendingnode);
		let sinlan = Math.sin(this.longascendingnode);

		pos = new Vector3(
			pos.x * coslan - pos.y * sinlan,
			pos.y * coslan + pos.x * sinlan,
			pos.z
		);

		vel = new Vector3(
			vel.x * coslan - vel.y * sinlan,
			vel.y * coslan + vel.x * sinlan,
			vel.z
		);

		return {r: pos, v: vel};
	}

	getPosFromE(E) {
		let pos;

		if (this.sma > 0) {
			// Elliptical case - perifocal coordinate system
			pos = new Vector3(
				this.sma * (Math.cos(E) - this.eccentricity), 
				this.semiminoraxis * Math.sin(E),
				0
			);
		} else {
			// Hyperbolic case - perifocal coordinate system
			pos = new Vector3(
				this.sma * (Math.cosh(E) - this.eccentricity), 
				this.semiminoraxis * Math.sinh(E), 
				0
			);
		}

		// Rotate by argument of periapsis
		let cosargp = Math.cos(this.argperiapsis);
		let sinargp = Math.sin(this.argperiapsis);

		pos = new Vector3(
			pos.x * cosargp - pos.y * sinargp,
			pos.y * cosargp + pos.x * sinargp,
			0
		);
		// Rotate by inclination
		let cosinc = Math.cos(this.inclination);
		let sininc = Math.sin(this.inclination);

		pos = new Vector3(
			pos.x,
			pos.y * cosinc,
			pos.y * sininc
		);
		// Rotate by longitude of ascending node
		let coslan = Math.cos(this.longascendingnode);
		let sinlan = Math.sin(this.longascendingnode);

		pos = new Vector3(
			pos.x * coslan - pos.y * sinlan,
			pos.y * coslan + pos.x * sinlan,
			pos.z
		);

		return pos;
	}
};

//Based on pykep.propagate_lagrangian
export function propagate(r0, v0, tof, mu) {
	if (tof === 0) {
		return {
			r: r0,
			v: v0
		};
	}
	if (tof < 0) {
		let res = propagate(r0, Vector3.mult(v0, -1), -tof, mu);

		return {
			r: res.r, 
			v: Vector3.mult(res.v, -1)
		};
	}

	let r = r0.norm;
	let v = v0.norm;

	let rvdot = Vector3.dot(r0, v0);

	let energy = (v * v / 2 - mu / r);
	let a = -mu / 2 / energy;

	let dM = Math.sqrt(Math.abs(mu / (a * a * a))) * tof;

	let C1 = rvdot / Math.sqrt(Math.abs(mu * a));
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

		let r2 = a + (r - a) * Math.cos(dE) + rvdot * Math.sqrt(a / mu) * Math.sin(dE);

		F = 1 - a / r * (1 - Math.cos(dE));
		G = a * rvdot / mu * (1 - Math.cos(dE)) + r * Math.sqrt(a / mu) * Math.sin(dE);
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
		// Set up:
		// f(dH)  = dot(r0, v0) * (cosh(dH) - 1) / sqrt(-mu * a) + (1 - r / a)sinh(dH) - dH - dM
		// f'(dH) = dot(r0, v0) * sinh(dH) / sqrt(-mu * a)       + (1 - r / a)cosh(dH) - 1

		// Newton's method is too unstable for hyperbolic orbits so I will do something else here
		// C2 > 0 since the orbit is hyperbolic, dM > 0 since tof > 0 in this part of the function
		// C1 can have any sign, but the absolute value of C1 is strictly less than C2

		let alpha = 0.5 * Math.log((C2 + C1) / (C2 - C1));
		let R = Math.sqrt(C2 * C2 - C1 * C1);

		// This part originally had 5 cases corresponding to different values of C1 and C2 when I was writing the code
		// but the above conditions removed all but this one

		let dH = -alpha;

		for (let i = 0; i < 15; i++) {
			dH = Math.asinh((C1 + dM + dH) / R) - alpha;
		}

		let r2 = a + (r - a) * Math.cosh(dH) + rvdot * Math.sqrt(-a / mu) * Math.sinh(dH);

		F = 1 - a / r * (1 - Math.cosh(dH));
		G = a * rvdot / mu * (1 - Math.cosh(dH)) + r * Math.sqrt(-a / mu) * Math.sinh(dH);
		Ft = -Math.sqrt(-mu * a) / (r2 * r) * Math.sinh(dH);
		Gt = 1 - a / r2 * (1 - Math.cosh(dH));
	}

	return {
		r: Vector3.add(Vector3.mult(r0, F),  Vector3.mult(v0, G)),
		v: Vector3.add(Vector3.mult(r0, Ft), Vector3.mult(v0, Gt))
	};
}
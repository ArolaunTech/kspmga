import { System } from './system.js';
import { Vector3 } from './vector.js';
import { realRootsQuartic } from './polynomial.js';

//Compute optimal MGA trajectories
export class MGAFinder {
	constructor(system) {
		this.system = system;
	}

	Fa(planet1, planet2, index, t1, t2, relvel) {
		let state1 = this.system.getState(planet1, t1);
		let state2 = this.system.getState(planet2, t2);

		let d1 = state1.r.norm;
		let d2 = state2.r.norm;

		let rdot = Vector3.dot(state1.r, state2.r);

		let rdiff = Vector3.sub(state1.r, state2.r);
		let c = rdiff.norm;

		let uc = Vector3.normalize(rdiff);
		let ui = Vector3.normalize(state1.r);

		let Kh = this.system.rootmu * c / (d1 * d2 + rdot);

		if (Kh === 0) {
			return 0;
		}

		let P = -2 * Vector3.dot(ui, state1.v);
		let Q = Vector3.dot(state1.v, state1.v) - relvel * relvel + 2 * Kh * Vector3.dot(ui, uc);
		let R = -2 * Kh * Vector3.dot(uc, state1.v);
		let S = Kh * Kh;

		if (!Math.isFinite(P) || !Math.isFinite(Q) || !Math.isFinite(R) || !Math.isFinite(S)) {
			return Infinity;
		}

		let vps = realRootsQuartic(P, Q, R, S);

		if (index >= vps.length) {
			return Infinity;
		}

		let vp = vps[index];
		let vel = Vector3.add(Vector3.mult(ui, vp), Vector3.mult(uc, Kh / vp));

		let h = Vector3.cross(state1.r, vel).norm;
		let p = h * h / this.system.rootmu;
		let sma = 1 / (2 / d1 - Vector3.dot(vel, vel) / this.system.rootmu);

		if (sma < 0) {
			return Infinity;
		}

		let cosdTheta = rdot / d1 / d2;
		let sindTheta = Math.sqrt(1 - cosdTheta * cosdTheta);

		if (vp < 0) {
			sindTheta *= -1;
		}

		let cosEa = 1 - (d1 * d2 - rdot) / (p * sma);
		let Ea = Math.acos(Ea);

		if (vp < 0) {
			Ea *= -1;
		}
		if ((d1 * d2 - rdot) > (d1 + d2) * p) {
			Ea *= -1;
		}

		let secsperradian = Math.sqrt(sma * sma * sma / this.system.rootmu);

		let ta = secsperradian * (Ea - Math.sin(Ea)) + d1 * d2 / Math.sqrt(this.system.rootmu * p) * sindTheta;

		return t2 - t1 - ta;
	}
};